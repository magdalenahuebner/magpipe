import pandas as pd
import numpy as np
import logging
import os
import pickle

from tqdm import tqdm
tqdm.pandas()

from scipy.stats import norm
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit

from src.stats_functions import z_score, get_pvalue, exp_func, random_fc
from src import visualisation as viz


def estimate_missing_fcs(dpoa, na_error):
    """
    Updates the fold change (FC) for phosphosites with specified NA error cases due to low signal 
    intensity either across control samples or within treatment samples. It uses Z-scores calculated
    from the signal distribution to estimate and adjust fold changes for these cases.

    Parameters:
    - dpoa (pd.DataFrame): DataFrame with initial DPOA results, including NA error cases, p-values, and fold change information.
    - na_error (str): Specifies the type of NA error case to update; can be 'control_missing' for phosphosites unquantified in 
                      control samples, or 'condition_missing' for those unquantified in a specific treatment.

    Returns:
    - pd.DataFrame: The updated dpoa DataFrame with fold changes adjusted for the specified NA error cases.
    """
    # Identify and extract rows in dpoa that match the specified na_error condition
    error_filter = dpoa['error'] == na_error
    if na_error == 'control_missing':
        error_filter &= (dpoa['n_cnd'] >= 2)  # Require at least 2 quantified conditions
    elif na_error == 'condition_missing':
        error_filter &= (dpoa['n_cnd'] == 0)  # No quantified conditions

    dpoa_error = dpoa.loc[error_filter].reset_index(drop=True)

    # Extract necessary columns for FC calculation
    dist_mean = dpoa_error['p_mean']  # Mean signal intensity distribution
    dist_std = dpoa_error['p_sd']  # Standard deviation of the signal intensity distribution
    fcs = dpoa_error['meansig_cnd' if na_error == 'control_missing' else 'meansig_ctr']  # The mean signal for condition or control
    ns = dpoa_error['n_cnd' if na_error == 'control_missing' else 'n_ctr']  # The sample size for condition or control

    # Calculate Z-scores to estimate how far each FC is from the mean of its distribution
    z_scores = [z_score(fcs.iloc[idx], dist_mean.iloc[idx], ns.iloc[idx], dist_std.iloc[idx]) for idx in dpoa_error.index]

    # Adjust FC calculation based on the na_error type
    if na_error == 'control_missing':
        fcs = [2 * norm.cdf(z) - 1 for z in z_scores]  # For missing control, adjust FC upwards
    elif na_error == 'condition_missing':
        fcs = [1 - 2 * norm.cdf(z) for z in z_scores]  # For missing condition, adjust FC downwards
    
    # Update the dpoa DataFrame with the new fold changes
    dpoa.loc[error_filter, 'fold_change'] = fcs
    logging.info(f"Updated {len(fcs)} fold changes where '{na_error}'.")

    return dpoa


def get_random_fcs(data, output_dir, force=False):
    """
    Generates or loads pre-generated random fold changes for phosphosites in the dataset. 
    If the output file exists and force is False, the existing random fold changes are loaded. 
    Otherwise, new random fold changes are generated, potentially overwriting existing files.

    Parameters:
    - data (pd.DataFrame): The dataset where rows represent phosphosites and columns represent samples.
    - output_dir (str): The directory to save or load the random fold changes file.
    - force (bool): If True, the function will generate and save new random fold changes regardless 
      of the existence of previous output files. Defaults to False.

    Returns:
    - pd.Series: A series containing the random fold changes for each phosphosite.
    """
    output_path = os.path.join(output_dir, 'error_model_fcs.pkl')

    if not force and os.path.exists(output_path):
        logging.info("Error model file exists. Loading from file.")
        with open(output_path, 'rb') as f:
            fcs = pickle.load(f)
    else:
        logging.info("Generating random fold changes.")
        fcs = data.progress_apply(random_fc, axis=1)
        with open(output_path, 'wb') as f:
            pickle.dump(fcs, f)
        logging.info(f"Random fold changes saved to {output_path}.")

    return fcs


def fit_em(data, output_dir, force=False, save_path=None):
    """
    Fits an error model to the dataset based on variability in fold changes across signal intensity bins. 
    The error model consists of a linear regression model fitted to the mean fold changes and an exponential 
    function fitted to the standard deviation of fold changes. The fold changes used for fitting the error model 
    are saved to or loaded from the specified directory.

    Parameters:
    - data (pd.DataFrame): The dataset containing signal intensity values for phosphosites (rows) 
      across samples (columns).
    - output_dir (str): The directory to save or load the fold changes used for fitting the error model.
    - force (bool): Forces the generation of random fold changes for phosphosites in the dataset if set to True. Default is False.
    - save_path (str, optional): Path to save the plot. If None, the plot is displayed directly.

    Returns:
    - reg_mean (LinearRegression object): The linear regression model fitted to mean fold changes.
    - popt (tuple): Parameters of the exponential function fitted to the standard deviation of fold changes.
    """
    # Preparing data for error model fitting (calculate means and random fold changes)
    logging.info("Preparing data for error model fitting.")
    signal_means = data.mean(axis=1)
    fcs = get_random_fcs(data, output_dir, force=force)  # Limiting step of the function

    # Invert approximately half the fold changes to make values symmetric about zero and make sure the values are sorted
    chosen_idx = fcs.sample(frac=0.5, random_state=1).index
    fcs[chosen_idx] = fcs[chosen_idx] * -1
    signal_means = signal_means.sort_values()
    fcs = fcs.reindex(signal_means.index.to_list())

    # Segment data into bins and calculate bin metrics for error model fitting
    bins = 100
    bin_size = int(len(signal_means) / bins)    
    bin_metrics = {'signal_means': [], 'fc_means': [], 'fc_stds': []}
    logging.info(f"Calculating bin metrics using {bins} bins.")
    for i in range(bins):
        bin_slice = slice(i * bin_size, (i + 1) * bin_size - 1)
        bin_metrics['signal_means'].append(np.mean(signal_means[bin_slice]))
        bin_metrics['fc_means'].append(np.mean(fcs[bin_slice]))
        bin_metrics['fc_stds'].append(np.std(fcs[bin_slice]))

    # Fit a linear model to the mean fold changes
    try:
        reg_mean = LinearRegression().fit(np.array(bin_metrics['signal_means']).reshape(-1, 1), np.array(bin_metrics['fc_means']).reshape(-1, 1))
        popt, _ = curve_fit(exp_func, np.array(bin_metrics['signal_means']), np.array(bin_metrics['fc_stds']))
    except Exception as e:
        logging.error(f"Error model fitting failed: {e}")
        return None, None  
    
    logging.info("Error model fitting complete.")

    # Plotting the error model
    viz.plot_error_model(x=signal_means, y=fcs, fit_data=bin_metrics, mean_fit=reg_mean, std_fit=popt, save_path=save_path)
    
    return reg_mean, popt


def estimate_sid_scores(dpoa, em_mean, em_std, na_error=None):
    """
    Estimates signal intensity dependent (SID) confidence scores for phosphosites with
    specified NA error cases, using a fitted error model to account for signal variability.

    Parameters:
    - dpoa (pd.DataFrame): DataFrame containing initial DPOA results, including NA error cases, fold changes, and signal intensity information.
    - em_mean (LinearRegression model): Linear regression model fitted to mean fold changes across signal intensity bins.
    - em_std (tuple): Parameters of the exponential function fitted to standard deviation of fold changes across signal intensity bins.
    - na_error (str, optional): Specifies the type of NA error case to address ('control_missing' or 'condition_missing'). Defaults to None.

    Returns:
    - pd.DataFrame: Updated dpoa DataFrame with 'sid_score' column containing the computed confidence scores.
    """
    logging.info(f"Estimating SID scores for NA error: {na_error}")

    # Filter dpoa rows based on the specified na_error condition
    error_filter = dpoa['error'] == na_error
    if not na_error:
        error_filter = dpoa['error'].isna()  # All normal cases with NA values
    elif na_error == 'control_missing':
        error_filter &= (dpoa['n_cnd'] >= 2)  # Require at least 2 quantified conditions for control missing
    elif na_error == 'condition_missing':
        error_filter &= (dpoa['n_cnd'] == 0)  # No quantified conditions for condition missing

    # Extract relevant data for SID score calculation
    fc_list = dpoa.loc[error_filter, 'fold_change'].to_list()
    meansig_list = dpoa.loc[error_filter, 'p_mean' if na_error == 'condition_missing' else 'meansig_cnd'].to_list()

    # Predict distribution mean using the error model and calculate standard deviation
    dist_mean = em_mean.predict(np.array(meansig_list).reshape(-1, 1))[:,0]
    dist_std = np.array([exp_func(x, *em_std) for x in meansig_list])

    # Calculate Z-scores and SID scores
    z_scores = [z_score(fc_list[idx], dist_mean[idx], 1, dist_std[idx]) for idx in range(len(fc_list))]
    sid_scores = [get_pvalue(z) for z in z_scores]

    # Update dpoa DataFrame with calculated SID scores
    dpoa.loc[error_filter, 'sid_score'] = sid_scores

    logging.info(f"SID scores estimated and updated for {len(sid_scores)} phosphosites with '{na_error}' error.")

    return dpoa

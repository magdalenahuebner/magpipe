import pandas as pd
import logging

from scipy.stats import norm

from src.stats_functions import z_score


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

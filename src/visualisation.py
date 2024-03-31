import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt

from scipy.stats import norm


def plot_qc_samples(df, p_quant=0.6, save_path=None):
    """
    Plots the fractions of quantified phosphosites for each sample and either displays the plot or saves it to a file.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data with samples as columns.
    - p_quant (float): The quantification threshold, highlighted on the plot.
    - save_path (str, optional): Path to save the plot. If None, the plot is displayed directly.
    """
    # Calculate the fraction of quantified phosphosites for each sample
    df_quant = 1 - df.isnull().mean()
    
    # Create a bar plot
    plt.figure(figsize=(10, 6))
    df_quant.plot(kind='bar')
    plt.axhline(y=p_quant, color='r', linestyle='--')
    plt.title('Fraction of Quantified Phosphosites per Sample')
    plt.xlabel('Samples')
    plt.ylabel('Fraction Quantified')
    plt.xticks(rotation=90, fontsize=8)  # Rotate and resize sample names for better readability
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels

    # Save or show the plot based on save_path
    if save_path:
        plt.savefig(save_path)
        logging.info(f"Plot saved to {save_path}")
    else:
        plt.show()


def plot_normalisation_boxplots(df, df_norm, save_path=None):
    """
    Generates boxplots for each sample in the DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data. Rows represent features, and columns represent samples.
    - title (str): Title of the plot.
    """
    # Create boxplots for the original data
    plt.figure(figsize=(12, 6))
    df.boxplot(rot=90, fontsize=8)
    plt.title('Before Normalisation')
    plt.ylabel('Area under the Peak (AUP)')  # Rotate and resize sample names for better readability
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels

    # Create boxplots for the normalised data
    plt.figure(figsize=(12, 6))
    df_norm.boxplot(rot=90, fontsize=8)
    plt.title('After Normalisation')
    plt.ylabel('Area under the Peak (AUP)')  # Rotate and resize sample names for better readability
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels

    if save_path:
        plt.savefig(save_path)
        logging.info(f"QC plots saved to {save_path}")
    else:
        plt.show()


def plot_phosphosite_signal_dist(data, phosphosite):
    """
    Plots the signal distribution for a given phosphosite across all samples, 
    with a normal distribution curve based on calculated mean and standard deviation.
    
    Parameters:
    - data (pd.DataFrame): DataFrame containing signal intensities with samples as columns and phosphosites as rows.
    - phosphosite (str): Phosphosite identifier.
    """
    # Filter data for the specified phosphosite
    phosphosite_data = data.loc[phosphosite]

    # Calculate mean and standard deviation of the phosphosite's distribution
    dist_mean = np.nanmean(phosphosite_data)
    dist_std = np.nanstd(phosphosite_data)

    # Plot histogram of phosphosite signal values
    plt.figure(figsize=(10, 6))
    plt.hist(phosphosite_data, bins=30, alpha=0.5, color='lightskyblue', density=True, label='Signal Intensity Distribution')

    # Plot normal distribution curve
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, dist_mean, dist_std)
    plt.plot(x, p, 'k', linewidth=2, label='Normal Distribution Curve')
    
    # Highlight the mean signal intensity
    plt.axvline(dist_mean, color='red', linestyle='dashed', linewidth=2, label=f'Mean Signal Intensity ({phosphosite})')
    plt.axvline(dist_mean + dist_std, color='steelblue', linestyle='dashed', linewidth=2, label=f'Standard Deviation ({phosphosite})')

    title = f'Signal Distribution for {phosphosite}'
    plt.title(title)
    plt.xlabel('Signal Intensity (log2(AUP))')
    plt.ylabel('Density')
    plt.legend()
    plt.show()


def plot_distance_to_signal_dist(dpoa, phosphosite, perturbagen, na_error='control_missing'):
    """
    Plot the signal distribution for a specific phosphosite and its signal intensity for a given condition.

    Parameters:
    - dpoa (pd.DataFrame): DataFrame containing DPOA results including 'phosphosite', 'perturbagen', 'p_mean', 'p_sd', 'meansig_cnd', and 'meansig_ctr'.
    - phosphosite (str): The phosphosite identifier.
    - perturbagen (str): The condition name.
    - signal_type (str): Type of signal to plot ('control' or 'condition').
    """
    # Extract distribution parameters
    p_mean = dpoa.loc[(dpoa['phosphosite'] == phosphosite) & (dpoa['perturbagen'] == perturbagen), 'p_mean'].values[0]
    p_sd = dpoa.loc[(dpoa['phosphosite'] == phosphosite) & (dpoa['perturbagen'] == perturbagen), 'p_sd'].values[0]

    # Define the signal intensity based on the NA error case
    if na_error == 'control_missing':
        signal = 'condition'
        signal_intensity = dpoa.loc[(dpoa['phosphosite'] == phosphosite) & (dpoa['perturbagen'] == perturbagen), 'meansig_cnd'].values[0]
    elif na_error == 'condition_missing':
        signal = 'control'
        signal_intensity = dpoa.loc[(dpoa['phosphosite'] == phosphosite) & (dpoa['perturbagen'] == perturbagen), 'meansig_ctr'].values[0]

    # Plot normal distribution curve
    x = np.linspace(p_mean - 4*p_sd, p_mean + 4*p_sd, 1000)
    p = norm.pdf(x, p_mean, p_sd)
    plt.figure(figsize=(10, 6))
    plt.plot(x, p, 'k', linewidth=2, label='Signal Distribution')

    # Highlight the mean signal intensity and the specific signal intensity of the condition of interest
    plt.axvline(p_mean, color='red', linestyle='dashed', linewidth=2, label=f'Mean Signal Intensity ({phosphosite})')
    plt.axvline(x=signal_intensity, color='steelblue', linestyle='dashed', linewidth=2, label=f'Mean Signal Intensity ({signal})')

    title = f'Signal Distribution for {phosphosite}, {perturbagen} ({signal}: {signal_intensity})'
    plt.title(title)
    plt.xlabel('Signal Intensity (log2(AUP))')
    plt.ylabel('Density')
    plt.legend()
    plt.show()

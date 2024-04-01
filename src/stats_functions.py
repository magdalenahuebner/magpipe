import pandas as pd
import numpy as np
import random

from scipy.stats import norm


def median_scaling(df, target_median=None):
    """
    Scales samples in the DataFrame to have the same median.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data to be scaled. Rows represent features, and columns represent samples.
    - target_median (float, optional): The target median value to scale each sample to. If None, uses the median of the medians across all samples.

    Returns:
    - pd.DataFrame: A DataFrame with samples scaled to the same median.
    """
    # Calculate the median of each sample
    sample_medians = df.median()
    
    # If no target_median is provided, calculate the median of the medians
    if target_median is None:
        target_median = sample_medians.median()

    # Scale each sample to the target median
    scaled_df = df.apply(lambda x: (x - x.median()) + target_median, axis=0)
    
    return scaled_df


def z_score(fc, dist_mean, n, dist_std):
    """
    Calculates the Z-score.

    Parameters:
    - fc: Fold change value.
    - dist_mean: Mean of the distribution.
    - n: Sample size.
    - dist_std: Standard deviation of the distribution.

    Returns:
    - Z-score.
    """
    return (fc - dist_mean) * np.sqrt(n) / dist_std


def get_pvalue(z):
    """
    Converts a Z-score to a corresponding p-value, reflecting the probability of observing such a value under the normal distribution.

    Parameters:
    - z: The Z-score to be converted.

    Returns:
    - float: The p-value corresponding to the Z-score.
    """
    return norm.cdf(z) if z < 0 else 1.0 - norm.cdf(z)
    

def exp_func(x, a, b, c):
    """
    Exponential function used for curve fitting, specifically fitting the standard deviation of fold changes.

    Parameters:
    - x: The independent variable, typically mean signal intensity.
    - a, b, c: Parameters of the exponential function to be optimized during curve fitting.

    Returns:
    - The value of the exponential function for a given x.
    """
    return a * b ** x + c


def random_fc(row):
    """
    Calculates an artificial fold change for a row by randomly dividing the measurements into two groups, then calculating the difference in their means.

    Parameters:
    - row (pd.Series): A row from the DataFrame, representing signal intensities for a single feature across different samples.

    Returns:
    - float: An artificial fold change calculated from the row data.
    """
    fcs = list()
    for _ in range(100):
        random.shuffle(row)
        g1 = row[:int(len(row) / 2)]
        g2 = row[int(len(row) / 2):]
        fcs.append(abs(np.mean(g1) - np.mean(g2)))
    fc = np.nanmean(fcs)
    return fc

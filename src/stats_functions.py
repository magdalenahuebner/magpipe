import pandas as pd
import numpy as np


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

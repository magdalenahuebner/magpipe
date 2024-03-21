import pandas as pd
import logging

import src.visualisation as viz
import src.stats_functions as stats


def filter_lq_samples(df, metadata, p_quant=0.6, save_path=None):
    """
    Removes low quality samples in the dataset based on a quantification threshold. Samples with less than the specified
    percentage of quantified phosphosites are removed.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data with columns as samples to be filtered.
    - metadata (pd.DataFrame): DataFrame containing the metadata with 'sampleID' used for renaming columns in data.
    - p_quant (float): The minimum percentage of quantified phosphosites per sample; defaults to 0.6 (60%).
    - save_path (str, optional): If provided, the quantification plot will be saved to this path.

    Returns:
    - tuple: A tuple containing two elements:
             1. The filtered DataFrame.
             2. The updated metadata DataFrame, corresponding to the filtered data.
    """
    assert 'sampleID' in metadata.columns, "Metadata must include 'sampleID' column."
   
    # Calculate and visualize the quantified phosphosites for each column 
    df_quant = 1 - df.isnull().mean()
    viz.plot_qc_samples(df, p_quant=p_quant, save_path=save_path)

    # Filter columns based on quantification threshold
    valid_cols = df_quant[df_quant >= p_quant].index
    filtered_df = df[valid_cols]
    filtered_mdat = metadata[metadata['sampleID'].isin(valid_cols)].reset_index(drop=True)    
    num_removed_samples = len(df.columns) - len(filtered_df.columns)
    logging.info(f"Removed {num_removed_samples} low-quality samples. Retained {len(filtered_df.columns)} samples.")
    
    return filtered_df, filtered_mdat


def filter_lq_phosphosites(df, metadata, condition_col, n_rep=2, n_cond=1, p_nas=0.5):
    """
    Removes phosphosites with poor quantification across all samples based on specified conditions.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the phosphoproteomics data with phosphosites as rows and samples as columns.
    - metadata (pd.DataFrame): DataFrame containing metadata.
    - condition_col (str): Column name in metadata DataFrame that holds the condition to group by.
    - n_rep (int): Minimum number of replicates required per phosphosite in at least n_cond conditions.
    - n_cond (int): Minimum number of conditions where the phosphosite should be present.
    - p_nas (float): Maximum fraction of NA (missing) values allowed per phosphosite across all samples.

    Returns:
    - tuple: A tuple containing two elements:
             1. filtered_data: The DataFrame after removing low-quality phosphosites.
             2. removed_psts: DataFrame of phosphosites removed due to poor quantification.
    """
    assert condition_col in metadata.columns, f"Metadata does not contain '{condition_col}' column."

    # Split data by groups specified in the condition column
    metadata = metadata[metadata['sampleID'].isin(df.columns)]  # align metadata with df
    grouped_data = {grp: df.loc[:, metadata.loc[metadata[condition_col] == grp, 'sampleID']] for grp in metadata[condition_col].unique()}

    # Apply the conditions across each group
    test_results = pd.DataFrame(index=df.index, columns=grouped_data.keys())
    for grp, grp_data in grouped_data.items():
        test_results[grp] = grp_data.apply(lambda x: x.count() >= n_rep, axis=1)

    # Determine phosphosites meeting the criteria
    valid_psts = (test_results.sum(axis=1) >= n_cond) | (df.isnull().mean(axis=1) <= p_nas)
    filtered_data = df.loc[valid_psts]
    removed_psts = df.loc[~valid_psts]  # Phosphosites removed due to poor quantification
    logging.info(f"Removed {len(removed_psts)} low-quality phosphosites. Retained {len(filtered_data)} phosphosites out of {len(df)} total phosphosites.")
    
    return filtered_data, removed_psts


def normalise_samples(df, save_path=None):
    """
    Normalizes phosphoproteomics data by scaling samples to have the same median and generates QC plots.

    Parameters:
    - df (pd.DataFrame): DataFrame containing phosphoproteomics data to be normalized.
    - save_path (str, optional): Path to save the QC plot. If None, plots are displayed directly.

    Returns:
    - pd.DataFrame: Normalized data frame.
    """
    # Applying median scaling normalization
    df_norm = stats.median_scaling(df)
    
    # Plotting QC after normalization
    viz.plot_normalisation_boxplots(df, df_norm, save_path=save_path)
    
    return df_norm

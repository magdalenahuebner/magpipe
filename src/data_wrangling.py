import numpy as np
import pandas as pd
import logging

import src.visualisation as viz


def rename_samples_based_on_metadata(df, metadata):
    """
    Renames the columns of a DataFrame based on the 'sampleID' column in the metadata DataFrame.
    Raises an error if the 'sampleID' column is missing in the metadata or if the column names in data do not match the index of metadata.

    Parameters:
    - df: DataFrame containing the data with columns (sample names) to be renamed.
    - metadata: DataFrame containing the metadata with 'sampleID' used for renaming columns in data.

    Returns:
    - The DataFrame with columns renamed.
    """
    assert 'sampleID' in metadata.columns, "Metadata must include 'sampleID' column."

    if metadata.columns[0] == 'sampleID':
        metadata.set_index(metadata.columns[0], drop=False, inplace=True)
    else:
        metadata.set_index(metadata.columns[0], inplace=True)

    # Check if data column names match metadata index
    if not set(df.columns).issubset(set(metadata.index)):
        raise ValueError("Column names of df do not match the index of metadata.")

    df = df.reindex(columns=metadata.index)
    df.rename(columns=metadata['sampleID'].to_dict(), inplace=True)
    metadata.reset_index(drop=True, inplace=True)
    return df, metadata


def filter_canonical_psts(df):
    """
    Filters the DataFrame to retain only rows with canonical phosphorylations (S, T, Y) in the specified column.

    Parameters:
    - df: DataFrame containing the data to be filtered, phosphosites (features) as row index formatted as ABC(S123).

    Returns:
    - Filtered DataFrame with only canonical phosphorylations.
    """
    initial_count = df.shape[0]
    pattern = r'\(S|\(T|\(Y'
    filtered_df = df[df.index.str.contains(pattern)]
    final_count = filtered_df.shape[0]

    removed_count = initial_count - final_count
    logging.info(f"Removed {removed_count} non-canonical phosphosite entries. Retained {final_count} entries.")

    return filtered_df


def merge_duplicates(df):
    """
    Merges duplicated entries in the DataFrame based on the index, taking the mean of duplicates.

    Parameters:
    - df: DataFrame with potential duplicated entries based on its index.

    Returns:
    - DataFrame with duplicates merged.
    """
    initial_count = df.shape[0]
    merged_df = df.groupby(df.index).mean()
    final_count = merged_df.shape[0]

    duplicates_count = initial_count - final_count
    logging.info(f"Merged {duplicates_count} duplicate entries. Total unique entries now: {final_count}.")
    
    return merged_df


def log2_transform_and_clean(df):
    """
    Applies a log2 transformation to the DataFrame and cleans the resulting -Inf values by replacing them with NaN.
    Checks for NaN values in the DataFrame and replaces them with 0, informing the user of this action.

    Parameters:
    - df: DataFrame containing the data to transform.

    Returns:
    - Transformed and cleaned DataFrame.
    """    
    if df.isnull().values.any():
        nan_count = df.isnull().sum().sum()
        logging.info(f"Detected {nan_count} NaN values in the DataFrame. Replacing NaN values with 0 before transformation.")
        df.fillna(0, inplace=True)
    
    transformed_df = np.log2(df)
    inf_count = np.isinf(transformed_df).sum().sum()
    # Calculate the percentage of -Inf values in relation to all values in the DataFrame
    inf_percentage = (inf_count / df.size) * 100
    transformed_df.replace(-np.inf, np.nan, inplace=True)
    
    logging.info(f"Replaced {inf_count} -Inf values with NaN after log2 transformation, which is approximately {inf_percentage:.2f}% of all values.")
    
    return transformed_df


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

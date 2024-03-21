import numpy as np
import pandas as pd
import logging


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

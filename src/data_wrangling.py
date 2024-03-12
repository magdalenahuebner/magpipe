import numpy as np
import pandas as pd


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
    if 'sampleID' not in metadata.columns:
        raise ValueError("Metadata does not contain 'sampleID' column.")

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
    return df


def filter_canonical_psts(df):
    """
    Filters the DataFrame to retain only rows with canonical phosphorylations (S, T, Y) in the specified column.

    Parameters:
    - df: DataFrame containing the data to be filtered, phosphosites (features) as row index formatted as ABC(S123).

    Returns:
    - Filtered DataFrame with only canonical phosphorylations.
    """
    pattern = r'\(S|\(T|\(Y'
    filtered_df = df[df.index.str.contains(pattern)]
    return filtered_df


def merge_duplicates(df):
    """
    Merges duplicated entries in the DataFrame based on the index, taking the mean of duplicates.

    Parameters:
    - df: DataFrame with potential duplicated entries based on its index.

    Returns:
    - DataFrame with duplicates merged.
    """
    return df.groupby(df.index).mean()


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
        print("NaN values detected in the DataFrame. Replacing NaN values with 0 before transformation.")
        df.fillna(0, inplace=True)

    transformed_df = np.log2(df).replace(-np.inf, np.nan)
    return transformed_df

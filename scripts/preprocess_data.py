# preprocess_data.py
import sys
import pandas as pd
from src.data_wrangling import *

def preprocess_data(aup_data_path, metadata_path, output_path):
    """
    Preprocesses LC/MS-MS phosphoproteomics experiment data by renaming columns according to metadata,
    filtering for canonical phosphorylations, removing duplicates, and applying log2 transformation.

    Parameters:
    - aup_data_path: Path to the CSV file containing Area Under the Peak (AUP) data from phosphoproteomics experiments.
                     Rows represent phosphosites (or features), and columns represent samples.
    - metadata_path: Path to the CSV file containing metadata for the samples. The first column must match the sample
                     names used in the AUP data file's columns.
    - output_path:   Path where the preprocessed AUP data CSV file will be saved.

    Returns:
    - None: The function saves the preprocessed data to a file specified by output_path.
    """
    # Load the AUP data and metadata
    aup = pd.read_csv(aup_data_path, sep='\t', index_col=0)
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=None)

    # Ensure the column names of AUP data match the first column of metadata for successful preprocessing
    if not set(aup.columns).issubset(set(metadata.iloc[:, 0])):
        raise ValueError("Column names of AUP data do not match the first column of metadata.")

    # Step 1: Rename columns based on 'sampleID' in metadata if 'sampleID' column is present
    aup = rename_samples_based_on_metadata(aup, metadata)

    # Step 2: Filter for canonical phosphorylations
    aup = filter_canonical_psts(aup)

    # Step 3: Remove duplicates
    aup = merge_duplicates(aup)

    # Step 4: Log2 transform and clean data
    aup = log2_transform_and_clean(aup)

    # Save the preprocessed AUP data
    aup.to_csv(output_path, sep='\t')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python preprocess_data.py <aup_data_path> <metadata_path> <output_path>")
        sys.exit(1)
    
    aup_data_path, metadata_path, output_path = sys.argv[1:]
    preprocess_data(aup_data_path, metadata_path, output_path)

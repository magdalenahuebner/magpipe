# perform_qc_and_normalise.py
import sys
import pandas as pd
import logging
import os

import src.magpipe.qc_and_normalisation as qcn

# Configure logging to display messages of level INFO and higher
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def perform_qc_and_normalise(aup_data_path, metadata_path, output_path, plot_dir):
    """
    Perform quality control (QC) and normalisation on LC/MS-MS phosphoproteomics experiment data.
    This involves excluding low-quality samples and phosphosites, followed by normalisation of the data.

    Parameters:
    - aup_data_path: Path to the TSV file containing Area Under the Peak (AUP) data from phosphoproteomics experiments.
                     Rows represent phosphosites (or features), and columns represent samples.
    - metadata_path: Path to the TSV file containing metadata for the samples. The first column must match the sample
                     names used in the AUP data file's columns.
    - output_path:   Path where the preprocessed AUP data TSV file will be saved.
    - plot_dir:      Directory where all generated plots will be saved

    Returns:
    - None: The function saves the preprocessed data to a file specified by output_path.
    """
    # Load the AUP data
    logging.info("Loading AUP data.")
    aup = pd.read_csv(aup_data_path, sep='\t', index_col=0)

    # Load metadata for renaming
    logging.info("Loading metadata... .")
    metadata = pd.read_csv(metadata_path, sep='\t')

    # Step 1: QC - Exclude low-quality samples
    logging.info("Performing QC: Excluding low-quality samples based on quantification threshold.")
    aup, metadata = qcn.filter_lq_samples(aup, metadata, p_quant=0.6, save_path=os.path.join(plot_dir, 'qc_plot_samples.png'))

    # Step 2: QC - Exclude low-quality phosphosites
    logging.info("Performing QC: Excluding low-quality phosphosites.")
    aup, _ = qcn.filter_lq_phosphosites(aup, metadata, condition_col='perturbagen', n_rep=2, n_cond=1, p_nas=0.5)

    # Step 3: Normalisation
    logging.info("Normalising sample peak areas to have the same median.")
    aup = qcn.normalise_samples(aup, save_path=os.path.join(plot_dir, 'normalisation_plot.png'))

    # Save the preprocessed AUP data
    logging.info(f"Saving preprocessed AUP data to {output_path}.")
    aup.to_csv(output_path, sep='\t')

if __name__ == "__main__":
    if len(sys.argv) != 5:
        logging.error("Usage: python preprocess_data.py <aup_data_path> <metadata_path> <output_path> <plot_dir>")
        sys.exit(1)
    
    aup_data_path, metadata_path, output_path, plot_dir = sys.argv[1:]
    perform_qc_and_normalise(aup_data_path, metadata_path, output_path, plot_dir)

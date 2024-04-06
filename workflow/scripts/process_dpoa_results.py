# process_dpoa_results.py
import sys
import pandas as pd
import logging
import os

import src.dpoa_processing as dpop

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_dpoa_results(aup_data_path, metadata_path, dpoa_path, results_path, output_dir, plot_dir):
    """
    Processes the results from Differential Phosphosite Occupancy Analysis (DPOA) by
    estimating missing fold changes for phosphosites with insufficient data and calculating
    Signal Intensity-Dependent Confidence Scores (SID-Scores).

    Parameters:
    - aup_data_path: Path to the TSV file containing the signal intensities.
    - metadata_path: Path to the TSV file containing metadata for the samples. The first column must match the sample
                     names used in the AUP data file's columns.
    - dpoa_path: Path to the TSV file containing DPOA results.
    - results_path: Path where the processed DPOA results will be saved.
    - output_dir: Directory to save or load pre-calculated or intermediate results.
    - plot_dir: Directory where all generated plots will be saved

    Returns:
    - None: Saves the processed DPOA results to the specified output path.
    """
    logging.info("Loading signal intensities data.")
    log2aup = pd.read_csv(aup_data_path, sep='\t', index_col=0)

    logging.info("Loading DPOA results.")
    dpoa = pd.read_csv(dpoa_path, sep='\t')

    logging.info("Loading metadata.")
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)

    # Step 1: Estimate missing fold changes
    # A: Cases where measurements are absent (mostly NA values) across control samples and all conditions have already been removed
    # B: Set cases where measurements are absent across control samples and a specific condition X to 0
    mask = (dpoa['error'] == 'control_missing') & (dpoa['n_cnd'] < 2)
    dpoa.loc[mask, 'fold_change'] = 0
    # C: Estimate missing fold changes for 'control_missing' cases
    dpoa = dpop.estimate_missing_fcs(dpoa, na_error='control_missing')
    # D: Estimate missing fold changes for 'condition_missing' cases
    dpoa = dpop.estimate_missing_fcs(dpoa, na_error='condition_missing')

    # Step 2: Calculate Signal Intensity-Dependent Confidence Scores (SID-Scores)
    # Fitting fold change error model to control samples
    controls = metadata.loc[metadata['perturbagen'] == 'control', 'sampleID']  # Identify control samples
    em_mean, em_std = dpop.fit_em(log2aup.loc[:, controls], output_dir=output_dir, force=False, save_path=os.path.join(plot_dir, 'fc_error_model_plot.png'))
    # Estimated SID-scores
    dpoa = dpop.estimate_sid_scores(dpoa, em_mean, em_std)
    # A: Cases where measurements are absent (mostly NA values) across control samples and all conditions have already been removed
    # B: Set cases where measurements are absent across control samples and a specific condition X to 1
    dpoa.loc[mask, 'sid_score'] = 1
    # C: Estimate sid-scores for 'control_missing' cases
    dpoa = dpop.estimate_sid_scores(dpoa, em_mean, em_std, na_error='control_missing')
    # D: Estimate sid-scores for 'condition_missing' cases
    dpoa = dpop.estimate_sid_scores(dpoa, em_mean, em_std, na_error='condition_missing')

    # Save the processed results
    logging.info(f"Saving processed DPOA results to {results_path}.")
    dpoa.to_csv(results_path, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        logging.error("Usage: python process_dpoa_results.py <aup_data_path> <metadata_path> <dpoa_path> <results_path> <output_dir> <plot_dir>")
        sys.exit(1)

    aup_data_path, metadata_path, dpoa_path, results_path, output_dir, plot_dir = sys.argv[1:]
    process_dpoa_results(aup_data_path, metadata_path, dpoa_path, results_path, output_dir, plot_dir)

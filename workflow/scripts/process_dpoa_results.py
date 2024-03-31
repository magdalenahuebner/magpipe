# process_dpoa_results.py
import sys
import pandas as pd
import logging

import src.dpoa_processing as dpop

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_dpoa_results(aup_data_path, dpoa_path, output_path):
    """
    Processes the results from Differential Phosphosite Occupancy Analysis (DPOA) by
    estimating missing fold changes for phosphosites with insufficient data.

    Parameters:
    - aup_data_path: Path to the CSV file containing the signal intensities.
    - dpoa_path: Path to the CSV file containing DPOA results.
    - output_path: Path where the processed DPOA results will be saved.

    Returns:
    - None: Saves the processed DPOA results to the specified output path.
    """
    logging.info("Loading signal intensities data.")
    data = pd.read_csv(aup_data_path, sep='\t', index_col=0)

    logging.info("Loading DPOA results.")
    dpoa = pd.read_csv(dpoa_path, sep='\t')

    # Step 1: Estimate missing fold changes
    # A: Cases where measurements are absent (mostly NA values) across control samples and all conditions have already been removed
    # B: Set cases where measurements are absent across control samples and a specific condition X to 0
    mask = (dpoa['error'] == 'control_missing') & (dpoa['n_cnd'] < 2)
    dpoa.loc[mask, 'fold_change'] = 0
    # C: Estimate missing fold changes for 'control_missing' cases
    dpoa = dpop.estimate_missing_fcs(dpoa, 'control_missing')
    # D: Estimate missing fold changes for 'condition_missing' cases
    dpoa = dpop.estimate_missing_fcs(dpoa, 'condition_missing')

    # Save the processed results
    logging.info(f"Saving processed DPOA results to {output_path}.")
    dpoa.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        logging.error("Usage: python process_dpoa_results.py <aup_data_path> <dpoa_path> <output_path>")
        sys.exit(1)

    aup_data_path, dpoa_path, output_path = sys.argv[1:]
    process_dpoa_results(aup_data_path, dpoa_path, output_path)

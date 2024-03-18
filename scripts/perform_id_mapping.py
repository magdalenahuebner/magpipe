import sys
import pandas as pd
import logging

import src.id_mapping as idm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def perform_id_mapping(aup_data_path, output_path):
    """
    Perform ID mapping on phosphoproteomics experiment data, updating protein IDs to HGNC symbols.

    Parameters:
    - aup_data_path: Path to the CSV file containing Area Under the Peak (AUP) data from phosphoproteomics experiments.
                     Rows represent phosphosites (or features), and columns represent samples.
    - output_path:   Path where the preprocessed AUP data CSV file will be saved.

    Returns:
    - None: The function saves the preprocessed data to a file specified by output_path.
    """
    # Load the AUP data
    logging.info("Loading AUP data.")
    aup = pd.read_csv(aup_data_path, sep='\t', index_col=0)

    # 1. Initial UniProt mapping
    logging.info("Performing initial UniProt ID mapping.")
    aup = idm.index_ids_to_hgnc_symbols(aup, database='uniprot_protmapper', phosphosites=True)

    # 2. Handle unmatched IDs with UniProt search
    logging.info("Handling unmatched IDs with UniProt search.")
    filtered_ids = [id_ for id_ in aup.index.to_list() if '_HUMAN' in id_]
    # Splitting the original DataFrame into two parts:
    aup_a = aup.loc[filtered_ids]  # Entries for updating
    aup_b = aup.drop(filtered_ids)  # Entries that don't need updating
    # Update the entries with the newly mapped HGNC symbols using UniProt search results.
    aup_a = idm.index_ids_to_hgnc_symbols(aup_a, database='uniprot_search', phosphosites=True)
    # Reassemble the DataFrame by concatenating the updated part with the unchanged part.
    aup = pd.concat([aup_a, aup_b])
    
    # 3. Ensure gene names correspond to currently approved HGNC symbols
    logging.info("Updating gene names to approved HGNC symbols.")
    aup = idm.index_ids_to_hgnc_symbols(aup, database='hgnc', phosphosites=True)

    # 4. Manually map remaining unmatched IDs
    logging.info("Processing manually mapped IDs.")
    # This assumes 'mapped_hgnc_symbols.tsv' is formatted as 'From\tTo'
    mapped_hgnc_symbols = pd.read_csv('resources/external/mapped_hgnc_symbols.tsv', sep='\t')
    mapping_dict = dict(zip(mapped_hgnc_symbols['From'], mapped_hgnc_symbols['To']))
    aup = idm.update_df_index_with_mappings(aup, mapping_dict, phosphosites=True)

    # Save the processed AUP data
    logging.info(f"Saving ID mapped AUP data to {output_path}.")
    aup.to_csv(output_path, sep='\t')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        logging.error("Usage: python perform_id_mapping.py <aup_data_path> <output_path>")
        sys.exit(1)
    
    aup_data_path, output_path = sys.argv[1:]
    perform_id_mapping(aup_data_path, output_path)

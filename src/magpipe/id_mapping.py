import pandas as pd
import logging
import re
import requests
import concurrent.futures

from UniProtMapper import ProtMapper


def extract_protein_ids(phosphosites, unique=True):
    """
    Extracts protein identifiers from a list of phosphosite strings. Each phosphosite string is expected
    to be formatted as 'PROTEIN(SITE)', where 'PROTEIN' is the protein identifier and 'SITE' is the phosphorylation site.

    Parameters:
    - phosphosites (list of str): A list of strings representing phosphosites.
    - unique (bool): If True (default), returns a list of unique protein identifiers. If False, returns all identifiers including duplicates.

    Returns:
    - list: A list of extracted protein identifiers. The list contains only unique identifiers if 'unique' is True.
    """
    protein_pattern = re.compile(r'([A-Za-z0-9_.\-]+)\(.*\)')
    protein_ids = [protein_pattern.match(phosphosite).group(1) for phosphosite in phosphosites if protein_pattern.match(phosphosite)]
    return list(dict.fromkeys(protein_ids)) if unique else protein_ids


def convert_phosphosite_ids(phosphosites, id_map):
    """
    Converts a list of phosphosites from old protein identifiers to new protein identifiers based on a provided mapping dictionary.
    Each phosphosite is represented as 'PROTEIN(SITE)' where 'PROTEIN' is the protein identifier and 'SITE' is the phosphorylation site.
    If a protein identifier in the phosphosites list does not exist in the id_map, it remains unchanged.

    Parameters:
    - phosphosites (list of str): A list of strings representing phosphosites.
    - id_map (dict): A dictionary where keys are old protein identifiers and values are the corresponding new protein identifiers.

    Returns:
    - list of str: A list of phosphosites with protein identifiers converted according to the id_map. If a protein identifiers is not found
      in the id_map, it remains unchanged in the output list.
    """
    protein_pattern = re.compile(r'([A-Za-z0-9_.\-]+)\((.*)\)')
    converted_phosphosites = []
    for phosphosite in phosphosites:
        match = protein_pattern.match(phosphosite)
        if match:
            protein_id, site = match.groups()
            new_id = id_map.get(protein_id, protein_id)
            converted_phosphosites.append(f'{new_id}({site})')
        else:
            converted_phosphosites.append(phosphosite)
    return converted_phosphosites


def make_request(url, params=None, headers=None, max_retries=3):
    """
    Makes a HTTP GET request with retry mechanism.
    """
    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()  # Raises HTTPError for bad responses
            return response.json()
        except requests.exceptions.HTTPError as e:
            logging.warning(f"HTTPError for URL {url}: {e}")
            retries += 1
        except requests.exceptions.RequestException as e:
            logging.error(f"RequestException for URL {url}: {e}")
            return None
    logging.error(f"Failed to fetch data after {max_retries} attempts for URL: {url}")
    return None


def fetch_approved_hgnc_symbol(protein_id):
    """
    Fetches the current approved HGNC (HUGO Gene Nomenclature Committee) symbol for a given gene identifier.
    The function searches the identifier as a current symbol, previous symbol, or alias symbol in HGNC.

    Parameters:
    - protein_id: A string representing a protein identifier (could be a symbol, previous symbol, alias symbol, etc.)

    Returns:
    - A dictionary containing the input identifier, match type (indicates if it's a current symbol, previous symbol, or alias symbol),
      and the current approved HGNC symbol. If no match is found, returns a dictionary indicating the identifier is unmatched.
    """
    api_url = "https://rest.genenames.org/search/"
    headers = {'Accept': 'application/json'}
    
    for symbol_type in ['symbol', 'prev_symbol', 'alias_symbol']:
        response_data = make_request(f'{api_url}{symbol_type}/{protein_id}', headers=headers)
        if response_data and response_data.get('response', {}).get('docs'):
            docs = response_data['response']['docs']
            return {'Input': protein_id, 'Match type': symbol_type, 'Approved symbol': docs[0].get('symbol')}
    return {'Input': protein_id, 'Match type': 'unmatched', 'Approved symbol': None}


def fetch_gene_name_from_uniprot(query):
    """
    Fetches the gene name for a given query from UniProt.
    
    This function queries the UniProt database to find the gene name associated with the given query. 
    The query can be a gene symbol, UniProt accession number, or any identifier recognized by UniProt.
    
    Parameters:
    - query (str): The query string used to search the UniProt database.
    
    Returns:
    - str or None: The gene name associated with the query if found; otherwise, None.
    """
    # Define the base URL for the UniProt search service
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # Perform the GET request
    response = requests.get(search_url, params={'query': query, 'format': 'json'})

    # Check if the response is successful
    if response.status_code == 200:
        results = response.json()
        # Access the first result if it exists
        if results.get('results'):
            first_result = results.get('results')[0]
            gene_name = first_result.get('genes', [{}])[0].get('geneName', {}).get('value', None)
        else:
            logging.error(f"\n{query} doesn't exist in UniProtKB database.")
            return None
    else:
        logging.error("\nFailed to retrieve data.")
        return None
    
    return gene_name


def get_uniprot_id_map(ids, from_id='UniProtKB_AC-ID', to_id='Gene_Name', tool='uniprot_protmapper', phosphosites=False, return_df=False):
    """
    Maps UniProt IDs to gene (protein) names using either the ProtMapper package or the UniProt search service.

    Parameters:
    - ids (list): A list of IDs to map.
    - from_id (str): The type of the input IDs (applicable only when tool is 'uniprot_protmapper').
    - to_id (str): The target ID type for mapping. Currently, only "Gene_Name" is supported for 'uniprot_search'.
    - tool (str): The mapping tool to use ('uniprot_protmapper' or 'uniprot_search').
    - phosphosites (bool): If True, preprocesses the IDs to extract protein IDs for phosphosites (applicable only when tool is 'protmapper').
    - return_df (bool): If True, returns a pandas DataFrame; otherwise, returns a dictionary.

    Returns:
    - pandas.DataFrame or dict: Mapped IDs in a DataFrame or dict format, depending on the value of return_df.
    """
    protein_ids = extract_protein_ids(ids, unique=True) if phosphosites else list(dict.fromkeys(ids))

    if tool == 'uniprot_protmapper':
        mapper = ProtMapper()
        result, failed = mapper.get(ids=protein_ids, from_db=from_id, to_db=to_id)

    elif tool == 'uniprot_search':
        gene_names = []
        for i, id in enumerate(protein_ids):
            gene_name = fetch_gene_name_from_uniprot(id)
            gene_names.append(gene_name)
            print(f"Processed {i+1}/{len(protein_ids)} identifiers.", end='\r')
        result = pd.DataFrame({'From': protein_ids, 'To': gene_names}).dropna()

    else:
        logging.error("Unsupported tool for UniProt ID mapping.")
    
    logging.info(f"\nMapped {len(result)} out of {len(protein_ids)} proteins from {from_id}s to {to_id}s.")
    if return_df:
        return result
    return dict(zip(result['From'], result['To']))


def get_hgnc_id_map(ids, phosphosites=False, return_df=False):
    """
    Fetches current approved HGNC symbols for a list of gene identifiers concurrently.
    
    Parameters:
    - ids (list): A list of gene identifiers.
    - phosphosites (bool): If True, preprocesses the IDs to extract protein IDs for phosphosites.
    - return_df (bool): If True, returns a pandas DataFrame; otherwise, returns a dictionary.

    Returns:
    - Either a DataFrame or a dictionary with mappings from input identifiers to approved HGNC symbols,
      depending on the value of return_df.
    """
    protein_ids = extract_protein_ids(ids, unique=True) if phosphosites else list(dict.fromkeys(ids))
    results = []

    logging.info(f"Starting to fetch symbols for {len(protein_ids)} identifiers...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        future_to_id = {executor.submit(fetch_approved_hgnc_symbol, pid): pid for pid in protein_ids}
        for i, future in enumerate(concurrent.futures.as_completed(future_to_id)):
            results.append(future.result())
            print(f"Processed {i+1}/{len(protein_ids)} identifiers", end='\r')

    df = pd.DataFrame(results)
    filtered_df = df[(df['Match type'] != 'symbol') & (df['Approved symbol'].notnull())]
    logging.info(f"\nMapped {len(filtered_df)} out of {len(protein_ids)} proteins to approved HGNC gene symbols.")
    if return_df:
        return df
    return filtered_df.set_index('Input')['Approved symbol'].to_dict() if not df.empty else {}


def update_df_index_with_mappings(df, mapping_dict, phosphosites=False):
    """
    Updates the DataFrame index with new identifiers based on a given mapping dictionary.
    It handles special processing for phosphosites if indicated.

    Parameters:
    - df (pd.DataFrame): The DataFrame whose index needs to be updated.
    - mapping_dict (dict): A dictionary with current identifiers as keys and new identifiers as values.
    - phosphosites (bool): Indicates if the identifiers are phosphosites, which requires special handling.

    Returns:
    - pd.DataFrame: The updated DataFrame with new identifiers in the index.
    """
    index = df.index.tolist()
    updated_index = convert_phosphosite_ids(index, mapping_dict) if phosphosites else [mapping_dict.get(id_, id_) for id_ in index]
    transformed_count = sum(1 for original, new in zip(index, updated_index) if original != new)
    id_type = 'phosphosites' if phosphosites else 'proteins'
    logging.info(f"Transformed {transformed_count} out of {len(index)} {id_type}.")
    updated_df = df.copy()
    updated_df.index = updated_index
    return updated_df


def index_ids_to_hgnc_symbols(df, database='hgnc', phosphosites=False):
    """
    Updates DataFrame index to HGNC symbols based on a specified database.

    Parameters:
    - df (pd.DataFrame): The input DataFrame with identifiers as its index.
    - database (str): The source database to use for mapping. Options include 'uniprot_protmapper', 'uniprot_search', and 'hgnc'.
    - phosphosites (bool): Indicates whether the identifiers are phosphosites, which requires special processing.

    Returns:
    - pd.DataFrame: A DataFrame identical to the input but with its index updated to HGNC symbols.
    """
    if database not in ['uniprot_protmapper', 'uniprot_search', 'hgnc']:
        raise ValueError("Invalid database specified. Choose from 'uniprot_protmapper', 'uniprot_search', 'hgnc'.")

    ids = df.index.tolist()
    mapping_dict = {}

    if database in ['uniprot_protmapper', 'uniprot_search']:
        mapping_dict = get_uniprot_id_map(ids, tool=database, phosphosites=phosphosites)
    elif database == 'hgnc':
        mapping_dict = get_hgnc_id_map(ids, phosphosites=phosphosites)

    return update_df_index_with_mappings(df, mapping_dict, phosphosites=phosphosites)

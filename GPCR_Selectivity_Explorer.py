import argparse
import pandas as pd
import json
import requests
import os
import logging
import shutil

from urllib.request import urlopen, urlretrieve
from Bio import Entrez
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from pymol import cmd, stored
from Bio.SeqUtils import seq1 as rename_AA
from typing import List, Tuple
from contextlib import contextmanager

# Set your email address (required by NCBI)
Entrez.email = "nour.aldin.kahlous@icm.uu.se"

# Configure logging to output to standard output with a specified format
logging.basicConfig(
    level=logging.INFO,  # Set to DEBUG for more granular logging
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]  # StreamHandler outputs to standard output
)


@contextmanager
def change_dir(destination):

    if os.path.exists(destination):
        shutil.rmtree(destination)  # Deletes the directory and its contents

    os.makedirs(destination)  # Creates the directory (including parent directories if needed)

    current_dir = os.getcwd()  # Get the current working directory
    try:
        os.chdir(destination)  # Change to the new directory
        yield  # Allows the code block to execute within the context
    finally:
        os.chdir(current_dir)  # Revert back to the original directory

def retrieve_json_file(protein_name: str) -> str:
    """
    Retrieves a JSON file containing residue information for a given protein from the GPCRdb API
    and saves it locally.

    Parameters
    ----------
    protein_name : str
        The name of the protein for which to retrieve the JSON file.

    Returns
    -------
    str
        The filename of the saved JSON file if successful, or None if the retrieval failed.
    """
    url = f"https://gpcrdb.org/services/residues/extended/{protein_name}/"
    headers = {"accept": "application/json"}

    logging.debug("Attempting to retrieve JSON file for protein: %s", protein_name)

    try:
        # Make the GET request
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()  # Raise an HTTPError for bad responses
    except requests.exceptions.RequestException as e:
        logging.error("Failed to retrieve JSON for %s. Error: %s", protein_name, e)
        return None

    # Save the JSON response to a file if the request was successful
    json_file = f"{protein_name}_table.json"
    try:
        with open(json_file, 'w') as f:
            json.dump(response.json(), f, indent=4)
        logging.debug("JSON file for %s saved as %s", protein_name, json_file)
        return json_file
    except IOError as e:
        logging.error("Failed to save JSON file for %s. Error: %s", protein_name, e)
        return None


def json_to_pandas2(json_file: str) -> pd.DataFrame:
    """
    Converts a JSON file of protein residue data into a Pandas DataFrame, extracting and transforming
    relevant fields to include specific residue numbering formats.

    Parameters
    ----------
    json_file : str
        The filename of the JSON file to convert.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame containing the transformed data, with additional columns for residue
        numbers and sequence information.
    """
    logging.info("Attempting to load JSON file: %s", json_file)

    try:
        # Open the JSON file and load the data
        with open(json_file, 'r') as f:
            data = json.load(f)
        logging.debug("JSON file %s successfully loaded", json_file)
    except (FileNotFoundError, IOError) as e:
        logging.error("Failed to load JSON file %s. Error: %s", json_file, e)
        return pd.DataFrame()  # Return empty DataFrame on failure

    try:
        # Convert the list of dictionaries to a Pandas DataFrame
        df = pd.DataFrame(data)
        logging.debug("Data successfully converted to DataFrame")

        # Split 'display_generic_number' to extract residue numbering information
        df[['residue_number', 'part_2']] = df['display_generic_number'].str.split('.', expand=True)
        df[['Sequence-based (BW)', 'Structure-based (GPCRdb)']] = df['part_2'].str.split('x', expand=True)
        df.drop(columns='part_2', inplace=True)

        # Create 'GPCRdb(A)' notation column
        df['GPCRdb(A)'] = df['residue_number'].astype(str) + 'x' + df['Structure-based (GPCRdb)'].astype(str)
        logging.debug("Residue number columns successfully extracted and transformed")

        # Add amino acid sequence column with file prefix
        prefix = "_".join(json_file.split('_')[:2])
        df[f"{prefix}_AA-seq"] = df['amino_acid'] + df['sequence_number'].astype(str)

        logging.info("Data transformation complete for file: %s", json_file)
        return df

    except KeyError as e:
        logging.error("Error in transforming data. Missing expected key: %s", e)
        return pd.DataFrame()  # Return empty DataFrame if transformation fails


def filter_by_gpcrdb(df: pd.DataFrame, value_list: List[str]) -> List[str]:
    """
    Filters a DataFrame by matching specified values in the 'protein_segment' column, 
    and returns a list of valid GPCRdb positions.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing protein data with a 'protein_segment' and 'GPCRdb(A)' column.
    value_list : List[str]
        A list of segment names to filter the DataFrame by.

    Returns
    -------
    List[str]
        A list of GPCRdb positions from the 'GPCRdb(A)' column that correspond to the specified 
        segments and exclude any non-valid entries.
    """
    logging.info("Filtering DataFrame by provided segment values")

    # Check if required columns are in the DataFrame
    if 'protein_segment' not in df.columns or 'GPCRdb(A)' not in df.columns:
        logging.error("Missing required columns 'protein_segment' or 'GPCRdb(A)' in DataFrame.")
        return []

    try:
        # Filter the DataFrame for rows where 'protein_segment' is in the value_list
        filtered_df = df[df['protein_segment'].isin(value_list)]
        logging.debug("DataFrame successfully filtered by 'protein_segment'")

        # Extract and clean up the GPCRdb positions
        positions = [
            item for item in filtered_df['GPCRdb(A)'] if item and item != 'NonexNone'
        ]
        logging.info("Filtered positions extracted successfully. Total valid positions: %d", len(positions))
        
        return positions

    except Exception as e:
        logging.error("Error occurred while filtering positions. Error: %s", e)
        return []


def fetch_species_positions(protein_names: List[str], positions: List[str]) -> Tuple[List[str], List[str]]:
    """
    Fetches protein data for a list of proteins, extracts their families, retrieves family alignment data,
    and returns a list of species and common GPCR positions.

    Parameters
    ----------
    protein_names : List[str]
        A list of protein names (e.g., ['5ht1a_human', 'adrb2_human']).
    positions : List[str]
        A list of GPCR positions to filter by.

    Returns
    -------
    Tuple[List[str], List[str]]
        A tuple containing:
        - full_species_list: Combined list of species (or indices) from the alignment data.
        - sorted_positions_list: List of common positions across the specified proteins.
    """
    full_species_list = []
    positions_list = None

    logging.info("Starting fetch for species and positions for proteins: %s", protein_names)

    for protein_name in protein_names:
        try:
            # Fetch the protein data
            protein_url = f"https://gpcrdb.org/services/protein/{protein_name}/"
            response = urlopen(protein_url)
            protein_data = json.loads(response.read().decode('utf-8'))
            logging.info("Protein data fetched successfully for: %s", protein_name)
            
            # Extract the family
            family = protein_data.get('family')
            if not family:
                logging.warning("Family data missing for protein: %s", protein_name)
                continue

            # Fetch the family alignment data
            alignment_url = f"https://gpcrdb.org/services/alignment/family_all/{family}"
            response = urlopen(alignment_url)
            alignment_data = json.loads(response.read().decode('utf-8'))
            logging.info("Alignment data fetched successfully for family: %s", family)

            # Convert the alignment data into a pandas DataFrame
            df = pd.DataFrame.from_dict(alignment_data, orient='index')
            species = list(df.index.values)
            full_species_list.extend(species)
            logging.debug("Species list extended with family alignment data for %s", protein_name)

        except Exception as e:
            logging.error("Error fetching data for %s: %s", protein_name, e)
            continue

        # Retrieve JSON file for the protein and filter positions
        json_file = retrieve_json_file(protein_name)
        if json_file:
            df = json_to_pandas2(json_file)
            plist = filter_by_gpcrdb(df, positions)
            logging.debug("Positions filtered for %s: %s", protein_name, plist)

            # Update positions list with the intersection of positions
            if positions_list is None:
                positions_list = plist
            else:
                positions_list = list(set(positions_list).intersection(plist))
                logging.info("Positions list updated by intersection for %s", protein_name)

    # Sort final positions list
    sorted_positions_list = sorted(positions_list) if positions_list else []
    logging.info("Final species and positions lists prepared")

    return full_species_list, sorted_positions_list


def fetch_gpcrdb_alignment(protein_list: List[str], gpcrdb_positions: List[str], batch_size: int = 10) -> pd.DataFrame:
    """
    Fetches the GPCRdb alignment data for a list of proteins in batches.

    Parameters
    ----------
    protein_list : List[str]
        List of protein names (e.g., ['adrb1_human', 'adrb2_human']).
    gpcrdb_positions : List[str]
        List of GPCRdb positions to retrieve (e.g., ['7x30', '7x31']).
    batch_size : int, optional
        Maximum number of proteins per request (default is 10).

    Returns
    -------
    pd.DataFrame
        Combined alignment data for all proteins, with one column for each GPCRdb position.
    """
    logging.info("Starting GPCRdb alignment data fetch for proteins in batches")
    full_df = pd.DataFrame()
    gpcrdb_positions_str = ",".join(gpcrdb_positions)

    # Process proteins in batches
    for i in range(0, len(protein_list), batch_size):
        batch = protein_list[i:i + batch_size]
        protein_str = ",".join(batch)
        url = f'https://gpcrdb.org/services/alignment/protein/{protein_str}/{gpcrdb_positions_str}/'
        
        try:
            # Fetch data from the GPCRdb API
            response = urlopen(url)
            alignment_data = json.loads(response.read().decode('utf-8'))
            logging.debug("Batch %d data fetched successfully for proteins: %s", i // batch_size + 1, batch)

            # Convert the alignment data into a DataFrame
            df = pd.DataFrame.from_dict(alignment_data, orient='index').reset_index()

            # Separate sequence into columns for each GPCRdb position
            split_sequences = df[0].str.split('', expand=True)
            split_sequences.drop(columns=[0, len(split_sequences.columns) - 1], inplace=True)
            logging.debug("Sequences split into individual positions for batch %d", i // batch_size + 1)

            # Combine protein names with the split sequences
            result_df = pd.concat([df['index'], split_sequences], axis=1)
            result_df.columns = ['GPCRdb(A)'] + gpcrdb_positions

            # Append batch data to the full DataFrame
            full_df = pd.concat([full_df, result_df], ignore_index=True)

        except Exception as e:
            logging.error("Error fetching data for batch %d: %s", i // batch_size + 1, e)
            continue

    logging.info("GPCRdb alignment data fetch complete")
    return full_df


def add_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds rows for the consensus amino acid sequence and consensus values at each position 
    using the frequency of each amino acid.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing sequence data, with the first column as protein names
        and the remaining columns as sequence positions.

    Returns
    -------
    pd.DataFrame
        DataFrame with additional rows for the consensus amino acid sequence and consensus values.
    """
    logging.info("Starting consensus calculation for alignment")

    try:
        # Extract sequences from DataFrame
        sequences = [
            Seq("".join(row[1:].values))  # Join amino acids, ignoring the first column (protein names)
            for _, row in df.iterrows()
        ]

        # Create a MultipleSeqAlignment object
        alignment = MultipleSeqAlignment([SeqRecord(seq, id=f"seq{i}") for i, seq in enumerate(sequences)])
        alignment_length = alignment.get_alignment_length()
        logging.debug("Alignment created with %d sequences and length %d", len(sequences), alignment_length)

        # Initialize amino acid count dictionary
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        counts = {aa: [0] * alignment_length for aa in alphabet}

        # Count amino acids at each position
        for i in range(alignment_length):
            column = alignment[:, i]
            for amino_acid in column:
                if amino_acid in counts:
                    counts[amino_acid][i] += 1

        # Determine consensus sequence
        consensus_list = []
        for i in range(alignment_length):
            position_counts = [(aa, counts[aa][i]) for aa in alphabet]
            max_count = max(count for _, count in position_counts)
            max_aas = [aa for aa, count in position_counts if count == max_count]
            consensus_list.append(max_aas[0] if len(max_aas) == 1 else 'X')

        consensus_row = pd.Series(["Seq consensus"] + consensus_list, index=df.columns)
        logging.info("Consensus sequence row generated")

        # Calculate consensus values
        consensus_values = []
        for col in df.columns[1:]:  # Skip protein names
            col_amino_acids = df[col].tolist()
            consensus_aa = consensus_row[col]
            consensus_count = col_amino_acids.count(consensus_aa)
            consensus_value = (consensus_count / len(col_amino_acids)) * 100
            consensus_values.append(f"{consensus_value:.1f}")
        
        consensus_values_row = pd.Series(["Seq consensus value"] + consensus_values, index=df.columns)
        logging.info("Consensus values row generated")

        # Concatenate the consensus rows to the original DataFrame
        df_with_consensus = pd.concat([df, consensus_row.to_frame().T, consensus_values_row.to_frame().T], ignore_index=True)
        logging.info("Consensus rows added to DataFrame")

        return df_with_consensus

    except Exception as e:
        logging.error("Error occurred while calculating consensus: %s", e)
        return df  # Return the original DataFrame if an error occurs



def process_alignment(df_with_consensus: pd.DataFrame, ref: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Processes the consensus DataFrame by transposing it, setting headers, and merging with reference residues.

    Parameters
    ----------
    df_with_consensus : pd.DataFrame
        Input DataFrame containing consensus sequence information.
    ref : str
        The reference protein name used to retrieve JSON file and align residue numbering.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        - df_ref_resides_table : DataFrame with reference residues and their GPCRdb numbering.
        - df_ref_consensus_with_residues_numbering : Merged DataFrame with consensus and reference residue numbering.
    """
    logging.info("Starting processing of alignment with reference: %s", ref)

    try:
        # Transpose the DataFrame
        df_transposed = df_with_consensus.T
        logging.debug("DataFrame transposed successfully")

        # Set the first row as the header and drop it from data
        df_transposed.columns = df_transposed.iloc[0]
        df_transposed = df_transposed.drop(df_transposed.index[0])
        logging.debug("Header set and first row dropped")

        # Rename index as 'GPCRdb(A)' and reset it
        df_transposed = df_transposed.rename_axis('GPCRdb(A)').reset_index()

        # Convert all data to string type
        df_transposed = df_transposed.astype(str)
        logging.info("Data type conversion to string completed")

        # Select only relevant columns
        columns_to_select = ['GPCRdb(A)', ref, 'Seq consensus', 'Seq consensus value']
        df_transposed = df_transposed[columns_to_select]

        # Retrieve and process JSON file for reference protein
        json_file = retrieve_json_file(ref)
        if json_file:
            df2 = json_to_pandas2(json_file)
            logging.info("Reference JSON file processed into DataFrame")

            # Select relevant columns from reference residues table
            df_ref_resides_table = df2[['GPCRdb(A)', f"{ref}_AA-seq"]]
            df_ref_resides_table.to_excel(f"{ref}_GPCRdb_residues_table_full_seq.xlsx")
            # Filter out invalid GPCRdb positions
            df_ref_resides_table = df_ref_resides_table[df_ref_resides_table['GPCRdb(A)'] != 'NonexNone']

            # Merge consensus DataFrame with reference residues table
            df_ref_consensus_with_residues_numbering = df_transposed.merge(df_ref_resides_table, on='GPCRdb(A)', how='left').dropna()
            logging.info("Merged consensus DataFrame with reference residues numbering")

            return df_ref_resides_table, df_ref_consensus_with_residues_numbering
        else:
            logging.warning("Failed to retrieve JSON file for reference protein: %s", ref)
            return pd.DataFrame(), pd.DataFrame()

    except Exception as e:
        logging.error("Error processing alignment: %s", e)
        return pd.DataFrame(), pd.DataFrame()



    
def process_alignment_custom(
    df_with_consensus: pd.DataFrame, ref: str, ref_residues_table: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Processes the consensus DataFrame by transposing it, ensuring consensus rows exist,
    setting headers, and merging with reference residues.

    Parameters
    ----------
    df_with_consensus : pd.DataFrame
        Input DataFrame containing sequence information (with or without consensus rows).
    ref : str
        The reference protein name used to retrieve Excel file and align residue numbering.
    ref_residues_table : str
        Path to Excel file containing reference residue numbering table.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        - df_ref_resides_table : DataFrame with reference residues and their GPCRdb numbering.
        - df_ref_consensus_with_residues_numbering : Merged DataFrame with consensus and reference residue numbering.
    """
    logging.info("Starting processing of alignment with reference: %s", ref)

    try:
        # ------------------------------------------------------------------
        # Ensure consensus rows exist
        # ------------------------------------------------------------------
        if not any(df_with_consensus.iloc[:, 0].str.contains("Seq consensus", case=False, na=False)):
            logging.info("Consensus rows not found – calculating consensus")

            df_with_consensus = add_consensus(df_with_consensus)
        else:
            logging.info("Consensus rows already present – skipping calculation")

        # ------------------------------------------------------------------
        # Transpose DataFrame
        # ------------------------------------------------------------------
        df_transposed = df_with_consensus.T
        logging.debug("DataFrame transposed successfully")

        # First row as header
        df_transposed.columns = df_transposed.iloc[0]
        df_transposed = df_transposed.drop(df_transposed.index[0])
        logging.debug("Header set and first row dropped")

        # Rename index and reset
        df_transposed = df_transposed.rename_axis('GPCRdb(A)').reset_index()

        # Convert everything to string
        df_transposed = df_transposed.astype(str)

        # Select relevant columns
        columns_to_select = ['GPCRdb(A)', ref, 'Seq consensus', 'Seq consensus value']
        df_transposed = df_transposed[columns_to_select]

        # ------------------------------------------------------------------
        # Load reference residues table
        # ------------------------------------------------------------------
        residues_table = pd.read_excel(ref_residues_table, header=0, dtype=str)
        if residues_table is not None:
            df2 = residues_table.iloc[:, 1:]
            logging.info("Reference residue table Excel file processed into DataFrame")

            # Select relevant columns
            df_ref_resides_table = df2[['GPCRdb(A)', f"{ref}_AA-seq"]]

            # Filter invalid GPCRdb positions
            df_ref_resides_table = df_ref_resides_table[
                df_ref_resides_table['GPCRdb(A)'] != 'NonexNone'
            ]

            # Merge consensus with reference residues
            df_ref_consensus_with_residues_numbering = (
                df_transposed.merge(df_ref_resides_table, on='GPCRdb(A)', how='left')
                .dropna()
            )

            logging.info("Merged consensus DataFrame with reference residues numbering")
            return df_ref_resides_table, df_ref_consensus_with_residues_numbering

        else:
            logging.warning("Failed to load reference residue table for: %s", ref)
            return pd.DataFrame(), pd.DataFrame()

    except Exception as e:
        logging.error("Error processing alignment: %s", e)
        return pd.DataFrame(), pd.DataFrame()
    

def GPRCs_Variance(
    ref: str, 
    df_ref_consensus_with_residues_numbering: pd.DataFrame, 
    target: str, 
    df_target_consensus_with_residues_numbering: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Processes and merges consensus DataFrames for reference and target GPCRs, 
    adjusting column names and combining consensus values and sequences.

    Parameters
    ----------
    ref : str
        Name of the reference protein.
    df_ref_consensus_with_residues_numbering : pd.DataFrame
        DataFrame containing consensus data for the reference protein.
    target : str
        Name of the target protein.
    df_target_consensus_with_residues_numbering : pd.DataFrame
        DataFrame containing consensus data for the target protein.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        - df_ref_consensus_with_residues_numbering : Processed DataFrame with renamed columns for reference protein.
        - df_target_consensus_with_residues_numbering : Processed DataFrame with renamed columns for target protein.
        - df_ref_target_consensus_with_residues_numbering : Merged DataFrame containing consensus data for both proteins.
    """
    try:
        logging.info("Renaming columns in reference and target DataFrames for merging.")

        # Rename columns in the reference DataFrame
        ref_rename_map = {
            "GPCRdb(A)" + ref: "GPCRdb(A)",
            "Seq consensus value": "Seq consensus value " + ref,
            "Seq consensus": "Seq consensus " + ref
        }
        df_ref_consensus_with_residues_numbering.rename(columns=ref_rename_map, inplace=True)

        # Rename columns in the target DataFrame
        target_rename_map = {
            "GPCRdb(A)" + target: "GPCRdb(A)",
            "Seq consensus value": "Seq consensus value " + target,
            "Seq consensus": "Seq consensus " + target
        }
        df_target_consensus_with_residues_numbering.rename(columns=target_rename_map, inplace=True)

        logging.info("Merging reference and target consensus DataFrames on 'GPCRdb(A)' column.")
        
        # Merge reference and target DataFrames on "GPCRdb(A)"
        df_ref_target_consensus_with_residues_numbering = pd.merge(
            df_ref_consensus_with_residues_numbering,
            df_target_consensus_with_residues_numbering,
            on="GPCRdb(A)"
        )

        logging.info("Merge successful. Combined DataFrame created.")
        
        return (
            df_ref_consensus_with_residues_numbering,
            df_target_consensus_with_residues_numbering,
            df_ref_target_consensus_with_residues_numbering
        )

    except Exception as e:
        logging.error("Error occurred while processing GPCRs variance: %s", e)
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()


def get_AF(protein_name: str) -> str:
    """
    Fetches the AlphaFold prediction PDB URL for a given protein using its UniProt accession number from GPCRdb
    and saves the PDB file using the protein name as the filename.

    Parameters
    ----------
    protein_name : str
        The name of the protein (e.g., 'adrb2_human').

    Returns
    -------
    str
        The local path where the PDB file is saved, or None if an error occurs.
    """
    try:
        # Step 1: Fetch UniProt accession number from GPCRdb
        protein_url = f"https://gpcrdb.org/services/protein/{protein_name}/"
        with urlopen(protein_url) as response:
            protein_data = json.loads(response.read().decode('utf-8'))
        accession = protein_data.get('accession')

        if not accession:
            logging.warning("Accession not found for protein %s", protein_name)
            return None
        logging.info("UniProt accession %s retrieved for protein %s", accession, protein_name)

        # Step 2: Fetch AlphaFold prediction data using the UniProt accession
        AF_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
        with urlopen(AF_url) as response:
            content = response.read().decode('utf-8')

        try:
            AF_data = json.loads(content)
        except json.JSONDecodeError as e:
            logging.error("Failed to parse AlphaFold JSON: %s", e)
            return None

        # Step 3: Extract PDB URL from AlphaFold data
        if AF_data and isinstance(AF_data, list) and len(AF_data) > 0:
            pdb_url = AF_data[0].get("pdbUrl")
            if pdb_url:
                # Step 4: Download and save PDB file
                file_name = f"{protein_name}.pdb"
                file_path = os.path.join(os.getcwd(), file_name)
                
                if os.path.exists(file_path):
                    os.remove(file_path)
                    logging.info("Existing PDB file %s deleted.", file_name)

                urlretrieve(pdb_url, file_path)
                logging.info("PDB file saved as %s", file_path)
                return file_path
            else:
                logging.warning("No PDB URL found for AlphaFold data for protein %s", protein_name)
                return None
        else:
            logging.warning("No AlphaFold data found for protein %s", protein_name)
            return None

    except Exception as e:
        logging.error("An error occurred while fetching AlphaFold data for %s: %s", protein_name, e)
        return None


def coexistence(rec: str, residue: str, cutoff: float) -> list:
    """
    Identifies residues within a specified cutoff distance of a given residue in a loaded structure.

    Parameters
    ----------
    rec : str
        The filename or path of the receptor structure to load.
    residue : int
        The residue number around which to search for nearby residues.
    cutoff : float
        The distance cutoff (in angstroms) within which residues are considered "nearby."

    Returns
    -------
    list
        A list of nearby residues in the format 'residue name and number', e.g., ['ASP45', 'GLY46'].
    """
    try:
        #logging.info("Starting coexistence function for receptor: %s, residue: %s, cutoff: %.2f", rec, residue, cutoff)
        
        # Initialize list to store selected residues
        stored.selected_residues = []
        
        # Load the receptor structure
        cmd.load(rec, "rec")
        logging.debug("Receptor %s loaded into PyMOL", rec)
        
        # Select residues within the specified cutoff distance
        selection_name = 'residues_within_cutoff'
        cmd.select(selection_name, f"br. all within {cutoff} of resi {residue}")
        logging.debug("Residues selected within cutoff distance of %.2f from residue %d", cutoff, residue)
        
        # Iterate over selected residues and store them
        cmd.iterate(f"{selection_name} and name CA", "stored.selected_residues.append(f'{resn}{resi}')")
        
        # Reinitialize to clean up the session
        cmd.reinitialize()
        #logging.info("Session reinitialized, selected residues stored")

        return stored.selected_residues

    except Exception as e:
        logging.error("Error in coexistence function: %s", e)
        return []


def fetch_and_combine(protein_list: List[str], segments: List[str]) -> pd.DataFrame:
    """
    Fetches GPCRdb positions for specified transmembrane (TM) and loop segments across a list of proteins,
    aligns the data, and combines results with consensus sequence information.

    Parameters
    ----------
    protein_list : List[str]
        A list of protein names to retrieve data for.
    segments : List[str]
        A list of segments (e.g., ['TM1', 'TM2', 'ECL1', 'ICL2']) to fetch positions for.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with consensus sequence information, or None if no data is available.
    """
    try:
        logging.info("Separating segments into TM and loop categories.")
        TM_segments = [seg for seg in segments if 'TM' in seg]
        loop_segments = [seg for seg in segments if 'ECL' in seg or 'ICL' in seg]

        # Fetch positions for TM segments
        df_tm = None
        if TM_segments:
            logging.info("Fetching positions for TM segments: %s", TM_segments)
            protein_list_tm, positions_tm = fetch_species_positions(protein_list, TM_segments)
            df_tm = fetch_gpcrdb_alignment(protein_list_tm, positions_tm)
            logging.debug("TM segments alignment data fetched successfully")

        # Fetch positions for loop segments
        df_loop = None
        if loop_segments:
            logging.info("Fetching positions for loop segments: %s", loop_segments)
            protein_list_loop, positions_loop = fetch_species_positions(protein_list, loop_segments)
            df_loop = fetch_gpcrdb_alignment(protein_list_loop, positions_loop)
            logging.debug("Loop segments alignment data fetched successfully")

        # Combine TM and loop DataFrames if both are present
        if df_tm is not None and df_loop is not None:
            df_combined = pd.merge(df_tm, df_loop, on="GPCRdb(A)")
            logging.info("TM and loop DataFrames merged successfully")
        elif df_tm is not None:
            df_combined = df_tm
            logging.info("Only TM segments data is available and used")
        elif df_loop is not None:
            df_combined = df_loop
            logging.info("Only loop segments data is available and used")
        else:
            logging.warning("No data available for the specified segments")
            return None

        # Add consensus sequence information if combined DataFrame exists
        if df_combined is not None:
            logging.info("Adding consensus sequence information to the combined DataFrame")
            df_with_consensus = add_consensus(df_combined)
            return df_with_consensus
        else:
            logging.warning("Combined DataFrame is empty, no consensus added")
            return None

    except Exception as e:
        logging.error("An error occurred during fetch and combine: %s", e)
        return None

def PDB_read_and_edit_noconsensus_coexistence2(
    ref: str, target: str, df_ref_target_consensus_with_residues_numbering: pd.DataFrame,
    conservation_cutoff: float, substitution_matrix_method: str
) -> None:
    """
    Analyzes and modifies PDB structures based on consensus information and substitution matrix.

    Parameters
    ----------
    ref : str
        Reference protein name.
    target : str
        Target protein name.
    df_ref_target_consensus_with_residues_numbering : pd.DataFrame
        DataFrame containing consensus sequence information for reference and target proteins.
    conservation_cutoff : float
        Conservation cutoff threshold to identify significant amino acids.
    substitution_matrix_method : str
        Name of the substitution matrix to use for calculating residue similarity.
    """
    logging.info("Starting PDB analysis with conservation cutoff: %s and substitution matrix: %s", conservation_cutoff, substitution_matrix_method)

    try:
        # Load substitution matrix
        substitution_matrix = substitution_matrices.load(substitution_matrix_method)
        dist_cutoff = 4.0
        mutations = []
        mutations_ref_to_target = []
        mutations_target_to_ref = []
        important_AA_index = []
        AA_ref_coexistence_list_indexes = []
        AA_target_coexistence_list_indexes = []

        # Load structures
        structure_r=get_AF(ref)
        structure_t=get_AF(target)
        structure_ref = PDBParser(QUIET=True).get_structure(ref, get_AF(ref))
        structure_target = PDBParser(QUIET=True).get_structure(target, get_AF(target))
        logging.info("Structures for %s and %s loaded successfully", ref, target)

        # Set up columns for mutation annotations
        df_ref_target_consensus_with_residues_numbering[f"{ref}_to_{target}"] = 'nan'
        df_ref_target_consensus_with_residues_numbering[f"{target}_to_{ref}"] = 'nan'
        df_ref_target_consensus_with_residues_numbering[f"{ref}_to_color_{target}"] = 'white'
        df_ref_target_consensus_with_residues_numbering[f"{target}_to_color_{ref}"] = 'white'
        df_ref_target_consensus_with_residues_numbering["Mutation"] = 'nan'
        df_ref_target_consensus_with_residues_numbering[f"{substitution_matrix_method}_score"] = 'nan'
        

        for AA_no, row in df_ref_target_consensus_with_residues_numbering.iterrows():
            AA_ref_consensus = row[f"Seq consensus {ref}"]
            AA_target_consensus = row[f"Seq consensus {target}"]
            AA_ref = row.iloc[4]
            AA_target = row.iloc[8]
            AA_ref_consensus_value = float(row[f'Seq consensus value {ref}'])
            AA_target_consensus_value = float(row[f'Seq consensus value {target}'])
            AA_GPRCdb_number = float(row['GPCRdb(A)'].replace('x', '.'))
            
            #rest bfactors
            structure_ref[0]['A'][int(AA_ref[1:])]["CA"].bfactor = 0
            structure_target[0]['A'][int(AA_target[1:])]["CA"].bfactor = 0
            # Assign b-factors to CA atoms for visualization
            structure_ref[0]['A'][int(AA_ref[1:])]["C"].bfactor = AA_GPRCdb_number
            structure_target[0]['A'][int(AA_target[1:])]["C"].bfactor = AA_GPRCdb_number
            structure_ref[0]['A'][int(AA_ref[1:])]["N"].bfactor=AA_ref_consensus_value
            structure_target[0]['A'][int(AA_target[1:])]["N"].bfactor=AA_target_consensus_value
            df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"Mutation"] = (f"{AA_ref[0]}↔{AA_target[0]}")
            df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{substitution_matrix_method}_score"]= substitution_matrix[AA_ref[0], AA_target[0]]
            if substitution_matrix[AA_ref[0], AA_target[0]] < 0:
                if (AA_ref_consensus_value >= conservation_cutoff and 
                    AA_target_consensus_value >= conservation_cutoff and 
                    AA_ref[0] == AA_ref_consensus[0] and 
                    AA_target[0] == AA_target_consensus[0]):

                    structure_ref[0]['A'][int(AA_ref[1:])]["CA"].bfactor = 100
                    structure_target[0]['A'][int(AA_target[1:])]["CA"].bfactor = 100
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_{target}"] = AA_ref + AA_target[0]
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_{ref}"] = AA_target + AA_ref[0]
                    # row[f"{ref}_to_{target}"] = AA_ref + AA_target[0]
                    # row[f"{target}_to_{ref}"] = AA_target + AA_ref[0]
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_color_{target}"] = 'red'
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_color_{ref}"] = 'red'
                    
                    mutations.append(AA_ref)
                    important_AA_index.append(AA_no)
                    mutations_ref_to_target.append(AA_ref + 'to' + AA_target)
                    mutations_target_to_ref.append(AA_target + 'to' + AA_ref)
                    AA_ref_coexistence_list=coexistence(structure_r,AA_ref[1:],dist_cutoff)
                    for AA in AA_ref_coexistence_list:
                        AA_index_ref=df_ref_target_consensus_with_residues_numbering.index[df_ref_target_consensus_with_residues_numbering[df_ref_target_consensus_with_residues_numbering.columns[4]]== rename_AA(AA[:3])+AA[3:]]                
                        if len(AA_index_ref)!=0:
                            AA_ref_coexistence_list_indexes.append(AA_index_ref[0])       
                    AA_target_coexistence_list=coexistence(structure_t,AA_target[1:],dist_cutoff)
                    for AA in AA_target_coexistence_list:
                        AA_index_target=df_ref_target_consensus_with_residues_numbering.index[df_ref_target_consensus_with_residues_numbering[df_ref_target_consensus_with_residues_numbering.columns[8]]== rename_AA(AA[:3])+AA[3:]]
                        if len(AA_index_target)!=0:
                            AA_target_coexistence_list_indexes.append(AA_index_target[0])  

                    
                    #logging.info("|11111| %s (%s) -> %s (%s)", AA_ref, AA_ref_consensus_value, AA_target, AA_target_consensus_value)

                elif AA_ref_consensus_value >= conservation_cutoff or AA_target_consensus_value >= conservation_cutoff:
                    if AA_ref_consensus_value >= conservation_cutoff and AA_ref[0] == AA_ref_consensus[0]:
                        structure_ref[0]['A'][int(AA_ref[1:])]["CA"].bfactor = 50
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_{target}"] = AA_ref + AA_target[0]
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_color_{target}"]= 'orange'
                        AA_ref_coexistence_list=coexistence(structure_r,AA_ref[1:],dist_cutoff)
                        for AA in AA_ref_coexistence_list:
                            AA_index_ref=df_ref_target_consensus_with_residues_numbering.index[df_ref_target_consensus_with_residues_numbering[df_ref_target_consensus_with_residues_numbering.columns[4]]== rename_AA(AA[:3])+AA[3:]]                
                            if len(AA_index_ref)!=0:
                                AA_ref_coexistence_list_indexes.append(AA_index_ref[0])  
                        #logging.info("|22222| %s (%s)", AA_ref, AA_ref_consensus_value)

                    if AA_target_consensus_value >= conservation_cutoff and AA_target[0] == AA_target_consensus[0]:
                        structure_target[0]['A'][int(AA_target[1:])]["CA"].bfactor = 50
                        row[f"{ref}_to_{target}"] = AA_ref + AA_target[0]
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_{ref}"] = AA_target + AA_ref[0]
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_color_{ref}"]= 'orange'
                        AA_target_coexistence_list=coexistence(structure_t,AA_target[1:],dist_cutoff)
                        for AA in AA_target_coexistence_list:
                            AA_index_target=df_ref_target_consensus_with_residues_numbering.index[df_ref_target_consensus_with_residues_numbering[df_ref_target_consensus_with_residues_numbering.columns[8]]== rename_AA(AA[:3])+AA[3:]]
                            if len(AA_index_target)!=0:
                                AA_target_coexistence_list_indexes.append(AA_index_target[0])  
                            
                    mutations.append(AA_ref)
                    important_AA_index.append(AA_no)
            else:
                structure_ref[0]['A'][int(AA_ref[1:])]["CA"].bfactor = 0
                structure_target[0]['A'][int(AA_target[1:])]["CA"].bfactor = 0
                row[f"{ref}_to_{target}"] = "Inconclusive"
                row[f"{target}_to_{ref}"] = "Inconclusive"
                logging.debug("Inconclusive AA conservation at %s -> %s", AA_ref, AA_target)

        logging.info("%s to %s: conservation_cutoff %s %s - No. of mutations: %d", ref, target, conservation_cutoff, substitution_matrix_method, len(mutations_ref_to_target))
        logging.info("%s to %s: conservation_cutoff %s %s - No. of mutations: %d", target, ref, conservation_cutoff, substitution_matrix_method, len(mutations_target_to_ref))
        output1 = (
        f"{ref} to {target}: conservation_cutoff {conservation_cutoff} {substitution_matrix_method} - No. of mutations: {len(mutations_ref_to_target)}\n"
    )
        output2 = (
        f"{target} to {ref}: conservation_cutoff {conservation_cutoff} {substitution_matrix_method} - No. of mutations: {len(mutations_target_to_ref)}\n"
    )
# Write the results to a file            
        with open(f"{ref}_to_{target}_mutations_number.txt" , "w") as file:
            file.write(output1)
        with open(f"{target}_to_{ref}_mutations_number.txt" , "w") as file:
            file.write(output2)
                        #### Make list of the coexistences AA indices 
        AA_ref_target_coexistence_list_indexs=(list(set(AA_ref_coexistence_list_indexes + AA_target_coexistence_list_indexes))) 
        AA_ref_target_coexistence_list_indexs_not_in_important=[]
        for i in AA_ref_target_coexistence_list_indexs:
            if i not in important_AA_index:
                AA_ref_target_coexistence_list_indexs_not_in_important.append(i)
                
                        ###Coexistance
        for AA_no in AA_ref_target_coexistence_list_indexs_not_in_important: 
            AA_ref=df_ref_target_consensus_with_residues_numbering.iloc[AA_no,4]
            AA_target=df_ref_target_consensus_with_residues_numbering.iloc[AA_no,8]
            AArefconsensus_value=float(df_ref_target_consensus_with_residues_numbering.iloc[AA_no,df_ref_target_consensus_with_residues_numbering.columns.get_loc('Seq consensus value '+ ref)])
            AAtargetconsensus_value=float(df_ref_target_consensus_with_residues_numbering.iloc[AA_no,df_ref_target_consensus_with_residues_numbering.columns.get_loc('Seq consensus value '+ target)])
            
            if AArefconsensus_value >=conservation_cutoff  and AAtargetconsensus_value >=conservation_cutoff:
                if  AA_ref[0] !=AA_target[0]:
                    structure_ref[0]['A'][int(AA_ref[1:])]["CA"].bfactor=25
                    structure_target[0]['A'][int(AA_target[1:])]["CA"].bfactor=25
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_{target}"] = AA_ref + AA_target[0]
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_{ref}"] = AA_target + AA_ref[0]
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_color_{target}"]='blue'
                    df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_color_{ref}"]='blue'
                    
            elif AArefconsensus_value >=conservation_cutoff  or AAtargetconsensus_value >=conservation_cutoff:
                if  AA_ref[0] !=AA_target[0]:
                    if AArefconsensus_value >=conservation_cutoff:
                        structure_ref[0]['A'][int(AA_ref[1:])]["CA"].bfactor=12.5
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_{ref}"] = AA_target + AA_ref[0]
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_color_{target}"]='cyan'
                    if AAtargetconsensus_value >=conservation_cutoff:
                        structure_target[0]['A'][int(AA_target[1:])]["CA"].bfactor=12.5
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{ref}_to_{target}"] = AA_ref + AA_target[0]
                        df_ref_target_consensus_with_residues_numbering.loc[AA_no, f"{target}_to_color_{ref}"]='cyan'
                        
        io = PDBIO()
        io.set_structure(structure_ref)
        io.save(f"{ref}_{conservation_cutoff}_{substitution_matrix_method}.pdb")
        io.set_structure(structure_target)
        io.save(f"{target}_{conservation_cutoff}_{substitution_matrix_method}.pdb")
        logging.info("Structures saved with conservation cutoff %s and substitution matrix %s", conservation_cutoff, substitution_matrix_method)

        filename_mutagenesis_all_positions = f"{ref}_{target}_{substitution_matrix_method}_{conservation_cutoff}_mutagenesis_all_positions.xlsx"
        df_ref_target_consensus_with_residues_numbering.to_excel(filename_mutagenesis_all_positions)
        logging.info("Analysis stored in: %s", filename_mutagenesis_all_positions)

    except Exception as e:
        logging.error("Error occurred in PDB_read_and_edit_noconsensus_coexistence2: %s", e)

def make_pymol_session(pdb1: str, pdb2: str, output_session: str) -> None:
    """
    Creates a PyMOL session by loading, aligning, and color-coding two PDB structures based on B-factors.
    
    Parameters
    ----------
    pdb1 : str
        The filename of the first PDB structure.
    pdb2 : str
        The filename of the second PDB structure.
    output_session : str
        The filename to save the PyMOL session.
    """
    logging.info("Starting PyMOL session creation with %s and %s", pdb1, pdb2)

    try:
        # Load the two PDB files from AlphaFold
        prot1 = pdb1
        prot2 = pdb2
        
        cmd.load(pdb1, prot1)
        cmd.load(pdb2, prot2)
        logging.info("PDB files %s and %s loaded as %s and %s", pdb1, pdb2, prot1, prot2)

        # Align the two structures
        cmd.align(prot1, prot2)
        logging.debug("Structures %s and %s aligned", prot1, prot2)

        # Color both structures white by element
        cmd.color("white", f"{prot1} and elem c")
        cmd.color("white", f"{prot2} and elem c")
        cmd.show("lines", "sidechain or name CA")
        logging.debug("Both structures colored white and displayed as lines")

        # Apply color coding based on B-factors
        # Select CA atoms with specific B-factors, show as sticks, and apply colors
        selections = [
            ("conserved_in_one", "br. name CA and b=12.5", "cyan"),
            ("conserved_in_both", "br. name CA and b=25", "blue"),
            ("conserved_and_different_in_one", "br. name CA and b=50", "orange"),
            ("conserved_and_different_in_both", "br. name CA and b=100", "red")
        ]

        for name, selection, color in selections:
            cmd.select(name, selection)
            cmd.show("sticks", name)
            cmd.color(color, f"{name} and elem c")
            logging.debug("Selection %s with B-factor criteria '%s' colored %s", name, selection, color)

        # Save the session
                
        # Export aligned PDB files
        aligned_pdb1 = "prot1_aligned.pdb"
        aligned_pdb2 = "prot2_aligned.pdb"
        logging.info("Aligned PDB files saved as %s and %s", aligned_pdb1, aligned_pdb2)
        
        cmd.save(aligned_pdb1, prot1)
        cmd.save(aligned_pdb2, prot2)
        cmd.save(output_session)
        logging.info("PyMOL session saved as %s", output_session)

    except Exception as e:
        logging.error("Error occurred in make_pymol_session: %s", e)

def main():
    parser = argparse.ArgumentParser(description=" Example to run: python GPCR_Selectivity_Explorer.py  -l1 adrb1_human adrb2_human adrb3_human -r1 adrb2_human -l2 drd1_human drd5_human -r2 drd1_human -sg TM1 TM2 TM3 TM4 TM5 TM6 TM7 ECL2 -c 80 -sb BLOSUM62")
    parser.add_argument("-l1", "--list1", nargs='+', default=None, help="Space-separated list for Reference 1 Receptor segments (e.g., adrb1_human adrb2_human adrb3_human)")
    parser.add_argument("-r1","--reference1",type=str,default = None, help = "Reference 1 Receptor Name in GPCRdb (e.g., adrb2_human)") 
    parser.add_argument("-l2", "--list2", nargs='+', default=None, help="Space-separated list for Reference 2 Receptor segments (e.g., drd1_human drd5_human)")
    parser.add_argument("-r2","--reference2",type=str,default = None, help = "Reference 2 Receptor Name in GPCRdb (e.g., drd1_human)")
    parser.add_argument("-sg", "--segments", nargs='+', default=None, help="Segment information for reference receptor residues from GPCRdb (e.g., TM1 TM2 ECL2 ICL1)")
    parser.add_argument("-c","--conservation_cutoff",type=int,default = 50,required=False,help = "Conservation Cutoff ")
    parser.add_argument("-sb", "--substitution_matrix_method", type=str, default="BLOSUM62", required=False, 
                    help=("Substitution matrix method. Options include: "
                    "BENNER22, BENNER6, BENNER74, BLOSUM45, BLOSUM50, BLOSUM62, "
                    "BLOSUM80, BLOSUM90, DAYHOFF, FENG, GENETIC, GONNET1992, "
                    "HOXD70, JOHNSON, JONES, LEVIN, MCLACHLAN, MDM78, "
                    "NUC.4.4, PAM250, PAM30, PAM70, RAO, RISLER, "
                    "SCHNEIDER, STR, TRANS"))

    parser.add_argument("-o","--outputfolder", type=str, default = 'outputs', required=False, help = "Output folder")
    
    parser.add_argument("--custom", action="store_true", required=False, help = "use custom alginments and residue tables")
    parser.add_argument("-r1f","--reference1_residues_table",type=str,default = None, help = "Reference 1 modifed residues table excel file from GPCRdb ")
    parser.add_argument("-r1A","--reference1_alignment",type=str,default = None, help = "Reference 1 modifed alignment excel file from GPCRdb")
    parser.add_argument("-r2f","--reference2_residues_table",type=str,default = None, help = "Reference 2 modifed residues table excel file from  GPCRdb ")
    parser.add_argument("-r2A","--reference2_alignment",type=str,default = None, help = "Reference 2 modifed alignment excel file from GPCRdb")


    args = parser.parse_args()
    if __name__ == "__main__":
        
        
        # Convert input files and output folder to absolute paths
        args.reference1_alignment = os.path.abspath(args.reference1_alignment) if args.reference1_alignment else None
        args.reference2_alignment = os.path.abspath(args.reference2_alignment) if args.reference2_alignment else None
        args.reference1_residues_table = os.path.abspath(args.reference1_residues_table) if args.reference1_residues_table else None
        args.reference2_residues_table = os.path.abspath(args.reference2_residues_table) if args.reference2_residues_table else None
        args.outputfolder = os.path.abspath(args.outputfolder)


        with change_dir(args.outputfolder):
            if not args.custom:
                df_with_consensus_r1 = fetch_and_combine(args.list1,args.segments)
                df_with_consensus_r2 = fetch_and_combine(args.list2,args.segments)
                
                df_with_consensus_r1.to_excel(f"{args.reference1}_GPCRdb_alignment.xlsx")
                df_with_consensus_r2.to_excel(f"{args.reference2}_GPCRdb_alignment.xlsx")
                
                df_resides_tabel_r1, df_consensus_with_residues_numbering_r1= process_alignment(df_with_consensus_r1,args.reference1)
                df_resides_tabel_r1.to_excel(f"{args.reference1}_GPCRdb_residues_table.xlsx")
                df_consensus_with_residues_numbering_r1.to_excel(f"{args.reference1}_consensus_with_residues_numbering.xlsx")

                df_resides_tabel_r2, df_consensus_with_residues_numbering_r2=process_alignment(df_with_consensus_r2,args.reference2)
                df_resides_tabel_r2.to_excel(f"{args.reference2}_GPCRdb_residues_table.xlsx")
                df_consensus_with_residues_numbering_r2.to_excel(f"{args.reference2}_consensus_with_residues_numbering.xlsx")   
                
                df_consensus_with_resdiues_numbering_r_t=GPRCs_Variance(args.reference1,df_consensus_with_residues_numbering_r1,args.reference2,df_consensus_with_residues_numbering_r2)[2]
                PDB_read_and_edit_noconsensus_coexistence2(args.reference1,args.reference2,df_consensus_with_resdiues_numbering_r_t,args.conservation_cutoff,args.substitution_matrix_method)
                make_pymol_session(args.reference1+'_'+str(args.conservation_cutoff)+'_'+args.substitution_matrix_method+'.pdb',args.reference2+'_'+str(args.conservation_cutoff)+'_'+args.substitution_matrix_method+'.pdb','pymol_'+f"{args.reference1}_{args.reference2}.pse")

            if args.custom:
                df_with_consensus_r1=pd.read_excel(args.reference1_alignment,index_col=0, header=0,dtype=str)
                df_with_consensus_r2=pd.read_excel(args.reference2_alignment,index_col=0, header=0,dtype=str)

                df_resides_tabel_r1, df_consensus_with_residues_numbering_r1= process_alignment_custom(df_with_consensus_r1,args.reference1,args.reference1_residues_table)
                df_resides_tabel_r1.to_excel(f"{args.reference1}_GPCRdb_residues_table2.xlsx")
                df_consensus_with_residues_numbering_r1.to_excel(f"{args.reference1}_consensus_with_residues_numbering2.xlsx")
                
                df_resides_tabel_r2, df_consensus_with_residues_numbering_r2=process_alignment_custom(df_with_consensus_r2,args.reference2,args.reference2_residues_table)
                df_resides_tabel_r2.to_excel(f"{args.reference2}_GPCRdb_residues_table2.xlsx")
                df_consensus_with_residues_numbering_r2.to_excel(f"{args.reference2}_consensus_with_residues_numbering2.xlsx")  
                
                df_consensus_with_resdiues_numbering_r_t=GPRCs_Variance(args.reference1,df_consensus_with_residues_numbering_r1,args.reference2,df_consensus_with_residues_numbering_r2)[2]
                PDB_read_and_edit_noconsensus_coexistence2(args.reference1,args.reference2,df_consensus_with_resdiues_numbering_r_t,args.conservation_cutoff,args.substitution_matrix_method)
                make_pymol_session(args.reference1+'_'+str(args.conservation_cutoff)+'_'+args.substitution_matrix_method+'.pdb',args.reference2+'_'+str(args.conservation_cutoff)+'_'+args.substitution_matrix_method+'.pdb','pymol_'+f"{args.reference1}_{args.reference2}.pse")

                
                
if __name__ == "__main__":
    main()

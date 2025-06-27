import os
import re
import pandas as pd

########################### Extract Specie and bait name ######################
def extract_species_bait(filepath):
    """
    Extracts the species and bait names from a filename formatted as 'species_bait_data.csv'.
    Example: 'human_src_data.csv' -> ('Human', 'SRC')
    Parameters:
        filename (str): The name or path of the data file
    Returns:
        tuple: (species_name (str), bait_name (str))
    Raises:
        ValueError: If the filename does not match the expected pattern
    """
    name = os.path.basename(filepath)
    name = os.path.splitext(name)[0]

    # Match expected pattern: species_bait_data or species_bait_table, etc.
    match = re.match(r'^([a-zA-Z]+)_([a-zA-Z0-9]+)_(data|table|)?$', name, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"❌ Filename '{filepath}' doesn't match expected pattern 'species_bait_data.csv'")

    species, bait = match.group(1), match.group(2)
    return species.capitalize(), bait.upper()

################## File Operations: Read, Merge, Save DataFrames ####################

def read_csv_file(file_path, header='infer'):
    """
    Read a CSV file and return a Pandas DataFrame.
    Parameters:
        file_path (str): Path to the CSV file to read
        header (int, optional): Row number to use as column headers
    Returns:
        DataFrame: Pandas DataFrame containing the read data
    """

    try:
        df = pd.read_csv(file_path, sep=',', header=header)
        print(f"File {file_path} loaded successfully.")
        return df
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None 


def save_csv_file(df, output_path):
    """
    Save a DataFrame to a CSV file.
    Parameters:
        df (DataFrame): DataFrame to save
        output_path (str): Path where the CSV file will be saved
    Returns:
        bool: True if successful, False otherwise
    """

    try:
        df.to_csv(output_path, index=False)
        print(f"✅File saved to {output_path}.")

    except Exception as e:
        print(f"❌Error saving file to {output_path}: {e}")
    
############################  Data Merge #######################################

def merge_files(file_list, output_path=None):
    """
    Merge multiple data files into one multi-indexed DataFrame.
    Each column gets a MultiIndex: (Species, Bait, OriginalColumn)

    Parameters:
        file_list (list): List of CSV file paths
        output_path (str): Output path to save the merged file

    Returns:
        DataFrame: Merged DataFrame with MultiIndex columns
    """
    all_dfs = []

    for path in file_list:
        species, bait = extract_species_bait(path)
        df = read_csv_file(path)
        if df is None:
            continue
        df.columns = pd.MultiIndex.from_product([[species], [bait], df.columns])
        all_dfs.append(df)

    # Concaténer horizontalement
    df_merged = pd.concat(all_dfs, axis=1)

    # Sauvegarde si souhaitée
    if output_path:
        df_merged.to_csv(output_path, index=False)
        print(f"✅ Fichier combiné sauvegardé : {output_path}")

    return df_merged
    

######################### reorganize dataframe columns ##########################

def group_columns_by_species(df):
    """
    Automatically reorder MultiIndex columns in a DataFrame by grouping all columns
    belonging to the same species together.

    Parameters:
        df (DataFrame): A pandas DataFrame with MultiIndex columns (species, bait, metric)

    Returns:
        DataFrame: Reordered DataFrame with species columns grouped
    """
    if not isinstance(df.columns, pd.MultiIndex):
        raise ValueError("❌ The DataFrame must have MultiIndex columns")

    # Get unique species in order of appearance
    species_list = []
    for col in df.columns:
        species = col[0]
        if species not in species_list:
            species_list.append(species)

    # Reorder columns by species
    ordered_cols = []
    for species in species_list:
        species_cols = [col for col in df.columns if col[0] == species]
        ordered_cols.extend(species_cols)

    return df[ordered_cols]

import pandas as pd




################## File Operations: Read, Merge, Save DataFrames ####################

def read_csv_file(file_path,header=None):

    """
    Read a CSV file and return a Pandas DataFrame.
    
    Parameters:
        file_path (str): Path to the CSV file to read
        header (int, optional): Row number to use as column headers
        
    Returns:
        DataFrame: Pandas DataFrame containing the read data
    """

    try:
        df = pd.read_csv(file_path, sep=',')
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
        print(f"File saved to {output_path}.")
        return True
    except Exception as e:
        print(f"Error saving file to {output_path}: {e}")
        return False
    
############################ Species Data Matching #######################################

def create_merged_species_dataset(human_path, mouse_path, output_path):
    """
    Merge human and mouse datasets into a single file with multi-level column headers.
    
    Creates a MultiIndex on columns:
    - Level 0: 'Human' or 'Mouse'
    - Level 1: Original column names
    
    Parameters:
        human_path (str): Path to the human CSV file
        mouse_path (str): Path to the mouse CSV file
        output_path (str): Path to save the combined file
    
    Returns:
        str: Path to the saved merged file
    """


    # Read the files
    df_Human = pd.read_csv(human_path)
    df_Mouse = pd.read_csv(mouse_path)

    # Check: same number of rows
    if df_Human.shape[0] != df_Mouse.shape[0]:
        raise ValueError("Files don't have the same number of rows.")

    # Apply MultiIndex to columns: (species, column_name)
    df_Human.columns = pd.MultiIndex.from_product([['Human'], df_Human.columns])
    df_Mouse.columns = pd.MultiIndex.from_product([['Mouse'], df_Mouse.columns])

    # Merge horizontally    
    df_combined = pd.concat([df_Human, df_Mouse], axis=1)

    # Export CSV with two header rows
    df_combined.columns = pd.MultiIndex.from_tuples(df_combined.columns)
    df_combined.to_csv(output_path, index=False)

    print(f"Combined file successfully saved: {output_path}")
    return output_path

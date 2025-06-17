import pandas as pd
import numpy as np
from file_op import read_csv_file, extract_species_and_bait


####################### Calculate Kapp #####################

def calculate_Kapp(df):
    """
    Calculate apparent binding constants (Kapp) for bait from pKapp values
    Parameters:
        df (DataFrame): DataFrame containing pKapp values for the bait
    Returns:
        DataFrame: Updated DataFrame with additional columns
    """

    df = df.copy()  # Create a copy to avoid modifying the original
    
    # Calculate Kapp
    if 'pKapp' not in df.columns:
        raise ValueError("Missing 'pKapp' column in dataframe")
    df['Kapp'] = 10 ** (-df['pKapp']) * 1e6
    return df
    

####################### Get Bait Concentration ####################
def get_bait_concentration(df, uniprot_id):
    """
    Retrieve the concentration (μM) of a protein based on its Uniprot ID.
    Parameters:
        df (DataFrame): DataFrame containing protein concentration data
        uniprot_id (str): Uniprot ID of the protein
    Returns:
        float: Concentration in μM, or None if not found
    """
    result = df[df['uniprot_id'] == uniprot_id]['[µM]'].values
    return result[0] if len(result) > 0 else None


#################### Complexome Calculations ##########################################
def calculate_complex_fraction(total_A, total_B, dissociation_constant):
    """
    Calculate the fraction of complex formed between two proteins A and B.
    Parameters:
        total_A (float): Total concentration of protein A (μM)
        total_B (float): Total concentration of protein B (μM)
        dissociation_constant (float): Kd value for the interaction (μM)
    Returns:
        float: Concentration of the formed complex (μM)
    """

    return ((total_A + total_B + dissociation_constant) - 
            np.sqrt((total_A + total_B + dissociation_constant)**2 - 4 * total_A * total_B)) / 2


def add_complexome_column(df, bait_conc):
    """
    Add a column to the DataFrame with complexome calculations for a specific bait.
    Parameters:
        df (DataFrame): DataFrame containing interaction data
        bait_conc (float): Total concentration of the bait protein (μM)
        bait (str): Name of the bait protein (e.g., 'SRC' or 'HCK')
    Returns:
        DataFrame: Updated DataFrame with new complexome column
    """
    df=df.copy()  # Create a copy to avoid modifying the original
    df['complexome'] = calculate_complex_fraction(bait_conc, df['[µM]'], df['Kapp'])
    return df


import pandas as pd
import numpy as np
from lib.file_op import extract_species_bait

####################### Get Bait Concentration ####################

def get_bait_concentration(df, filepath):
    """
    Extract bait name from filepath, find the corresponding Uniprot ID in the DataFrame,
    and retrieve its concentration.
    Searches for the bait in the 'gene_id' column using a case- and space-insensitive match.

    Parameters:
        df (DataFrame): DataFrame with 'gene_id' and 'concentration_uM' columns
        filepath (str): Filepath formatted as 'species_bait_data.csv'

    Returns:
        float: Concentration in μM, or None if not found
    """

    try:
        _, bait = extract_species_bait(filepath)

        if 'gene_id' not in df.columns or 'concentration_uM' not in df.columns:
            print(" Required columns ('gene_id', 'concentration_uM') not found.")
            return None

        # Match bait gene name, case-insensitive
        row = df[df['gene_id'].str.upper().str.strip() == bait.strip().upper()]

        if row.empty:
            print(f"❌ Bait '{bait}' not found in gene_id column.")
            return None

        return row['concentration_uM'].values[0]

    except Exception as e:
        print(f"❌ Error in get_bait_concentration: {e}")
        return None

####################### Calculate Kapp #####################

def calculate_Kapp(df):
    """
    Calculate apparent binding constants (Kapp) for bait from pKapp values
    Parameters:
        df (DataFrame): DataFrame containing pKapp values for the bait
    Returns:
        DataFrame: Updated DataFrame with new 'Kapp' column
    """

    df = df.copy()  # Create a copy to avoid modifying the original
    
    # Calculate Kapp
    if 'pKapp' not in df.columns:
        raise ValueError("Missing 'pKapp' column in dataframe")
    df['Kapp'] = 10 ** (-df['pKapp']) * 1e6
    return df
    

#################### Complexome Calculations ##########################################

def calculate_complexome(df, filepath):
    """
    Adds a 'complexome' column to the DataFrame using bait information extracted from the filepath.
    The complexome is calculated as the concentration of the complex formed between the bait and each partner.

    Parameters:
        df (DataFrame): DataFrame with 'gene_id', 'uniprot_id', 'concentration_uM', and 'Kapp' columns
        filepath (str): Filepath from which to extract the bait name (e.g., 'human_src_data.csv')

    Returns:
        DataFrame: Updated DataFrame with new 'complexome' column
    """
    df = df.copy()

    try:
        _, bait = extract_species_bait(filepath)

        # Search bait by gene name
        if 'gene_id' not in df.columns:
            print("❌ 'gene_id' column not found.")
            return df

        bait_row = df[df['gene_id'].str.upper().str.strip() == bait.strip().upper()]
        if bait_row.empty:
            print(f"❌ Bait '{bait}' not found in 'gene_id' column.")
            return df

        bait_conc = bait_row['concentration_uM'].values[0]

        # Apply the formula directly (vectorized)
        A = bait_conc
        B = df['concentration_uM']
        Kd = df['Kapp']
        df['complexome'] = ((A + B + Kd) - np.sqrt((A + B + Kd) ** 2 - 4 * A * B)) / 2

        # Apply the formula directly (vectorized)
        #df['complexome'] = ((bait_conc + df['concentration_uM'] + df['Kapp']) -
                            #np.sqrt((bait_conc + df['concentration_uM'] + df['Kapp'])**2 - 4 * bait_conc * df['concentration_uM'])) / 2
        
    except Exception as e:
        print(f"❌ Error in add_complexome_column: {e}")

    return df
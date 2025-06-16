import pandas as pd
import numpy as np

####################### Calculate Kapp and delta_pKapp #####################

def calculate_Kapp(df, gene_id=None):

    """
    Calculate apparent binding constants (Kapp) for SRC and HCK from pKapp values,
    and compute the difference (delta_pKapp).
    
    Parameters:
        df (DataFrame): DataFrame containing pKapp values for SRC and HCK
        gene_id (str, optional): If provided, display results for this specific gene
    
    Returns:
        DataFrame: Updated DataFrame with additional columns (SRC_Kapp_sign, HCK_Kapp_sign, delta_pKapp)
    """

    df = df.copy()  # Create a copy to avoid modifying the original
    
    # Calculate Kapp for SRC and HCK
    df['SRC_Kapp_sign'] = 10 ** (-df['SRC_pKapp_sign']) * 1000000
    df['HCK_Kapp_sign'] = 10 ** (-df['HCK_pKapp_sign']) * 1000000
    
    # Calculate the difference delta_pKapp
    df['delta_pKapp'] = df['SRC_pKapp_sign'] - df['HCK_pKapp_sign']
    
    # If a gene_id is specified, display results for this gene
    if gene_id:
        gene_info = df[df['Gene_ID'] == gene_id][['Gene_ID', 'Uniprot_ID', 'SRC_pKapp_sign', 'SRC_Kapp_sign', 'HCK_pKapp_sign', 'HCK_Kapp_sign', 'delta_pKapp']]
        if not gene_info.empty:
            print(gene_info)
        else:
            print(f"Gene ID {gene_id} not found in the data.")
    
    return df


####################### Get Bait Concentration ####################

def get_protein_concentration(df, uniprot_id):

    """
    Retrieve the concentration (μM) of a protein based on its Uniprot ID.
    
    Parameters:
        df (DataFrame): DataFrame containing protein concentration data
        uniprot_id (str): Uniprot ID of the protein
    
    Returns:
        float: Concentration in μM, or None if not found
    """

    result = df[df['Uniprot_ID'] == uniprot_id]['Concentration_µM'].values
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


def add_complexome_column(df, total_bait_conc, bait_name):
    """
    Add a column to the DataFrame with complexome calculations for a specific bait.
    
    Parameters:
        df (DataFrame): DataFrame containing interaction data
        total_bait_conc (float): Total concentration of the bait protein (μM)
        bait_name (str): Name of the bait protein (e.g., 'SRC' or 'HCK')
    
    Returns:
        DataFrame: Updated DataFrame with new complexome column
    """

    df[f'{bait_name}_Complexome'] = calculate_complex_fraction(
        total_bait_conc, 
        df['Concentration_µM'], 
        df[f'{bait_name}_Kapp_sign']
    )
    return df

def calculate_complexome_difference(df, bait1, bait2):
    """
    Calculate the difference in complexome values between two baits.
    
    Parameters:
        df (DataFrame): DataFrame containing complexome data
        bait1 (str): Name of the first bait protein
        bait2 (str): Name of the second bait protein
    
    Returns:
        DataFrame: Updated DataFrame with delta_complexome column
    """

    df['delta_complexome'] = df[f'{bait1}_Complexome'] - df[f'{bait2}_Complexome']
    return df


###################### Log Scale Transformations ##########################################

def apply_log_transformation(df, columns, epsilon=1e-5):
    """
    Apply log10 transformation to specified columns with special handling for delta columns.
    
    For delta_* columns: applies sign(x) * log10(abs(x))
    For other columns: applies log10(x), replacing values ≤ 0 with epsilon
    
    Parameters:
        df (DataFrame): DataFrame containing columns to transform
        columns (list): List of column names to transform
        epsilon (float): Small value to replace zeros/negative values
    
    Returns:
        DataFrame: DataFrame with additional log-transformed columns
    """
     
    df = df.copy()
    for col in columns:
        if col not in df.columns:
            print(f"Column '{col}' not found.")
            continue

        df[col] = pd.to_numeric(df[col], errors='coerce')

        if col.startswith("delta_"):
            # For delta columns, preserve sign but log-transform the absolute value
            safe_abs = np.where(df[col].values == 0, epsilon, np.abs(df[col].values))
            signed_log = np.sign(df[col].values) * np.log10(safe_abs)
            df[f'log_{col}'] = signed_log
        else:
            # # General case: safe log10 (positive values only)
            safe_vals = np.where(df[col].values <= 0, epsilon, df[col].values)
            df[f'log_{col}'] = np.log10(safe_vals)

    return df       

###################### Calculate Complexome Ratio ############################

def calculate_intra_complexome_ratio(df, bait1, bait2, epsilon=1e-5):
    """
    Calculate the log ratio of complexome between two baits within the same species.
    
    Parameters:
        df (DataFrame): DataFrame containing log-transformed complexome data
        bait1 (str): Name of the first bait protein
        bait2 (str): Name of the second bait protein
        epsilon (float): Small value for numerical stability
    
    Returns:
        DataFrame: Updated DataFrame with intra-species complexome ratio column
    """

    df = df.copy()

    col1 = f'log_{bait1}_Complexome'
    col2 = f'log_{bait2}_Complexome'

    # Vérification des colonnes
    if col1 not in df.columns or col2 not in df.columns:
        raise ValueError(f"Missing columns: {col1} or {col2}")

    # Calcul des différences
    df[f'ratio_log_intra'] = df[col1] - df[col2]
    return df


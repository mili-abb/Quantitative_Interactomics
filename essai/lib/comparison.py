import os
from file_op import extract_species_and_bait, create_merged_dataset
import pandas as pd
import numpy as np

############################Comparison between interactomes ######################################
### D'abord merge les tableaux pour une même espèce, puis calculer les différences entre les baits

############################Calculate deltapKapp ################################
def calculate_deltapKapp(df, bait1, bait2):
    """Calculate the difference delta_pKapp""" 

    df = df.copy()  # Create a copy to avoid modifying the original
    #df['delta_pKapp'] = df[f'{bait1}_pKapp'] - df[f'{bait2}_pKapp']
    ####tableau multiindexé donc on utilise xs pour extraire les colonnes
    df[('delta', 'pKapp')] = df.xs('pKapp', level=2, axis=1)[(slice(None), bait1)] - \
                         df.xs('pKapp', level=2, axis=1)[(slice(None), bait2)].values
    return df

################### delta complexome #############################

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

#################### au sein d'une même espèce, calculer le ratio complexe entre deux baits ##########################
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

    col1 = f'{bait1}_log_complexome'
    col2 = f'{bait2}_log_complexome'

    # Vérification des colonnes
    if col1 not in df.columns or col2 not in df.columns:
        raise ValueError(f"Missing columns: {col1} or {col2}")

    # Calcul des différences
    df[f'ratio_log_intra'] = df[col1] - df[col2]
    return df
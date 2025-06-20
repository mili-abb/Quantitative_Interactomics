import pandas as pd
import numpy as np
import itertools

############################Comparison between interactomes ######################################

############################Calculate deltapKapp ################################
def calculate_deltapKapp(df):
    """
    For each species, compute delta_pKapp between the first and second bait 
    (as they appear in the MultiIndex level 1).

    The result is a new column: (species, 'BAIT1_minus_BAIT2', 'delta_pKapp')

    Parameters:
        df (DataFrame): MultiIndex DataFrame with columns (Species, Bait, Measurement)

    Returns:
        DataFrame: Updated DataFrame with delta_pKapp columns
    """
    df = df.copy()

    species_list = df.columns.get_level_values(0).unique()

    for species in species_list:
        baits = df[species].columns.get_level_values(0).unique()

        if len(baits) < 2:
            print(f"⚠️ Not enough baits for {species}, skipping...")
            continue

        bait1, bait2 = baits[:2]  # Only the first two baits

        try:
            delta_col = (species, f'{bait1}_minus_{bait2}', 'delta_pKapp')
            df[delta_col] = df[(species, bait1, 'pKapp')] - df[(species, bait2, 'pKapp')]
        except KeyError as e:
            print(f"❌ Missing pKapp for {species}/{bait1} or {bait2}: {e}")
        
    return df


################### delta complexome #############################

def calculate_complexome_difference(df):
    """
    For each species, compute delta_complexome between the first two baits.

    The result is added as a new column:
    (species, 'BAIT1_minus_BAIT2', 'delta_complexome')

    Parameters:
        df (DataFrame): MultiIndex DataFrame with (species, bait, column)

    Returns:
        DataFrame: Updated DataFrame with delta_complexome columns
    """
     
    df = df.copy()

    species_list = df.columns.get_level_values(0).unique()

    for species in species_list:
        baits = df[species].columns.get_level_values(0).unique()

        if len(baits) < 2:
            print(f"⚠️ Not enough baits for {species}, skipping...")
            continue

        bait1, bait2 = baits[:2]

        try:
            delta_col = (species, f'{bait1}_minus_{bait2}', 'delta_complexome')
            df[delta_col] = df[(species, bait1, 'complexome')] - df[(species, bait2, 'complexome')]
        except KeyError as e:
            print(f"❌ Missing complexome for {species}/{bait1} or {bait2}: {e}")

    return df

###################### Log Scale Transformations ##########################################

def apply_log_transformation(df, target_metric_keywords=("delta_pKapp", "delta_complexome"), epsilon=1e-5):
    """
    Apply log10 transformation to selected multi-index columns:
    - log10 on: 'concentration_uM', 'complexome'
    - signed log10 on: 'delta_complexome'
    - skip: delta_pKapp and unrelated metrics

    Parameters:
        df (DataFrame): MultiIndex column DataFrame
        epsilon (float): Value to avoid log(0)

    Returns:
        DataFrame: Updated DataFrame with log-transformed columns
    """

    df = df.copy()

    for col in df.columns:
        # col is a tuple like ('Human', 'SRC_minus_HCK', 'delta_pKapp',)
        if len(col) != 3:
            continue

        species, subkey, metric = col

        # Skip delta_pKapp
        if metric == "delta_pKapp":
            continue

        col_data = pd.to_numeric(df[col], errors='coerce')

        # Special case: signed log for delta_complexome
        if metric == "delta_complexome":
            safe_abs = np.where(col_data == 0, epsilon, np.abs(col_data))
            signed_log = np.sign(col_data) * np.log10(safe_abs)
            log_col = (species, subkey, "log_delta_complexome")
            df[log_col] = signed_log

        # Standard log for concentration or complexome
        elif metric in ["concentration_uM", "complexome"]:
            safe_vals = np.where(col_data <= 0, epsilon, col_data)
            log_vals = np.log10(safe_vals)
            log_col = (species, subkey, f"log_{metric}")
            df[log_col] = log_vals
        
        # Skip everything else (like pKapp, Kapp, etc.)
        else:
            continue

    return df
      

#################### au sein d'une même espèce, calculer le ratio complexe entre deux baits ##########################
###################### Calculate Complexome Ratio ############################

def calculate_intraspecies_log_complexome_ratio(df, epsilon=1e-5):
    """
     Automatically calculates the log ratio of complexome between two baits
    within each species and adds a 'ratio_log_intra' column for each species.

    Assumes MultiIndex columns: (species, bait, metric)
    Applies only if 'log_complexome' exists for both baits.

    Parameters:
        df (DataFrame): MultiIndex column DataFrame
        epsilon (float): Small value for numerical stability

    Returns:
        DataFrame: Updated DataFrame with intra-species complexome ratio columns
    """

    df = df.copy()

    # Get all species
    species_list = list({col[0] for col in df.columns if len(col) == 3})

    for species in species_list:
        # Extract all baits for this species that have log_complexome
        baits = [col[1] for col in df.columns if col[0] == species and col[2] == 'log_complexome']


        if len(baits) != 2:
            print(f"⚠️ Skipping species '{species}': expected 2 baits with 'log_complexome', found {len(baits)}")
            continue

        bait1, bait2 = baits

        col1 = (species, bait1, 'log_complexome')
        col2 = (species, bait2, 'log_complexome')
        ratio_col = (species, f'{bait1}_vs_{bait2}', 'ratio_log_complexome')

        df[ratio_col] = df[col1] - df[col2]

    return df
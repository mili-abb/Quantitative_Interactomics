
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib_venn import venn2

import seaborn as sns
import os
import subprocess
import math
from scipy.optimize import curve_fit
from scipy import stats

try:
    import uncertainties.unumpy as unp
    import uncertainties as unc
except:
    from pip._internal import main as pipmain
    pipmain(['install','uncertainties'])
    import uncertainties.unumpy as unp
    import uncertainties as unc

########################### Plot Complexome Analysis ####################################

def plot_complexome_analysis(df, species, bait_name, pKapp_col, log_complexome_col, save_fig=True):
    """
    Plot complexome analysis using a DataFrame with multi-level headers.
    
    Creates a 3x3 grid of bar charts showing binding affinity, concentration, and complexome
    for a specific bait protein, sorted by different criteria.
    
    Parameters:
        df (DataFrame): DataFrame with multi-level columns (species, column_name)
        species (str): Species to analyze ('Human' or 'Mouse')
        bait_name (str): Name of the bait protein ('SRC' or 'HCK')
        pKapp_col (str): Column name containing pKapp values
        log_complexome_col (str): Column name containing log-transformed complexome values
        save_fig (bool): Whether to save the figure to disk
    
    Returns:
        None: Displays and optionally saves the figure
    """

    # Determine column origin
    origin = 'Human' if species == 'Human' else 'Mouse'

    # Filter out rows with missing data for the complexome column
    df_filtered = df.dropna(subset=[(origin, log_complexome_col)])

    # Columns to display
    cols_to_show = [
        (origin, "Gene_ID"),
        (origin, "Uniprot_ID"),
        (origin, "log_Concentration_µM"),
        (origin, pKapp_col),
        (origin, log_complexome_col)
    ]
    print(f"\n{species} {bait_name} Complexome partners:\n")
    try:
        print(df_filtered[cols_to_show].head())
    except KeyError as e:
        print(f"Missing columns: {e}")
        return
  
    # Plot parameters
    sort_criteria = [(origin, pKapp_col), (origin, 'log_Concentration_µM'), (origin, log_complexome_col)]
    data_columns = [(origin, pKapp_col), (origin, "Concentration_µM"), (origin, f"{bait_name}_Complexome")]
    plot_titles = ['Affinity of Binding Partners', 'Concentration of Binding Partners', 'Complexome']
    y_labels = ['pKapp (µM)', '[Partner] (µM)', '[Complex] (µM)']
    colors = ['#016c59', '#67a9cf', '#bdc9e1']

    # Create figure and axes
    #fig, axs = plt.subplots(3, 3, figsize=(4, 4))

    fig, axs = plt.subplots(3, 3, figsize=(5, 4.5))
    fig.suptitle(f"{species} {bait_name} Complexome Analysis", fontsize=8, fontweight='bold')

    # Generate plots
    for i, sort_by in enumerate(sort_criteria):
        sorted_df = df_filtered.sort_values(by=sort_by, ascending=False).reset_index(drop=True)
        x_indices = np.arange(1, len(sorted_df) + 1)

        for j, (data_col, color, ylabel, title) in enumerate(zip(
                data_columns, colors, y_labels, plot_titles)):

            axs[i, j].bar(x_indices, sorted_df[data_col], color=color, width=0.6, alpha=0.7)
            axs[i, j].set_ylabel(ylabel, fontsize=6)
            axs[i, j].set_xlabel('Binding Partners', fontsize=6)
            axs[i, j].tick_params(labelsize=5)

            if i == 0:
                axs[0, j].set_title(title, fontsize=6, fontweight='bold')
            axs[i, j].grid(axis='y', linestyle='--', alpha=0.6)

            # Set log scale on Y axis except for pKapp
            if j > 0:
                axs[i, j].set_yscale('log')
            if j == 0:
                axs[i, j].set_ylim(3.5, sorted_df[(origin, pKapp_col)].max() + 0.5)

    plt.subplots_adjust(top=0.88, hspace=0.4, wspace=0.3)
    plt.tight_layout(rect=[0, 0, 1, 1])

    # Save figure if requested
    if save_fig:
        output_dir = "/home/gogllab/Desktop/Quantitative_Interactomics/output/Complexome_Plots"
        os.makedirs(output_dir, exist_ok=True)
        file_name = os.path.join(output_dir, f"{species}_{bait_name}_Complexome_Analysis.png")
        fig.savefig(file_name, dpi=300, bbox_inches='tight')
        print(f"Figure saved: {file_name}")

    plt.show()
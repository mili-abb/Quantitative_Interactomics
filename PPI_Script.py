import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import os
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



########################## Venn Diagram Visualization ####################################

#### REVOIR venn diagram human vs mouse (comptabiliser les genes ids matched qu'il y a chez human ET mouse !) #### 

def plot_venn_partners(df, output_file_prefix, bait=None, species=None, interspecies=False):
    """
    Plot a Venn diagram comparing interaction partners.

    Parameters:
    - df: pandas DataFrame containing protein interaction data with multi-indexed columns (species, feature)
    - output_file_prefix: str, prefix for saving the output PNG file
    - bait: str, either 'SRC' or 'HCK' for inter-species comparison (required if interspecies=True)
    - species: str, either 'Human' or 'Mouse' for intra-species comparison (required if interspecies=False)
    - interspecies: bool, whether to compare partners between species (True) or between baits within a species (False)
    """
    
    if interspecies:
        if bait not in ['SRC', 'HCK']:
            print("⚠️ Please specify bait='SRC' or 'HCK' for inter-species comparison.")
            return

        partners_human = set(df[df[('Human', f'{bait}_Complexome')] > 0][('Human', 'Gene_ID')])
        partners_mouse = set(df[df[('Mouse', f'{bait}_Complexome')] > 0][('Mouse', 'Gene_ID')])

        inter = partners_human & partners_mouse
        union = partners_human | partners_mouse
        jaccard = len(inter) / len(union) if union else 0

        plt.figure(figsize=(4, 4))
        venn = venn2([partners_human, partners_mouse], set_labels=('Human', 'Mouse'))

        if venn.set_labels[0]:
            venn.set_labels[0].set_fontsize(12)
            venn.set_labels[0].set_position((-0.7, 0))
        if venn.set_labels[1]:
            venn.set_labels[1].set_fontsize(12)
            venn.set_labels[1].set_position((0.7, 0))

        for text in venn.subset_labels:
            if text: text.set_fontsize(16)

        # Add textual annotations below the diagram
        plt.title(f"{bait} Partners: Human vs Mouse", fontsize=11)
        plt.text(-1.0, -0.9, f"Total Human: {len(partners_human)}", fontsize=10)
        plt.text(-1.0, -1.05, f"Total Mouse: {len(partners_mouse)}", fontsize=10)
        plt.text(-1.0, -1.20, f"Jaccard: {jaccard:.2f}", fontsize=10)

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.25)  # Add more space below
        plt.savefig(f'output/Venn_Diagrams/{output_file_prefix}.png', dpi=600)
        plt.show()

    else:
        if species not in ['Human', 'Mouse']:
            print("⚠️ Please specify species='Human' or 'Mouse' for intra-species comparison.")
            return

        partners_src = set(df[df[(species, 'SRC_Complexome')] > 0][(species, 'Gene_ID')])
        partners_hck = set(df[df[(species, 'HCK_Complexome')] > 0][(species, 'Gene_ID')])

        inter = partners_src & partners_hck
        union = partners_src | partners_hck
        jaccard = len(inter) / len(union) if union else 0

        plt.figure(figsize=(4, 4))
        
        venn = venn2([partners_src, partners_hck], set_labels=('SRC', 'HCK'))

        if venn.set_labels[0]:
            venn.set_labels[0].set_fontsize(12)
            venn.set_labels[0].set_position((-0.7, 0))
        if venn.set_labels[1]:
            venn.set_labels[1].set_fontsize(12)
            venn.set_labels[1].set_position((0.7, 0))

        for text in venn.subset_labels:
            if text: text.set_fontsize(16)

        # Add textual annotations below the diagram
        plt.title(f"{species}: SRC vs HCK", fontsize=11)
        plt.text(-1.0, -0.9, f"Total SRC: {len(partners_src)}", fontsize=10)
        plt.text(-1.0, -1.05, f"Total HCK: {len(partners_hck)}", fontsize=10)
        plt.text(-1.0, -1.20, f"Jaccard: {jaccard:.2f}", fontsize=10)

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.25)  # Add more space below
        plt.savefig(f'output/Venn_Diagrams/{output_file_prefix}.png', dpi=600)
        plt.show()


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



######################### Linear Regression Plot ################################
def analyze_and_plot_correlation(csv_file, output_file_prefix, species_x, col_x, species_y, col_y, graph_type, title, annotate_genes=False):

    """
    Process data, perform linear regression analysis, and create correlation plots.
    
    Parameters:
        csv_file (str): Path to CSV file with multi-level headers
        output_file_prefix (str): Prefix for the output image file
        species_x (str): Species for x-axis data ('Human' or 'Mouse')
        col_x (str): Column name for x-axis data
        species_y (str): Species for y-axis data ('Human' or 'Mouse')
        col_y (str): Column name for y-axis data
        title (str): Plot title
        graph_type (str): Type of data being plotted ('affinity', 'delta_pKapp', 'log_complexome', etc.)
        annotate_genes (bool): Whether to annotate points with gene IDs
    
    Returns:
        None: Displays and saves the plot
    """

 # Load data with multi-level headers
    data = pd.read_csv(csv_file, header=[0, 1])

    # Access columns with (species, column_name) tuples
    x_col = (species_x, col_x)
    y_col = (species_y, col_y)

    # Remove rows with missing values
    data_clean = data.dropna(subset=[x_col, y_col])
    if data_clean.empty:
        print(f"No valid data for correlation: {species_y} {col_y} vs {species_x} {col_x}.")
        return

    # Extract values for analysis
    x = data_clean[x_col].values
    y = data_clean[y_col].values
    n = len(y)

    # Define linear function for fitting
    def linear_function(x, a, b):
        return a * x + b

    # Perform curve fitting
    popt, pcov = curve_fit(linear_function, x, y)
    a, b = popt[0], popt[1]
    # Calculate R²
    r2 = 1.0 - (sum((y - linear_function(x, a, b)) ** 2) / ((n - 1.0) * np.var(y, ddof=1)))
    # Calculate uncertainty
    a_unc, b_unc = unc.correlated_values(popt, pcov)

    # Generate points for plotting the fit
    px = np.linspace(min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1, 1000)
    py = a_unc * px + b_unc
    nom = unp.nominal_values(py)
    std = unp.std_devs(py)

    # Calculate prediction bands    
    def calculate_prediction_bands(x, xd, yd, p, func, conf=0.95):
        alpha = 1.0 - conf
        N = xd.size
        var_n = len(p)
        q = stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
        se = np.sqrt(1. / (N - var_n) * np.sum((yd - func(xd, *p)) ** 2))
        sx = (x - xd.mean()) ** 2
        sxd = np.sum((xd - xd.mean()) ** 2)
        yp = func(x, *p)
        dy = q * se * np.sqrt(1.0 + (1.0 / N) + (sx / sxd))
        lpb, upb = yp - dy, yp + dy
        return lpb, upb

    lower_pred_band, upper_pred_band = calculate_prediction_bands(px, x, y, popt, linear_function)

    # Create the plot
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(px, nom, c='black', label='Fit')
    ax.fill_between(px, nom - 1.96 * std, nom + 1.96 * std, color='lightgray', alpha=0.2, label='95% CI')

    # Set point color based on graph type
    point_color = {
        'affinity_intra': '#016c59',
        'affinity_inter': '#016c59',
        'delta_pKapp': '#016c59',
        'log_complexome': '#67a9cf',
        'delta_complexome': '#67a9cf'
    }.get(graph_type, "#4F61C5")

    ax.plot(x, y, 'o', color=point_color, markersize=3, markeredgecolor='black', markeredgewidth=0.5)

# Annotate points with gene IDs if requested
    if annotate_genes:
        gene_col = (species_x, 'Gene_ID')
        if gene_col in data_clean.columns:
            for i, row in data_clean.iterrows():
                ax.annotate(row[gene_col], (row[x_col], row[y_col]), fontsize=4, alpha=0.7)

    # Axes config + labels
    label_map = {
        'affinity_intra': (r'p$\it{K}_{\rm app}$ (HCK)', r'p$\it{K}_{\rm app}$ (SRC)', [3.5, 6.5], [3.5, 6.5]),
        'affinity_inter': (r'p$\it{K}_{\rm app}$ (Mouse)', r'p$\it{K}_{\rm app}$ (Human)', [3.5, 6.5], [3.5, 6.5]),
        'delta_pKapp': (r'Δp$\it{K}_{\rm app}$ (Mouse)', r'Δp$\it{K}_{\rm app}$ (Human)', [-1.0, 0.5], [-1.0, 0.5]),
        'log_complexome': ('log[Complexe] (μM, HCK)', 'log[Complexe] (μM, SRC)', [-6, 0], [-6, 0]),
        'log_ratio_complexome': ('Ratio[complexe](Mouse)', 'Ratio[complexe](Human)', [-2, 1], [-2, 0])
    }
    
    if graph_type in label_map:
        xlab, ylab, xlim, ylim = label_map[graph_type]
        ax.set_xlabel(xlab, fontsize=15)
        ax.set_ylabel(ylab, fontsize=15)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if 'delta' in graph_type:
            ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
            ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

    # Insert regression statistics
    def extract_bait(col_name):
        for bait in ['SRC', 'HCK']:
            if bait in col_name:
                return bait
        return ''

    bait_x = extract_bait(col_x)
    bait_y = extract_bait(col_y)

    stat_text = (
        f"{species_y} {bait_y} vs {species_x} {bait_x}\n"
        f"$y = {a:.2f}x + {b:.2f}$\n"
        f"$R^2 = {r2:.2f}$\n"
        f"$a = {a_unc:.2uP}$\n"
        f"$b = {b_unc:.2uP}$"
    )
    ax.text(0.05, 0.95, stat_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout()

    # Save the plot
    os.makedirs("output", exist_ok=True)
    output_path = f'output/Lin_Reg_Plots/{output_file_prefix}.png'
    plt.savefig(output_path, dpi=600)
    print(f"Plot saved to: {output_path}")
    plt.show()

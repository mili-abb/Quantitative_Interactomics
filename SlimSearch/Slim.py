import pandas as pd
import os
import glob
import re
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import numpy as np
import random


############# Concatenate motif files into a single file #############

def concat_motif_files(file_paths, output_path):
    """
    Concatenates multiple TSV files containing the same columns, corrects the InstanceId 
    column to make it continuous, and writes the result to a single output file.
    
    Parameters:
        file_paths (list): List of paths to the files to concatenate.
        output_path (str): Path of the output file.
    """

    # Read all files into a list of DataFrames
    dfs = [pd.read_csv(fp, sep="\t") for fp in file_paths]
    
    # Concatenate the DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)

    # Recalculate the 'InstanceId' column continuously
    if 'InstanceId' in combined_df.columns:
        combined_df['InstanceId'] = range(1, len(combined_df) + 1)
    
    # Write the result to a TSV file
    combined_df.to_csv(output_path, sep="\t", index=False)
    print(f"Fichier combiné sauvegardé sous : {output_path}")


############## Calculate frequency of each motif in the proteome ################

def proteome_motif_frequencies(motif_files, species="Human", output_file="frequencies.txt"):
    """
    Calculates the frequency of proteins containing a motif for each given file
    and saves the results in a tabulated .txt file.

    Parameters:
        motif_files (list): List of paths to motif files (SlimSearch format).
        species (str): "Human" or "Mouse" - used to adjust proteome size.
        output_file (str): Name of the tabulated output .txt file.

    Returns:
        pd.DataFrame: Table containing the calculated frequencies.
    """

    # Proteome size by species
    proteome_sizes = {
        "Human": 20417,
        "Mouse": 17181
    }

    results = []

    for filepath in motif_files:

        # Read the file
        df = pd.read_csv(filepath, sep="\t")

        motif_name = os.path.splitext(os.path.basename(filepath))[0].replace("SLIM_H_", "").replace("SLIM_M_", "")

        # Count unique proteins
        num_proteins = df['ProteinAcc'].nunique()

        # Calculate frequency
        freq = num_proteins / proteome_sizes[species]

        results.append({
            "Motif": motif_name,
            "Unique_Protein_Count": num_proteins,
            "Frequency": round(freq, 4)
        })

    # Create final table
    result_df = pd.DataFrame(results)

    ## Save as tabulated .txt
    result_df.to_csv(output_file, sep="\t", index=False)

    print(f"Results saved to: {output_file}")
    return result_df



################# Motif frequency in SRC/HCK partners (solo, intersection, union) #####################

# Function to generate a table of SRC/HCK protein partners with their motifs
def generate_partners_motifs_table(df, motif_folder, species="Human", baits=["SRC", "HCK"], output_path="OUTPUT"):
    """
    Generates a detailed table of SRC/HCK protein partners with their motifs.
    
    Parameters:
        df (DataFrame): Dataset containing complexome data
        motif_folder (str): Path to the folder containing motif files
        species (str): "Human" or "Mouse"
        baits (list): List of bait proteins to analyze
        output_path (str): Directory to save output files
    
    Returns:
        tuple: (all_partners_info, partners_by_bait, intersection_partners)
    """
    
    os.makedirs(output_path, exist_ok=True)
    
    ## Get all motif files for the species
    motif_files = sorted([
        f for f in glob.glob(os.path.join(motif_folder, f"SLIM_{species[0]}_motif*.txt"))
        if not os.path.basename(f).split("motif")[-1].split(".")[0].startswith("3_")
    ])
    
    ## Get partners for each bait once
    partners_by_bait = {}
    uniprot_by_gene = {}  # Dictionary to store gene to UniProt mapping
    all_partners_info = {}  # Dictionary to store all partner information
    
    col_geneid = (species, "Gene_ID")
    col_uniprot = (species, "Uniprot_ID")
    
    # Get partners for each bait and store gene ID and UniProt ID mapping
    for bait in baits:
        col_complexome = (species, f"{bait}_Complexome")
        
        # Get partners as gene IDs
        partners_df = df[df[col_complexome] > 0][[col_geneid, col_uniprot]].dropna(subset=[col_geneid])
        partners_gene_ids = set(partners_df[col_geneid].astype(str).str.strip().str.upper())
        partners_by_bait[bait] = partners_gene_ids
        
        # Create gene ID to UniProt ID mapping
        for _, row in partners_df.iterrows():
            gene_id = str(row[col_geneid]).strip().upper()
            uniprot_id = str(row[col_uniprot]).strip() if pd.notna(row[col_uniprot]) else "NA"
            uniprot_by_gene[gene_id] = uniprot_id
            
            # Initialize entry in all_partners_info
            if gene_id not in all_partners_info:
                all_partners_info[gene_id] = {
                    "Uniprot_ID": uniprot_id,
                    "SRC_Partner": False,
                    "HCK_Partner": False,
                    "Intersection": False,
                    "Motifs": set()
                }
            
            # Mark as partner of this bait
            all_partners_info[gene_id][f"{bait}_Partner"] = True
    
    # Calculate intersection of partners
    intersection_partners = partners_by_bait[baits[0]]
    for bait in baits[1:]:
        intersection_partners &= partners_by_bait[bait]
    
    # Mark intersection partners
    for gene_id in intersection_partners:
        if gene_id in all_partners_info:
            all_partners_info[gene_id]["Intersection"] = True
    
    # For each motif, update partner information
    for motif_file in motif_files:
        motif_name = re.search(r"(motif\d+)", os.path.basename(motif_file)).group(1)
        motif_df = pd.read_csv(motif_file, sep="\t")
        
        # Clean GeneName to avoid matching bugs
        motif_genes = set(motif_df["GeneName"].dropna().astype(str).str.strip().str.upper())
        
        #  # Update motif information in all_partners_info
        for gene_id in motif_genes:
            if gene_id in all_partners_info:
                all_partners_info[gene_id]["Motifs"].add(motif_name)
    
    # Create comprehensive output file
    comprehensive_file = os.path.join(output_path, f"{species}_comprehensive_partners_motifs.tsv")
    with open(comprehensive_file, 'w') as f:
        f.write("Gene_ID\tUniprot_ID\tSRC_Partner\tHCK_Partner\tIntersection\tMotifs\n")
        
        for gene_id, info in sorted(all_partners_info.items()):
            motifs_str = ";".join(sorted(info["Motifs"])) if info["Motifs"] else "None"
            f.write(f"{gene_id}\t{info['Uniprot_ID']}\t{info['SRC_Partner']}\t{info['HCK_Partner']}\t{info['Intersection']}\t{motifs_str}\n")
    
    print(f"Comprehensive partners and motifs data saved to: {comprehensive_file}")
    
    return all_partners_info, partners_by_bait, intersection_partners

 

####### Utility function to get partners of a bait protein
def get_protein_partners(df, bait, species):
    """
    Retrieves the interaction partners of a specific bait protein.
    
    Parameters:
        df (DataFrame): Dataset containing complexome data
        bait (str): Bait protein name (e.g., "SRC", "HCK")
        species (str): Species name ("Human" or "Mouse")
        
    Returns:
        set: Set of gene IDs of the partners
    """

    col_complexome = (species, f"{bait}_Complexome")
    col_geneid = (species, "Gene_ID")
    return set(df[df[col_complexome] > 0][col_geneid].dropna().astype(str).str.strip().str.upper())

####### Main function  to analyze motif enrichment in SRC, HCK, union and intersection of partners

def analyze_motif_enrichment(df, motif_folder, species="Human", baits=["SRC", "HCK"],  output_path="OUTPUT", output_file=None):
    
    """
    Analyzes motif enrichment in protein interaction partners for specified baits,
    including individual baits, their union, and their intersection.
    Also saves lists of intersection partners with each motif to separate files.

    Parameters:
        df (DataFrame): Dataset containing complexome data
        motif_folder (str): Path to folder containing motif files
        species (str): "Human" or "Mouse"
        baits (list): List of bait proteins to analyze
        output_path (str): Directory to save output files
        output_file (str, optional): Name of the output file. If None, defaults to "{species[0]}_partners_freq.txt"

    Returns:
        DataFrame: Results of motif frequency analysis for partners
    """

    results = []
    os.makedirs(output_path, exist_ok=True)

    if output_file is None:
        output_file = f"{species[0]}_partners_freq.txt"

    # All motif files for the species (e.g., SLIM_H_motif1.txt → SLIM_H_motif9.txt)
    motif_files = sorted([
        f for f in glob.glob(os.path.join(motif_folder, f"SLIM_{species[0]}_motif*.txt"))
        if not os.path.basename(f).split("motif")[-1].split(".")[0].startswith("3_")
    ])
    
    # Get partners for each bait once
    partners_by_bait = {}
    for bait in baits:
        partners_by_bait[bait] = get_protein_partners(df, bait, species)
    
    
    # Calculate intersection of partners
    intersection_partners = set(partners_by_bait[baits[0]])
    for bait in baits[1:]:
        intersection_partners &= set(partners_by_bait[bait])
    
    # For each motif, analyze the partners
    for motif_file in motif_files:
        motif_name = re.search(r"(motif\d+)", os.path.basename(motif_file)).group(1)
        motif_df = pd.read_csv(motif_file, sep="\t")
        # Clean GeneName to avoid matching bugs
        motif_genes = set(motif_df["GeneName"].dropna().astype(str).str.strip().str.upper())
        
        
        # Analysis for each individual bait
        for bait in baits:
            partners = set(partners_by_bait[bait])
            total_partners = len(partners)
            with_motif = len(partners & motif_genes)
            freq = with_motif / total_partners if total_partners > 0 else 0

            results.append({
                "Motif": motif_name,
                "Species": species,
                "Bait": bait,
                "Total_partners": total_partners,
                "Partners_with_motif": with_motif,
                "Frequency": freq
            })
        
        # Calculate and analyze the UNION
        union_partners = set()
        for bait in baits:
            union_partners.update(partners_by_bait[bait])
        
        total_union = len(union_partners)
        with_motif_union = len(union_partners & motif_genes)
        freq_union = with_motif_union / total_union if total_union > 0 else 0
        
        results.append({
            "Motif": motif_name,
            "Species": species,
            "Bait": "UNION",
            "Total_partners": total_union,
            "Partners_with_motif": with_motif_union,
            "Frequency": freq_union
        })
        
        # Calculate and analyze the INTERSECTION

        inter_partners = set(intersection_partners)
        total_intersection = len(inter_partners)
        intersection_with_motif = inter_partners & motif_genes
        with_motif_intersection = len(intersection_with_motif)
        freq_intersection = with_motif_intersection / total_intersection if total_intersection > 0 else 0
            
        results.append({
            "Motif": motif_name,
            "Species": species,
            "Bait": "INTERSECTION",
            "Total_partners": total_intersection,
            "Partners_with_motif": with_motif_intersection,
            "Frequency": freq_intersection
        })
    

    # Create final table
    result_df = pd.DataFrame(results)
    
    # Full path to output file
    full_output_path = os.path.join(output_path, output_file)
    
    # Save as tabulated .txt
    result_df.to_csv(full_output_path, sep="\t", index=False)
    
    print(f"Results saved to: {full_output_path}")
    return result_df


########################## Calculate P-Values for motif frequencies ##########################

def compute_motif_enrichment_pvalues(proteome_freq_df, partners_freq_df, species="Human"):

    """
    Calculates the p-value for motif enrichment in interaction partners (SRC, HCK)
    compared to its frequency in the whole proteome using Fisher's exact test.

    Parameters:
        proteome_freq_df (DataFrame): Contains columns ['Motif', 'Unique_Protein_Count', 'Frequency']
        partners_freq_df (DataFrame): Contains columns ['Motif', 'Species', 'Bait', 'Total_Partners', 'Partners_With_Motif', 'Frequency']
        species (str): "Human" or "Mouse"

    Returns:
        DataFrame: Merged table enriched with a 'p_value' column
    """

    proteome_sizes = {
        "Human": 20417,
        "Mouse": 17181
    }

    # Ensure motif IDs are strings and matching across both dataframes before continuing
    proteome_freq_df["Motif"] = proteome_freq_df["Motif"].astype(str)
    partners_freq_df["Motif"] = partners_freq_df["Motif"].astype(str)

    # Rename columns to match the expected column names in the merge
    proteome_freq_df = proteome_freq_df.rename(columns={
        "Unique_Protein_Count": "Unique_Protein_Count",
        "Frequency": "Frequency_Proteome"
    })
    
    partners_freq_df = partners_freq_df.rename(columns={
        "Total_partners": "Total_partners",
        "Partners_with_motif": "Partners_with_motif",
        "Frequency": "Frequency_Partners"
    })
    merged = pd.merge(proteome_freq_df, partners_freq_df, on="Motif")

    all_results = []

    for _, row in merged.iterrows():
        # Fisher's exact test contingency table:
        # a = partners with motif
        # b = partners without motif
        # c = non-partners with motif
        # d = non-partners without motif
        a = row["Partners_with_motif"]
        b = row["Total_partners"] - a
        total_with_motif = row["Unique_Protein_Count"]
        c = total_with_motif - a
        d = proteome_sizes[species] - (a + b + c)

        # Prevent negative values
        d = max(0, d)

        contingency = [[a, b], [c, d]]
        _, p_value = fisher_exact(contingency, alternative="greater")

        row["p_value"] = p_value
        all_results.append(row)

    return pd.DataFrame(all_results)

############## Benjamini-Hochberg FDR correction ######################################

def apply_fdr_correction(df, p_col="p_value", alpha=0.05):
    """
    Applies Benjamini-Hochberg FDR correction to a column of p-values.
    
    Parameters:
        df (DataFrame): Contains p-values.
        p_col (str): Name of the column containing p-values.
        alpha (float): Significance threshold for FDR.

    Returns:
        DataFrame: With additional columns 'p_adj' and 'significant'
    
    """
    df_copy = df.copy()
    pvals = df_copy[p_col].values
    rejected, pvals_corrected, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
    df_copy["p_value"] = pvals
    df_copy["p_adj"] = pvals_corrected
    df_copy["significant"] = rejected
    return df_copy


############### Extract top significant motifs #######################################

def extract_top_significant_motifs(input_file, output_file, alpha=0.05):

    """
    Reads a file containing adjusted p-value (FDR) results
    and extracts significant motifs (p_adj < alpha), sorted by p_adj.

    Parameters:
        input_file (str): Path to the input file (.txt or .csv).
        output_file (str): Path to the filtered output file.
        alpha (float): Significance threshold for p_adj (default = 0.05).

    Returns:
        pd.DataFrame: DataFrame containing sorted significant rows.
    
    """

    # Read the file (expected with tab separator)
    df = pd.read_csv(input_file, sep="\t")

    # Filter significant results
    significant_motifs = df[df["p_adj"] < alpha].copy()

    # Sort by increasing significance (p_adj)
    significant_motifs = significant_motifs.sort_values(by="p_adj")

    # Save
    significant_motifs.to_csv(output_file, sep="\t", index=False)
    print(f"Significant motifs saved to : {output_file}")

    return significant_motifs

    ################ Calculate Fold Change ################

def calculate_fold_change (df):

    """
    Calculates Fold Change (FC) and log2 Fold Change (log2FC)
    between frequency in partners and frequency in the proteome.

    Parameters:
        df (pd.DataFrame): Must contain 'Frequency_Partners' and 'Frequency_Proteome'

    Returns:
        pd.DataFrame: With two additional columns 'FoldChange' and 'log2FC'
    """

    df = df.copy()
    epsilon = 1e-10  # to avoid division by 0

    # Calculate raw FC
    df["FoldChange"] = (df["Frequency_Partners"] + epsilon) / (df["Frequency_Proteome"] + epsilon)

    # Calculate log2(FC)
    df["log2FC"] = np.log2(df["FoldChange"])

    return df


############### Generate Volcano plots ################################################################

#def create_volcano_plots(df, species, output_dir=None, alpha=0.05, annotate_top=9):

    """
    Creates separate volcano plots for each group (SRC, HCK, UNION, INTERSECTION).

    Parameters:
        df (pd.DataFrame): Must contain 'Frequency_Partners', 'p_adj', 'significant', 'Motif', 'Bait'
        species (str): For plot titles
        output_dir (str): Folder to save PNG files (optional)
        alpha (float): Significance threshold (p_adj < alpha)
        annotate_top (int): Number of points to annotate (by p_value significance)
    """
    
    # List of groups to analyze
    groups = ["SRC", "HCK", "UNION", "INTERSECTION"]
    
    # For each group, create a volcano plot

    for group_name in groups:
        # Filter for current group
        group_df = df[df["Bait"] == group_name].copy()
        
        # Skip if group doesn't exist in the data
        if group_df.empty:
            print(f"No data for group {group_name}. Volcano plot not generated.")
            continue
        
        # Calculate -log10(p_adj)
        group_df["-log10(p_adj)"] = -np.log10(group_df["p_adj"].replace(0, np.nan))

        # Define colors
        colors = group_df["significant"].map({True: "red", False: "gray"})

        # Create plot
        plt.figure(figsize=(4, 4))
        plt.scatter(group_df["Frequency_Partners"], group_df["-log10(p_adj)"],
                    c=colors, edgecolor='black', s=80, alpha=0.7)


        # Title and axes
        plt.title(f"Volcano Plot - {group_name} ({species})", fontsize=14)
        plt.xlabel("Motif Frequency in Partners", fontsize=12)
        plt.ylabel("-log10(p-value adj)", fontsize=12)
        plt.legend()

        # Annotations: most extreme significant points
        significant_motifs = group_df[group_df["significant"]].copy()
        if annotate_top > 0 and not significant_motifs.empty:
            # Sort by significance (higher -log10(p_adj) means more significant)
            top_motifs = significant_motifs.sort_values("-log10(p_adj)", ascending=False).head(annotate_top)
            for _, row in top_motifs.iterrows():
                plt.text(row["Frequency_Partners"], row["-log10(p_adj)"] + 0.2,
                        f"{row['Motif']}",
                        fontsize=12, ha='center', va='bottom', clip_on=False)

        plt.tight_layout()
        
        # Save the plot if an output directory is specified
        if output_dir:
            # Ensure directory exists
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"volcano_{species}_{group_name}.png")
            plt.savefig(output_path, dpi=300)
            print(f"Volcano plot saved: {output_path}")
        
        plt.show()

def create_volcano_plots(df, species, output_dir=None, alpha=0.05, annotate_top=9):
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    groups = ["SRC", "HCK", "UNION", "INTERSECTION"]

    for group_name in groups:
        group_df = df[df["Bait"] == group_name].copy()
        if group_df.empty:
            print(f"No data for group {group_name}. Volcano plot not generated.")
            continue

        group_df["-log10(p_adj)"] = -np.log10(group_df["p_adj"].replace(0, np.nan))

        colors = group_df["significant"].map({True: "red", False: "gray"})

        plt.figure(figsize=(5, 5))
        ax = plt.gca()

        scatter = ax.scatter(group_df["Frequency_Partners"], group_df["-log10(p_adj)"],
                             c=colors, edgecolor='black', s=80, alpha=0.7, linewidth=0.5)

        # Ajuster les limites de l'axe y
        y_max = group_df["-log10(p_adj)"].max()
        ax.set_ylim(0, y_max + 1.5)

        # Titres et axes
        ax.set_title(f"{group_name} ({species})", fontsize=14, weight='bold')
        ax.set_xlabel("Motif Frequency in Partners", fontsize=12)
        ax.set_ylabel("-log10(p-value adj)", fontsize=12)
        ax.tick_params(labelsize=12)

        # Supprimer la légende vide
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(fontsize=10)

        # Annoter les motifs significatifs
        significant_motifs = group_df[group_df["significant"]].copy()
        if annotate_top > 0 and not significant_motifs.empty:
            top_motifs = significant_motifs.sort_values("-log10(p_adj)", ascending=False).head(annotate_top)
            for _, row in top_motifs.iterrows():
                x = row["Frequency_Partners"]
                y = row["-log10(p_adj)"]
                dx = -0.03 if x > 0.8 else 0.03  # décalage horizontal
                dy = 0.5  # décalage vertical
                ax.annotate(row["Motif"],
                            xy=(x, y),
                            xytext=(x + dx, y + dy),
                            textcoords='data',
                            fontsize=10,
                            ha='center',
                            va='bottom',
                            arrowprops=dict(arrowstyle='-', lw=0.4, color='gray', alpha=0.7),
                            clip_on=False)

        plt.tight_layout()

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"volcano_{species}_{group_name}.png")
            plt.savefig(output_path, dpi=300)
            print(f"Volcano plot saved: {output_path}")

        plt.show()




################ Fonction API UNIPROT récupérer la séquence des partenaires ########################

def fetch_uniprot_sequences(input_file, output_file=None, batch_size=100):
    """
    Enrichit un fichier de partenaires avec les séquences complètes des protéines
    en interrogeant l'API UniProtKB.
    
    Parameters:
        input_file (str): Chemin vers le fichier TSV contenant les données des partenaires
        output_file (str): Chemin pour le fichier de sortie (par défaut: input_file + '_with_sequences.tsv')
        batch_size (int): Nombre d'identifiants UniProt à traiter par lot (pour optimiser les requêtes API)
        
    Returns:
        pd.DataFrame: DataFrame contenant les données enrichies avec les séquences
    """
    import pandas as pd
    import requests
    import time
    import os

    # Définir le fichier de sortie si non spécifié
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_with_sequences.tsv"
    
    # Lire le fichier d'entrée
    df = pd.read_csv(input_file, sep="\t")
    print(f"Fichier chargé avec {len(df)} entrées")
    
    # Extraire les identifiants UniProt valides (pas 'NA')
    df['Uniprot_ID'] = df['Uniprot_ID'].astype(str).str.strip()
    valid_uniprots = df[df['Uniprot_ID'] != "NA"]['Uniprot_ID'].unique().tolist()
    print(f"Identifiants Uniprot valides à traiter: {len(valid_uniprots)}")
    
    # Initialiser un dictionnaire pour stocker les séquences
    sequences = {}
    
    # Traitement par lots pour éviter de surcharger l'API
    for i in range(0, len(valid_uniprots), batch_size):
        batch = valid_uniprots[i:i+batch_size]
        print(f"Traitement du lot {i//batch_size + 1}/{(len(valid_uniprots) + batch_size - 1)//batch_size}...")
        
        # Construire une requête avec OR entre les accession
        ids_query = " OR ".join([f"accession:{uid}" for uid in batch])
        url = f"https://rest.uniprot.org/uniprotkb/search?query={ids_query}&format=tsv&fields=accession,sequence"
        print(f"URL requête: {url}")  # pour test
        
        # Faire la requête avec gestion des erreurs et des tentatives
        max_retries = 3
        retry_delay = 5  # secondes
        
        for attempt in range(max_retries):
            try:
                response = requests.get(url)
                response.raise_for_status()  # Génère une exception pour les codes d'erreur HTTP
                print(response.text[:300])  # aperçu pour débogage

                # Traiter la réponse
                lines = response.text.strip().split('\n')
                if len(lines) <= 1:  # Seulement l'en-tête, pas de données
                    print("⚠️ Aucune donnée récupérée pour ce lot.")
                    break
                    
                # Sauter l'en-tête et traiter chaque ligne
                for line in lines[1:]:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        accession = parts[0].strip()
                        sequence = parts[-1].strip()
                        sequences[accession] = sequence
                
                break  # Sortir de la boucle de tentatives si la requête a réussi
                
            except requests.exceptions.RequestException as e:
                if attempt < max_retries - 1:
                    print(f"Erreur lors de la requête: {e}. Nouvelle tentative dans {retry_delay}s...")
                    time.sleep(retry_delay)
                    retry_delay *= 2  # Augmenter le délai pour les tentatives suivantes
                else:
                    print(f"Échec après {max_retries} tentatives pour le lot {i//batch_size + 1}. Erreur: {e}")
        
        # Pause pour éviter de surcharger l'API
        time.sleep(1)
    
    # Ajouter une colonne de séquences au DataFrame
    df['Sequence'] = df['Uniprot_ID'].apply(lambda x: sequences.get(x, ""))
    
    # Sauvegarder le DataFrame enrichi
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Données enrichies avec séquences sauvegardées dans: {output_file}")
    print(f"Séquences récupérées: {sum(1 for seq in df['Sequence'] if seq)}/{len(df)}")
    
    return df

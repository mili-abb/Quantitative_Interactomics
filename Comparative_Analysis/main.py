'''This script serves as the main entry point for the Comparative Analysis project.
It imports necessary modules and functions from the project, allowing for the execution of various analyses related to
complexome data, including species bait extraction, CSV file operations, calculations of Kapp and complexome,
comparison of complexomes, plotting Venn diagrams, and performing linear regression analysis.'''

__version__ = "0.1.1"  # Version of the script
__author__ = "Milissa Abboute"
__email__ = "abboutemilissapro@gmail.com"
__date__ = "2025-06-27"  # Date of the last update

    ##-- IMPORT LIBRARIES --##
import os
from lib.file_op import (
    extract_species_bait,
    read_csv_file,
    save_csv_file,
    merge_files,
    group_columns_by_species
)
from lib.calcul import (
    calculate_Kapp,
    get_bait_concentration,
    calculate_complexome,
    log_scale
)
from lib.comparison import (
    calculate_deltapKapp,
    calculate_complexome_difference,
    log_delta,
    calculate_intraspecies_log_complexome_ratio
)
from lib.Venn_diag import plot_venn_diagram
from lib.complexome_plot import plot_complexome_analysis
from lib.lin_reg import analyze_and_plot_correlation

    
def main():

    print("Script started...")

    # === Step 1: Process each input file individually ===
    input_files = {
        "human_src": "input/human_src_data.csv",
        "human_hck": "input/human_hck_data.csv",
        "mouse_src": "input/mouse_src_data.csv",
        "mouse_hck": "input/mouse_hck_data.csv"
    }

    processed_paths = []

    for key, filepath in input_files.items():
        print(f"\n📂 Processing file: {filepath}")
        output_path = f"output/{key}_data.csv"
        df_raw = read_csv_file(filepath, header=0)
        species, bait = extract_species_bait(filepath)
        print(f"➡ Species: {species}, Bait: {bait}")
        get_bait_concentration(df_raw, filepath)

        df = calculate_Kapp(df_raw)
        df = calculate_complexome(df, filepath)
        df = log_scale(df)
        save_csv_file(df, output_path)
        processed_paths.append(output_path)

    # === Step 2: Merge and calculate delta values ===
    print("\n🔀 Merging processed files...")
    merged_path = "output/merged.csv"
    df_merged = merge_files(processed_paths, merged_path)

    print("🔧 Calculating delta and ratio metrics...")
    df = calculate_deltapKapp(df_merged)
    df = calculate_complexome_difference(df)
    df = log_delta(df)
    df = calculate_intraspecies_log_complexome_ratio(df)
    df = group_columns_by_species(df)
    save_csv_file(df, merged_path)

    # === Step 3: Generate all plots ===
    print("\n📊 Generating Venn diagrams...")
    plot_venn_diagram(df, output="output/venn_diagrams")

    print("\n📈 Generating complexome plots...")
    plot_complexome_analysis(df, output="output/complexome_plots")

    print("\n📉 Generating linear regression plots...")
    analyze_and_plot_correlation(df, output_dir="output/Lin_reg_plots", annotate_genes=False)

    print("\n✅ Full analysis pipeline complete!")


if __name__ == "__main__":
    main()

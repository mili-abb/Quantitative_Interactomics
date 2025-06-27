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
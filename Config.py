"""
PPI Analysis Configuration File
===============================

This module contains all configuration parameters, constants, and default settings
for the Protein-Protein Interaction analysis pipeline.

Sections:
- File paths and directories
- Visualization parameters (colors, sizes, DPI)
- Biochemical constants
- Plot configuration
- Analysis parameters
"""

import os
from pathlib import Path

# =============================================================================
# FILE PATHS AND DIRECTORIES
# =============================================================================

# Base directories
BASE_DIR = "/home/gogllab/Desktop/Quantitative_Interactomics"
DATA_DIR = os.path.join(BASE_DIR, "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

# Output subdirectories
COMPLEXOME_PLOTS_DIR = os.path.join(OUTPUT_DIR, "Complexome_Plots")
LINEAR_REG_PLOTS_DIR = os.path.join(OUTPUT_DIR, "Lin_Reg_Plots")
VENN_DIAGRAMS_DIR = os.path.join(OUTPUT_DIR, "Venn_Diagrams")

# Default input file paths
DEFAULT_HUMAN_DATA = os.path.join(DATA_DIR, "Human.csv")
DEFAULT_MOUSE_DATA = os.path.join(DATA_DIR, "Mouse.csv")
DEFAULT_MERGED_DATA = os.path.join(OUTPUT_DIR, "Matching.csv")

# =============================================================================
# BAIT AND SPECIES CONFIGURATION
# =============================================================================
# Bait proteins and species for analysis

BAIT_PROTEINS = ['SRC', 'HCK']
SPECIES = ['Human', 'Mouse']

# =============================================================================
# VISUALIZATION PARAMETERS
# =============================================================================

# Color schemes
COLORS = {
    # Primary analysis colors
    'affinity': '#016c59',          # Dark teal for binding affinity
    'concentration': '#67a9cf',      # Blue for concentration
    'complexome': '#bdc9e1',        # Light blue for complexome
    
    # Graph type specific colors
    'affinity_intra': '#016c59',
    'affinity_inter': '#016c59', 
    'delta_pKapp': '#016c59',
    'log_complexome': '#67a9cf',
    'delta_complexome': '#67a9cf',
    
    # Default and utility colors
    'default_point': '#4F61C5',
    'fit_line': 'black',
    'confidence_band': 'lightgray',
    'grid_lines': 'gray',
    'reference_lines': 'gray'
}

# Color palette for multiple series
COLOR_PALETTE = ['#016c59', '#67a9cf', '#bdc9e1', '#4F61C5', '#8B0000', '#228B22']


# =============================================================================
# PLOT CONFIGURATION
# =============================================================================

# Figure sizes (width, height in inches)
FIGURE_SIZES = {
    'correlation_plot': (4, 4),
    'complexome_analysis': (5, 4.5),
    'venn_diagram': (4, 4),
    'default': (6, 4)
}

# DPI settings for different output types
DPI_SETTINGS = {
    'screen': 100,
    'print': 300,
    'publication': 600,
    'default': 300
}

# Font sizes
FONT_SIZES = {
    'title': 11,
    'suptitle': 8,
    'subtitle': 6,
    'axis_label': 15,
    'axis_label_small': 6,
    'tick_label': 5,
    'legend': 9,
    'annotation': 4,
    'venn_label': 12,
    'venn_number': 16,
    'stats_text': 9,
    'info_text': 10
}

# Line and marker properties
LINE_STYLES = {
    'fit_line': '-',
    'grid': '--',
    'reference': '--'
}

MARKER_PROPERTIES = {
    'size': 3,
    'edge_color': 'black',
    'edge_width': 0.5,
    'alpha': 0.7
}

# Grid and layout
GRID_PROPERTIES = {
    'alpha': 0.6,
    'axis': 'y',
    'linestyle': '--'
}

LAYOUT_PROPERTIES = {
    'tight_layout': True,
    'subplots_adjust': {
        'top': 0.88,
        'hspace': 0.4,
        'wspace': 0.3,
        'bottom': 0.25  # For Venn diagrams
    }
}


# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# Statistical analysis
CONFIDENCE_LEVEL = 0.95
STATISTICAL_ALPHA = 0.05

# Linear regression parameters
LINEAR_REG_PARAMS = {
    'confidence_interval': 0.95,
    'prediction_bands': True,
    'show_uncertainty': True
}

# Complexome analysis parameters
COMPLEXOME_PARAMS = {
    'min_complex_threshold': 0,  # Minimum complexome value to consider
    'log_transform': True,
    'epsilon': LOG_EPSILON
}

# =============================================================================
# PLOT AXIS LIMITS AND LABELS
# =============================================================================

# Predefined axis limits for different plot types
AXIS_LIMITS = {
    'affinity_intra': {'x': [3.5, 6.5], 'y': [3.5, 6.5]},
    'affinity_inter': {'x': [3.5, 6.5], 'y': [3.5, 6.5]},
    'delta_pKapp': {'x': [-1.0, 0.5], 'y': [-1.0, 0.5]},
    'log_complexome': {'x': [-6, 0], 'y': [-6, 0]},
    'log_ratio_complexome': {'x': [-2, 1], 'y': [-2, 0]}
}

# Axis labels with proper formatting
AXIS_LABELS = {
    'pKapp_hck': r'p$\it{K}_{\rm app}$ (HCK)',
    'pKapp_src': r'p$\it{K}_{\rm app}$ (SRC)',
    'pKapp_mouse': r'p$\it{K}_{\rm app}$ (Mouse)',
    'pKapp_human': r'p$\it{K}_{\rm app}$ (Human)',
    'delta_pKapp_mouse': r'Δp$\it{K}_{\rm app}$ (Mouse)',
    'delta_pKapp_human': r'Δp$\it{K}_{\rm app}$ (Human)',
    'log_complex_hck': 'log[Complexe] (μM, HCK)',
    'log_complex_src': 'log[Complexe] (μM, SRC)',
    'ratio_complex_mouse': 'Ratio[complexe](Mouse)',
    'ratio_complex_human': 'Ratio[complexe](Human)',
    'concentration_um': '[Partner] (µM)',
    'complex_concentration': '[Complex] (µM)',
    'binding_partners': 'Binding Partners'
}

# Plot titles
PLOT_TITLES = {
    'affinity_partners': 'Affinity of Binding Partners',
    'concentration_partners': 'Concentration of Binding Partners',
    'complexome': 'Complexome',
    'complexome_analysis': '{} {} Complexome Analysis',
    'venn_intra': '{}: SRC vs HCK',
    'venn_inter': '{} Partners: Human vs Mouse'
}

# =============================================================================
# FILE FORMAT SETTINGS
# =============================================================================

# Supported file formats
SUPPORTED_FORMATS = ['png', 'pdf', 'svg', 'eps']
DEFAULT_FORMAT = 'png'

# CSV settings
CSV_SETTINGS = {
    'separator': ',',
    'index': False,
    'encoding': 'utf-8'
}

# =============================================================================
# COLUMN MAPPING AND NAMING CONVENTIONS
# =============================================================================

# Standard column names used throughout the pipeline
STANDARD_COLUMNS = {
    'gene_id': 'Gene_ID',
    'uniprot_id': 'Uniprot_ID',
    'concentration': 'Concentration_µM',
    'log_concentration': 'log_Concentration_µM',
    'src_pkapp': 'SRC_pKapp_sign',
    'hck_pkapp': 'HCK_pKapp_sign',
    'src_kapp': 'SRC_Kapp_sign',
    'hck_kapp': 'HCK_Kapp_sign',
    'delta_pkapp': 'delta_pKapp',
    'src_complexome': 'SRC_Complexome',
    'hck_complexome': 'HCK_Complexome',
    'delta_complexome': 'delta_complexome',
    'ratio_log_intra': 'ratio_log_intra'
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def create_output_directories():
    """Create all necessary output directories if they don't exist."""
    directories = [
        OUTPUT_DIR,
        COMPLEXOME_PLOTS_DIR,
        LINEAR_REG_PLOTS_DIR,
        VENN_DIAGRAMS_DIR
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        
def get_color(graph_type, default='default_point'):
    """Get color for a specific graph type."""
    return COLORS.get(graph_type, COLORS.get(default, '#4F61C5'))

def get_figure_size(plot_type, default='default'):
    """Get figure size for a specific plot type."""
    return FIGURE_SIZES.get(plot_type, FIGURE_SIZES.get(default, (6, 4)))

def get_dpi(output_type='default'):
    """Get DPI setting for output type."""
    return DPI_SETTINGS.get(output_type, DPI_SETTINGS['default'])

def get_axis_limits(graph_type):
    """Get axis limits for a specific graph type."""
    return AXIS_LIMITS.get(graph_type, {'x': None, 'y': None})

def get_axis_labels(label_type):
    """Get formatted axis label."""
    return AXIS_LABELS.get(label_type, label_type)

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

def validate_species(species):
    """Validate species name."""
    if species not in SPECIES:
        raise ValueError(f"Species must be one of {SPECIES}, got {species}")
    return True

def validate_bait(bait):
    """Validate bait protein name."""
    if bait not in BAIT_PROTEINS:
        raise ValueError(f"Bait must be one of {BAIT_PROTEINS}, got {bait}")
    return True

def validate_file_format(format_type):
    """Validate file format."""
    if format_type not in SUPPORTED_FORMATS:
        raise ValueError(f"Format must be one of {SUPPORTED_FORMATS}, got {format_type}")
    return True

# =============================================================================
# CONFIGURATION SUMMARY
# =============================================================================

def print_config_summary():
    """Print a summary of current configuration."""
    print("PPI Analysis Configuration Summary")
    print("=" * 50)
    print(f"Base directory: {BASE_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Supported species: {', '.join(SPECIES)}")
    print(f"Bait proteins: {', '.join(BAIT_PROTEINS)}")
    print(f"Default DPI: {get_dpi()}")
    print(f"Default format: {DEFAULT_FORMAT}")
    print(f"Log epsilon: {LOG_EPSILON}")
    print("=" * 50)

if __name__ == "__main__":
    # Test configuration
    print_config_summary()
    create_output_directories()
    print("Configuration loaded successfully!")
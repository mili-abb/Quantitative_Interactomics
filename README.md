# Protein - Protein Interaction Analysis pipeline 

This Python-based pipeline was developed to analyze data from **Native Hold-Up (nHU)** experiments, particularly focusing on the quantitative interactomes of the **SRC** and **HCK** kinases across human and mouse samples. 
It includes modules for data preprocessing, affinity and complexome calculation, comparative analyses between baits and species, and rich visualization capabilities.

## 📂 Project Structure
```.
├── PPI_Script.py               # Main script with all analysis functions
├── SRC_HCK_interactome_analysis.ipynb  # Example notebook showcasing usage
├── data/                       # Your input CSV files
├── output/
│   ├── Venn_Diagrams/          
│   ├── Complexome_Plots/       # Bar plots of binding & complexome data
│   └── Lin_Reg_Plots/          # Correlation plots (linear regression)
└── README.md                   # This documentation
```

## 🧪 Requirements

You need Python 3.12.2 and the following packages:
- pandas
- numpy
- matplotlib
- matplotlib-venn
- scipy
- uncertainties

To install dependencies:
```pip install -r requirements.txt```
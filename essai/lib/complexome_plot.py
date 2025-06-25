import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

def plot_complexome_analysis(df, output, pKapp_metric="pKapp", complexome_metric="log_complexome", save_fig=True):
    """
    Automatically plots complexome analyses for each species and bait in a multi-indexed DataFrame.
    
    Parameters:
        df (DataFrame): DataFrame with MultiIndex columns (species, bait, metric)
        output : Directory where plots will be saved
        pKapp_metric (str): Name of the metric column for affinity (default: 'pKapp')
        complexome_metric (str): Name of the log-transformed complexome column

    Returns:
        None
    """
    os.makedirs(output, exist_ok=True)
    species_list = sorted({col[0] for col in df.columns if len(col) == 3})
    
    for species in species_list:
        baits = sorted({col[1] for col in df.columns if col[0] == species and col[2] == complexome_metric})
        
        for bait in baits:
            try:
                print(f"\n🔍 Plotting: {species} - {bait} Complexome Analysis")
                
                # Filter valid rows
                df_filtered = df.dropna(subset=[(species, bait, complexome_metric)])

                if df_filtered.empty:
                    print(f"⚠️ No data available for {species} {bait}, skipping...")
                    continue

                # Plot setup
                fig, axs = plt.subplots(3, 3, figsize=(5, 4.5))
                fig.suptitle(f"{species} {bait} Complexome Analysis", fontsize=9, fontweight='bold')

                sort_by_metrics = [pKapp_metric, 'log_concentration_uM', complexome_metric]
                y_axes = [pKapp_metric, 'concentration_uM', "complexome"]
                titles = ['Affinity of Binding Partners', 'Concentration of Binding Partners', 'Complexome']
                y_labels = ['pKapp (µM)', '[Partner] (µM)', '[Complex] (µM)']
                colors = ['#016c59', '#67a9cf', '#bdc9e1']

                for i, sort_key in enumerate(sort_by_metrics):
                    sort_col = (species, bait, sort_key)
                    sorted_df = df_filtered.sort_values(by=sort_col, ascending=False)
                    x = np.arange(1, len(sorted_df) + 1)

                    for j, (y_key, color, ylabel, title) in enumerate(zip(y_axes, colors, y_labels, titles)):
                        y_col = (species, bait, y_key)
                        axs[i, j].bar(x, sorted_df[y_col], color=color, width=0.6, alpha=0.7)
                        axs[i, j].set_ylabel(ylabel, fontsize=6)
                        axs[i, j].set_xlabel('Binding Partners', fontsize=6)
                        axs[i, j].tick_params(labelsize=5)

                        if i == 0:
                            axs[0, j].set_title(title, fontsize=6, fontweight='bold')
                        axs[i, j].grid(axis='y', linestyle='--', alpha=0.6)


                        if j > 0:
                            axs[i, j].set_yscale('log')
                        if j == 0:
                            axs[i, j].set_ylim(3.5, sorted_df[(species, bait, pKapp_metric)].max() + 0.5)
                            axs[i, j].yaxis.set_major_locator(ticker.MultipleLocator(0.5))

                plt.subplots_adjust(top=0.88, hspace=0.4, wspace=0.3)
                plt.tight_layout(rect=[0, 0, 1, 1])


                if save_fig:
                    file_name = os.path.join(output, f"{species}_{bait}_Complexome_Analysis.png")
                    fig.savefig(file_name, dpi=300, bbox_inches='tight')
                    print(f"✅ Figure saved: {file_name}")

                    plt.show()
                    plt.close(fig)

            except Exception as e:
                print(f"❌ Error plotting {species} {bait}: {e}")
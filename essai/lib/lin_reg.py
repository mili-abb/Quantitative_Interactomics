import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy import stats
import uncertainties.unumpy as unp
import uncertainties as unc

def analyze_and_plot_correlation(df, output_dir="output/Lin_Reg_Plots", target_metrics=("pKapp", "delta_pKapp", "log_complexome", "delta_complexome", "ratio_log_complexome"), annotate_genes=False):
    """
    Automatise la génération de régressions linéaires inter-espèces pour chaque bait et métrique donnée.

    Paramètres :
        df (DataFrame): DataFrame à colonnes MultiIndex (species, bait, metric)
        output_dir (str): Répertoire de sauvegarde
        target_metrics (tuple): Métriques ciblées pour les corrélations
        annotate_genes (bool): Annoter ou non avec les Gene_ID
    """
    os.makedirs(output_dir, exist_ok=True)

    # Identifier espèces et baits
    species_list = sorted({col[0] for col in df.columns if len(col) == 3})
    bait_list = sorted({col[1] for col in df.columns if len(col) == 3})
    metric_list = sorted({col[2] for col in df.columns if len(col) == 3})

    label_map = {
        "pKapp": (r"p$\it{K}_{\rm app}$", [3.5, 6.5]),
        "delta_pKapp": (r"Δp$\it{K}_{\rm app}$", [-1, 0.5]),
        "log_complexome": (r"log[Complexe] (μM)", [-6, 0]),
        "delta_complexome": (r"Δlog[Complexe] (μM)", [-1, 1]),
        "ratio_log_complexome": (r"Ratio log[Complexe]", [-2.5, 1.5])
    }

    color_map = {
        'pKapp': '#016c59',
        'delta_pKapp': '#016c59',
        'log_complexome': '#67a9cf',
        'delta_complexome': '#67a9cf',
        'ratio_log_complexome': '#bdc9e1'
    }

    for bait in bait_list:
        for metric in metric_list:
            if metric not in target_metrics:
                continue

            valid_species = [sp for sp in species_list if (sp, bait, metric) in df.columns]
            if len(valid_species) < 2:
                continue

            for i in range(len(valid_species)):
                for j in range(i + 1, len(valid_species)):
                    sp_x, sp_y = valid_species[i], valid_species[j]
                    x_col = (sp_x, bait, metric)
                    y_col = (sp_y, bait, metric)

                    data_clean = df.dropna(subset=[x_col, y_col])
                    if data_clean.empty:
                        print(f"⚠️ Pas de données valides pour {bait} - {metric} : {sp_x} vs {sp_y}")
                        continue

                    x = data_clean[x_col].values
                    y = data_clean[y_col].values
                    n = len(y)

                    # Régression linéaire
                    def linear_function(x, a, b):
                        return a * x + b

                    try:
                        popt, pcov = curve_fit(linear_function, x, y)
                        a, b = popt[0], popt[1]

                        r2 = 1.0 - (sum((y - linear_function(x, a, b))**2) / ((n - 1.0) * np.var(y, ddof=1)))
                        a_unc, b_unc = unc.correlated_values(popt, pcov)

                        # Prédiction
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

                        lpb , upb = calculate_prediction_bands(px, x, y, popt, linear_function)

                        # Create the plot
                        fig, ax = plt.subplots(figsize=(4, 4))
                        ax.plot(px, nom, c='black', label='Fit')
                        ax.fill_between(px, nom - 1.96 * std, nom + 1.96 * std, color='lightgray', alpha=0.2, label='95% CI')
                        ax.plot(x, y, 'o', color=color_map.get(metric, "#4F61C5"), markersize=3, markeredgecolor='black', markeredgewidth=0.5)

                        # Annotations
                        if annotate_genes and (sp_x, bait, 'gene_id') in data_clean.columns:
                            for _, row in data_clean.iterrows():
                                ax.annotate(row[(sp_x, bait, 'gene_id')], (row[x_col], row[y_col]), fontsize=4)

                        xlab = f"{label_map[metric][0]} ({sp_x})"
                        ylab = f"{label_map[metric][0]} ({sp_y})"
                        xlim, ylim = label_map[metric][1], label_map[metric][1]

                        ax.set_xlabel(xlab, fontsize=13)
                        ax.set_ylabel(ylab, fontsize=13)
                        ax.set_xlim(xlim)
                        ax.set_ylim(ylim)

                        if "delta" in metric:
                            ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
                            ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

                        stat_text = (
                            f"$y = {a:.2f}x + {b:.2f}$\n"
                            f"$R^2 = {r2:.2f}$\n"
                            f"$a = {a_unc:.2uP}$\n"
                            f"$b = {b_unc:.2uP}$"
                        )
                        ax.text(0.05, 0.95, stat_text, transform=ax.transAxes,
                                fontsize=9, verticalalignment='top',
                                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

                        plt.tight_layout()
                        filename = f"{output_dir}/{bait}_{metric}_{sp_x}_vs_{sp_y}.png"
                        plt.savefig(filename, dpi=600)
                        plt.show()
                        plt.close()
                        print(f"✅ Plot enregistré : {filename}")

                    except Exception as e:
                        print(f"❌ Erreur pour {bait}, {metric}, {sp_x} vs {sp_y} : {e}")
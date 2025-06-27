import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy import stats
import uncertainties.unumpy as unp
import uncertainties as unc


def linear_function(x, a, b):
    return a * x + b

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
    return yp - dy, yp + dy

def perform_regression(x, y):
    popt, pcov = curve_fit(linear_function, x, y)
    a, b = popt
    a_unc, b_unc = unc.correlated_values(popt, pcov)
    r2 = 1.0 - (sum((y - linear_function(x, a, b)) ** 2) / ((len(y) - 1.0) * np.var(y, ddof=1)))
    return a, b, a_unc, b_unc, r2, popt

def generate_plot(ax, x, y, px, py_nom, std, color, label_x, label_y, stat_text, xlim, ylim):
    ax.plot(px, py_nom, c='black', label='Fit')
    ax.fill_between(px, py_nom - 1.96 * std, py_nom + 1.96 * std, color='lightgray', alpha=0.2, label='95% CI')
    ax.plot(x, y, 'o', color=color, markersize=3, markeredgecolor='black', markeredgewidth=0.5)
    ax.set_xlabel(label_x, fontsize=13)
    ax.set_ylabel(label_y, fontsize=13)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.text(0.03, 0.95, stat_text, transform=ax.transAxes, fontsize=7,
            verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.6, pad=0.3))

def analyze_and_plot_correlation(df, output_dir="output/Lin_Reg_Plots",
                                  target_metrics=("pKapp", "delta_pKapp", "log_complexome", "ratio_log_complexome"),
                                  annotate_genes=False):

    os.makedirs(output_dir, exist_ok=True)

    species_list = sorted({col[0] for col in df.columns if len(col) == 3})
    bait_list = sorted({col[1] for col in df.columns if len(col) == 3})
    target_metrics = list(target_metrics)

    color_map = {
        'pKapp': '#016c59',
        'delta_pKapp': '#016c59',
        'log_complexome': '#67a9cf',
        'ratio_log_complexome': '#67a9cf'
    }

    label_map = {
        'pKapp': (r'p$\it{K}_{\rm app}$', [3.5, 6.5], [3.5, 6.5]),
        'delta_pKapp': (r'\u0394p$\it{K}_{\rm app}$', [-1.0, 0.5], [-1.0, 0.5]),
        'log_complexome': ('log[Complex]', [-6, 0], [-6, 0]),
        'ratio_log_complexome': ('Complex ratio', [-2, 1], [-2, 0])
    }

    ordered_metrics = [m for m in target_metrics if 'pKapp' in m] + \
                      [m for m in target_metrics if 'complexome' in m]

    for metric in ordered_metrics:
        # Intra-species
        for species in species_list:
            valid_baits = [bait for bait in bait_list if (species, bait, metric) in df.columns]
            for i in range(len(valid_baits)):
                for j in range(i + 1, len(valid_baits)):
                    bait_x, bait_y = valid_baits[i], valid_baits[j]
                    x_col, y_col = (species, bait_x, metric), (species, bait_y, metric)
                    if x_col not in df.columns or y_col not in df.columns:
                        continue
                    data_clean = df.dropna(subset=[x_col, y_col])
                    if data_clean.empty:
                        continue
                    x, y = data_clean[x_col].values, data_clean[y_col].values
                    try:
                        a, b, a_unc, b_unc, r2, popt = perform_regression(x, y)
                        px = np.linspace(min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1, 1000)
                        py = a_unc * px + b_unc
                        nom = unp.nominal_values(py)
                        std = unp.std_devs(py)

                        fig, ax = plt.subplots(figsize=(4, 4))
                        generate_plot(ax, x, y, px, nom, std, color_map.get(metric, "#4F61C5"),
                                      f"{label_map[metric][0]} ({bait_x})", f"{label_map[metric][0]} ({bait_y})",
                                      f"{species}: {bait_y} vs {bait_x}\n$y = {a:.2f}x + {b:.2f}$\n$R^2 = {r2:.2f}$\n"
                                      f"$a = {a_unc:.2uP}$\n$b = {b_unc:.2uP}$",
                                      *label_map[metric][1:])

                        if "delta" in metric:
                            ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
                            ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

                        if annotate_genes and (species, bait_x, 'gene_id') in data_clean.columns:
                            for _, row in data_clean.iterrows():
                                ax.annotate(row[(species, bait_x, 'gene_id')], (row[x_col], row[y_col]), fontsize=4)

                        plt.tight_layout()
                        filename = f"{output_dir}/{metric}_intra_{species}_{bait_y}_vs_{bait_x}.png"
                        plt.savefig(filename, dpi=600)
                        plt.show()
                        plt.close()
                    except Exception as e:
                        print(f"❌ Error for {species}, {metric}, {bait_x} vs {bait_y} : {e}")

        # Inter-species
        for bait in bait_list:
            valid_species = [sp for sp in species_list if (sp, bait, metric) in df.columns]
            for i in range(len(valid_species)):
                for j in range(i + 1, len(valid_species)):
                    sp_x, sp_y = valid_species[i], valid_species[j]
                    x_col, y_col = (sp_y, bait, metric), (sp_x, bait, metric)
                    if x_col not in df.columns or y_col not in df.columns:
                        continue
                    data_clean = df.dropna(subset=[x_col, y_col])
                    if data_clean.empty:
                        continue
                    x, y = data_clean[x_col].values, data_clean[y_col].values
                    try:
                        a, b, a_unc, b_unc, r2, popt = perform_regression(x, y)
                        px = np.linspace(min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1, 1000)
                        py = a_unc * px + b_unc
                        nom = unp.nominal_values(py)
                        std = unp.std_devs(py)
                        lpb, upb = calculate_prediction_bands(px, x, y, popt, linear_function)

                        fig, ax = plt.subplots(figsize=(4, 4))
                        generate_plot(ax, x, y, px, nom, std, color_map.get(metric, "#4F61C5"),
                                      f"{label_map[metric][0]} ({sp_y})", f"{label_map[metric][0]} ({sp_x})",
                                      f"{sp_x} vs {sp_y} - {bait}\n$y = {a:.2f}x + {b:.2f}$\n$R^2 = {r2:.2f}$\n"
                                      f"$a = {a_unc:.2uP}$\n$b = {b_unc:.2uP}$",
                                      *label_map[metric][1:])

                        if "delta" in metric:
                            ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
                            ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
                            
                        if annotate_genes and (species, bait_x, 'gene_id') in data_clean.columns:
                            for _, row in data_clean.iterrows():
                                ax.annotate(row[(species, bait_x, 'gene_id')], (row[x_col], row[y_col]), fontsize=4)

                        plt.tight_layout()
                        filename = f"{output_dir}/{bait}_{metric}_{sp_x}_vs_{sp_y}.png"
                        plt.savefig(filename, dpi=600)
                        plt.show()
                        plt.close()
                    except Exception as e:
                        print(f"❌ Error for {bait}, {metric}, {sp_x} vs {sp_y} : {e}")

    print("✅ All plots generated and saved successfully.")

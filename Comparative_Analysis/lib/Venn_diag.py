import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

########################## Venn Diagram Visualization ####################################

#### REVOIR venn diagram human vs mouse (comptabiliser les genes ids matched qu'il y a chez human ET mouse !) #### 

def plot_venn_diagram (df, output):
    """
    Automatically generates Venn diagrams:
    - Intra-species: comparing baits within each species
    - Inter-species: comparing species for each bait

    Parameters:
        df (DataFrame): MultiIndex column DataFrame (species, bait, metric)
        output_dir (str): Directory to save figures
    """
    os.makedirs(output, exist_ok=True)

    # Extract all species and baits from MultiIndex
    species_list = sorted({col[0] for col in df.columns if len(col) == 3})
    bait_list = sorted({col[1] for col in df.columns if len(col) == 3 and col[2] == "complexome"})

    # Intra-species Venn (SRC vs HCK or any bait pair)
    for species in species_list:
        baits_in_species = sorted({col[1] for col in df.columns if col[0] == species and col[2] == "complexome"})
        if len(baits_in_species) < 2:
            continue

        for i in range(len(baits_in_species)):
            for j in range(i+1, len(baits_in_species)):
                bait1, bait2 = baits_in_species[i], baits_in_species[j]

                try:
                    partners_1 = set(df[df[(species, bait1, 'complexome')] > 0][(species, bait1, 'gene_id')])
                    partners_2 = set(df[df[(species, bait2, 'complexome')] > 0][(species, bait2, 'gene_id')])
                    inter = partners_1 & partners_2
                    union = partners_1 | partners_2
                    jaccard = len(inter) / len(union) if union else 0

                    plt.figure(figsize=(4, 4))
                    venn = venn2([partners_1, partners_2], set_labels=(bait1, bait2))
                    for text in venn.set_labels:
                        if text: text.set_fontsize(12)
                        if venn.set_labels[0]:
                            venn.set_labels[0].set_position((-0.7, 0))
                        if venn.set_labels[1]:
                            venn.set_labels[1].set_position((0.7, 0))
                    for text in venn.subset_labels:
                        if text: text.set_fontsize(16)

                    plt.title(f"{species}: {bait1} vs {bait2}", fontsize=11)
                    plt.text(-1.0, -0.9, f"Total {bait1}: {len(partners_1)}", fontsize=10)
                    plt.text(-1.0, -1.05, f"Total {bait2}: {len(partners_2)}", fontsize=10)
                    plt.text(-1.0, -1.20, f"Jaccard: {jaccard:.2f}", fontsize=10)
                    plt.tight_layout()
                    plt.subplots_adjust(bottom=0.25)

                    outfile = os.path.join(output, f"{species}_{bait1}_vs_{bait2}.png")
                    plt.savefig(outfile, dpi=600)
                    plt.show()
                    plt.close()
                    print(f"✅ Saved intra-species Venn: {outfile}")
                except Exception as e:
                    print(f"❌ Error for {species} {bait1} vs {bait2}: {e}")

    # Inter-species Venn (Human vs Mouse or others) for each bait
    for bait in bait_list:
        bait_species = [sp for sp in species_list if (sp, bait, "complexome") in df.columns]
        if len(bait_species) < 2:
            continue

        for i in range(len(bait_species)):
            for j in range(i+1, len(bait_species)):
                sp1, sp2 = bait_species[i], bait_species[j]

                try:
                    partners_1 = set(df[df[(sp1, bait, 'complexome')] > 0][(sp1, bait, 'gene_id')])
                    partners_2 = set(df[df[(sp2, bait, 'complexome')] > 0][(sp2, bait, 'gene_id')])
                    inter = partners_1 & partners_2
                    union = partners_1 | partners_2
                    jaccard = len(inter) / len(union) if union else 0

                    plt.figure(figsize=(4, 4))
                    venn = venn2([partners_1, partners_2], set_labels=(sp1, sp2))
                    for text in venn.set_labels:
                        if text: text.set_fontsize(12)
                        if venn.set_labels[0]:
                            venn.set_labels[0].set_position((-0.7, 0))
                        if venn.set_labels[1]:
                            venn.set_labels[1].set_position((0.7, 0))
                    for text in venn.subset_labels:
                        if text: text.set_fontsize(16)

                    plt.title(f"{bait}: {sp1} vs {sp2}", fontsize=11)
                    plt.text(-1.0, -0.9, f"{sp1}: {len(partners_1)}", fontsize=10)
                    plt.text(-1.0, -1.05, f"{sp2}: {len(partners_2)}", fontsize=10)
                    plt.text(-1.0, -1.20, f"Jaccard: {jaccard:.2f}", fontsize=10)
                    plt.tight_layout()
                    plt.subplots_adjust(bottom=0.25)

                    outfile = os.path.join(output, f"{bait}_{sp1}_vs_{sp2}.png")
                    plt.savefig(outfile, dpi=600)
                    plt.show()
                    plt.close()
                    print(f"✅ Saved inter-species Venn: {outfile}")
                except Exception as e:
                    print(f"❌ Error for {bait} {sp1} vs {sp2}: {e}")

from matplotlib_venn import venn2
import matplotlib.pyplot as plt


    
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
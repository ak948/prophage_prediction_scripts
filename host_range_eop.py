"""
Read eop.tsv file
Row 1 = Phages, Column 1 = Strains
Use seaborn to produce heatmap with log scale
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors

eop_data = pd.read_csv(f"~/phylo_eop.tsv", sep="\t", index_col=0, header=0)

sns.set_theme(style="white")
plt.figure(figsize=(12, 6))
ax = sns.heatmap(
    eop_data,
    cmap="magma_r", fmt=".2%",
    linewidths=.4, linecolor="black",
    norm=colors.LogNorm(),
    cbar_kws={"label": "EOP"}
)
# Take isolation_strain_eop.tsv and do ax.text(
isolation_strain_data = pd.read_csv(f"~/isolation_strain_eop.tsv", sep="\t")
phage_labels = list(eop_data.columns)
strain_labels = list(eop_data.index)

for _, row in isolation_strain_data.iterrows():
    phage = row["Phage"]
    strain = row["Isolation Strain"]

    if phage in phage_labels and strain in strain_labels:
        x = phage_labels.index(phage) + 0.5
        y = strain_labels.index(strain) +0.5

        ax.text(x, y, "I", horizontalalignment='center', verticalalignment='center', color='white', fontsize=8)
plt.xlabel("Phages")
plt.ylabel("Strains")
plt.title("Phage Host Ranges")
plt.tight_layout()
plt.savefig("flights_heatmap_default.svg",format="svg", dpi=300)
#plt.show()
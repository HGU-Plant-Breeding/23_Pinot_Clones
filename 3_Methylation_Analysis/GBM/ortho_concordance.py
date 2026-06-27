#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 1. Load Data
print("Loading data...")
df_ortho = pd.read_csv("one_to_one_orthologs.tsv", sep="\t")
df_class = pd.read_csv("20-13_classification.tsv", sep="\t")

# 2. Create mapping dictionaries
class_dict = dict(zip(df_class['Gene_ID'], df_class['Classification']))
cg_frac_dict = dict(zip(df_class['Gene_ID'], pd.to_numeric(df_class['CG_frac'], errors='coerce')))
chg_frac_dict = dict(zip(df_class['Gene_ID'], pd.to_numeric(df_class['CHG_frac'], errors='coerce')))

# 3. Map values to HapA and HapB
df_ortho['HapA_Class'] = df_ortho['HapA_GeneID'].map(class_dict).fillna('Missing')
df_ortho['HapB_Class'] = df_ortho['HapB_GeneID'].map(class_dict).fillna('Missing')

df_ortho['HapA_CG_frac'] = df_ortho['HapA_GeneID'].map(cg_frac_dict)
df_ortho['HapB_CG_frac'] = df_ortho['HapB_GeneID'].map(cg_frac_dict)

df_ortho['HapA_CHG_frac'] = df_ortho['HapA_GeneID'].map(chg_frac_dict)
df_ortho['HapB_CHG_frac'] = df_ortho['HapB_GeneID'].map(chg_frac_dict)

cat_order =['UM', 'gbM', 'teM', 'Unclassified', 'Missing']

# 4. Concordance Matrix
print("\n--- Concordance Matrix ---")
matrix = pd.crosstab(df_ortho['HapA_Class'], df_ortho['HapB_Class'])
matrix = matrix.reindex(index=cat_order, columns=cat_order, fill_value=0)
print(matrix)

# 5. Extract Valid Pairs
valid_states = ['UM', 'gbM', 'teM']
valid_mask = df_ortho['HapA_Class'].isin(valid_states) & df_ortho['HapB_Class'].isin(valid_states)
df_valid = df_ortho[valid_mask].copy()

# Categorize discordant types for coloring
def label_discordance(row):
    states = {row['HapA_Class'], row['HapB_Class']}
    if len(states) == 1:
        return "Concordant"
    elif states == {'UM', 'gbM'}:
        return "UM vs gbM"
    elif states == {'UM', 'teM'}:
        return "UM vs teM"
    elif states == {'gbM', 'teM'}:
        return "gbM vs teM"
    else:
        return "Other"

df_valid['Pair_Type'] = df_valid.apply(label_discordance, axis=1)

concordant = df_valid[df_valid['Pair_Type'] == "Concordant"]
discordant = df_valid[df_valid['Pair_Type'] != "Concordant"]

discordant.to_csv("discordant_epialleles.tsv", sep="\t", index=False)
print("\nSaved discordant pairs to 'discordant_epialleles.tsv'")

# 6. Plotting
print("Generating plots...")
sns.set_theme(style="whitegrid")

# --- Plot A: Heatmap ---
plt.figure(figsize=(7, 6))
sns.heatmap(matrix, annot=True, fmt="d", cmap="Blues", cbar=False, linewidths=.5)
plt.title("Ortholog Methylation State Concordance", fontweight="bold", pad=15)
plt.ylabel("HapA Classification")
plt.xlabel("HapB Classification")
plt.tight_layout()
plt.savefig("ortho_heatmap.png", dpi=300)
plt.close()


# --- Plot Setup for Panel Scatter ---
color_map = {
    "UM vs gbM": "#4A90D9",  # Blue
    "UM vs teM": "#E57373",  # Red
    "gbM vs teM": "#9C27B0"  # Purple
}

# Create a 1x2 panel figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
fig.suptitle("Allele-Specific Methylation: HapA vs HapB", fontweight="bold", fontsize=15, y=0.98)

# -- Subplot 1: CG Context --
ax1.scatter(concordant['HapA_CG_frac'], concordant['HapB_CG_frac'], 
            alpha=0.15, color='gray', s=10, label='Concordant')

for ptype, color in color_map.items():
    subset = discordant[discordant['Pair_Type'] == ptype]
    if not subset.empty:
        ax1.scatter(subset['HapA_CG_frac'], subset['HapB_CG_frac'], 
                    alpha=0.8, color=color, s=25, edgecolor='black', linewidth=0.5, label=ptype)

ax1.plot([0, 1], [0, 1], 'k--', lw=1.5, zorder=0)
ax1.set_title("CG Context", fontweight="bold", fontsize=13)
ax1.set_xlabel("HapA CG Fraction", fontsize=12)
ax1.set_ylabel("HapB CG Fraction", fontsize=12)
ax1.set_xlim(-0.02, 1.02)
ax1.set_ylim(-0.02, 1.02)
ax1.set_aspect('equal') # Forces the plot to be a perfect square

# -- Subplot 2: CHG Context --
ax2.scatter(concordant['HapA_CHG_frac'], concordant['HapB_CHG_frac'], 
            alpha=0.15, color='gray', s=10, label='Concordant')

for ptype, color in color_map.items():
    subset = discordant[discordant['Pair_Type'] == ptype]
    if not subset.empty:
        ax2.scatter(subset['HapA_CHG_frac'], subset['HapB_CHG_frac'], 
                    alpha=0.8, color=color, s=25, edgecolor='black', linewidth=0.5, label=ptype)

ax2.plot([0, 1], [0, 1], 'k--', lw=1.5, zorder=0)
ax2.set_title("CHG Context", fontweight="bold", fontsize=13)
ax2.set_xlabel("HapA CHG Fraction", fontsize=12)
ax2.set_ylabel("HapB CHG Fraction", fontsize=12)
ax2.set_xlim(-0.02, 1.02)
ax2.set_ylim(-0.02, 1.02)
ax2.set_aspect('equal') # Forces the plot to be a perfect square

# -- Combine Legend Outside --
# Extract handles and labels from the first axis (without altering them yet)
handles, labels = ax1.get_legend_handles_labels()

# Create the legend first
leg = fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.92, 0.5), 
                 fontsize=11, markerscale=1.5, frameon=True, title="Methylation State")

# NOW change the alpha specifically on the proxy artists inside the legend!
for lh in leg.legend_handles:
    lh.set_alpha(1)

# Adjust spacing so the plots don't get squished and leave room for the legend
plt.subplots_adjust(left=0.05, right=0.9, wspace=0.2)

# Save the panel
plt.savefig("ortho_scatter_panel.png", dpi=300, bbox_inches='tight')
plt.close()

print("Saved 'ortho_heatmap.png' and 'ortho_scatter_panel.png'")

#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

print("Loading data...")
# Load datasets
df_class = pd.read_csv("20-13_classification.tsv", sep="\t")
amps = pd.read_csv("significant_AMPs.tsv", sep="\t")

# Create dictionaries for mapping
class_dict = dict(zip(df_class['Gene_ID'], df_class['Classification']))
cg_dict = dict(zip(df_class['Gene_ID'], pd.to_numeric(df_class['CG_frac'], errors='coerce')))
chg_dict = dict(zip(df_class['Gene_ID'], pd.to_numeric(df_class['CHG_frac'], errors='coerce')))

# Map to the 1,240 AMPs
amps['HapA_Body'] = amps['HapA_GeneID'].map(class_dict).fillna("Unclassified")
amps['HapB_Body'] = amps['HapB_GeneID'].map(class_dict).fillna("Unclassified")

amps['HapA_CG'] = amps['HapA_GeneID'].map(cg_dict)
amps['HapB_CG'] = amps['HapB_GeneID'].map(cg_dict)

amps['HapA_CHG'] = amps['HapA_GeneID'].map(chg_dict)
amps['HapB_CHG'] = amps['HapB_GeneID'].map(chg_dict)

# Function to categorize the biological outcome of the gene body
def categorize_body_outcome(row):
    states = {row['HapA_Body'], row['HapB_Body']}
    if "Unclassified" in states:
        return "Involves Unclassified"
    elif states == {'gbM'}:
        return "Concordant gbM\n(Insulated)"
    elif states == {'UM'}:
        return "Concordant UM\n(Insulated)"
    elif states == {'teM'}:
        return "Concordant teM"
    elif states == {'UM', 'teM'}:
        return "UM vs teM\n(Silencing Spread)"
    elif states == {'UM', 'gbM'}:
        return "UM vs gbM"
    else:
        return "Other"

amps['Outcome'] = amps.apply(categorize_body_outcome, axis=1)

# Set up the Figure (1x3 Panel)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle("Gene-Body Methylation States of the 1,240 Asymmetrically Methylated Promoters (AMPs)", 
             fontweight="bold", fontsize=15, y=1.02)

# Define unified colors
color_map = {
    "Concordant gbM\n(Insulated)": "#4A90D9",  # Blue
    "Concordant UM\n(Insulated)": "gray",      # Gray
    "UM vs teM\n(Silencing Spread)": "#E57373",# Red
    "Concordant teM": "#9C27B0",               # Purple
    "UM vs gbM": "#81C784",                    # Green
    "Involves Unclassified": "#EEEEEE"         # Light gray for bars
}

# -------------------------------------------------------------
# PANEL A & B: Scatter Plots of the AMP Gene Bodies
# -------------------------------------------------------------
plot_data = amps[amps['Outcome'] != "Involves Unclassified"].copy()

for outcome, color in color_map.items():
    if outcome == "Involves Unclassified": continue
    
    subset = plot_data[plot_data['Outcome'] == outcome]
    if not subset.empty:
        # Use transparency to show density for the large concordant groups
        alpha_val = 0.4 if "Concordant UM" in outcome else 0.8
        
        # Plot on CG axis
        ax1.scatter(subset['HapA_CG'], subset['HapB_CG'], 
                    alpha=alpha_val, color=color, s=40, edgecolor='black', 
                    linewidth=0.5, label=outcome)
        
        # Plot on CHG axis
        ax2.scatter(subset['HapA_CHG'], subset['HapB_CHG'], 
                    alpha=alpha_val, color=color, s=40, edgecolor='black', 
                    linewidth=0.5)

# Format Panel A (CG)
ax1.plot([0, 1], [0, 1], 'k--', lw=1.5, zorder=0)
ax1.set_title("A. Gene-Body CG Fraction", fontweight="bold", fontsize=13)
ax1.set_xlabel("HapA Gene Body CG Fraction", fontsize=11)
ax1.set_ylabel("HapB Gene Body CG Fraction", fontsize=11)
ax1.set_xlim(-0.02, 1.02)
ax1.set_ylim(-0.02, 1.02)
ax1.set_aspect('equal')

# Format Panel B (CHG)
ax2.plot([0, 1],[0, 1], 'k--', lw=1.5, zorder=0)
ax2.set_title("B. Gene-Body CHG Fraction", fontweight="bold", fontsize=13)
ax2.set_xlabel("HapA Gene Body CHG Fraction", fontsize=11)
ax2.set_ylabel("HapB Gene Body CHG Fraction", fontsize=11)
ax2.set_xlim(-0.02, 1.02)
ax2.set_ylim(-0.02, 1.02)
ax2.set_aspect('equal')

# -------------------------------------------------------------
# PANEL C: Bar Chart Summarizing Outcomes
# -------------------------------------------------------------
counts = amps['Outcome'].value_counts()

bar_order =[
    "Concordant UM\n(Insulated)", 
    "Concordant gbM\n(Insulated)", 
    "UM vs teM\n(Silencing Spread)",
    "Concordant teM",
    "UM vs gbM",
    "Involves Unclassified"
]

bar_counts =[counts.get(cat, 0) for cat in bar_order]
bar_colors =[color_map[cat] for cat in bar_order]

bars = ax3.bar(bar_order, bar_counts, color=bar_colors, edgecolor='black', linewidth=0.8)

for bar in bars:
    yval = bar.get_height()
    ax3.text(bar.get_x() + bar.get_width()/2, yval + 10, int(yval), 
             ha='center', va='bottom', fontsize=10, fontweight='bold')

ax3.set_title("C. Categorical Outcomes of AMPs", fontweight="bold", fontsize=13)
ax3.set_ylabel("Number of AMP Pairs", fontsize=11)
ax3.set_xticklabels(bar_order, rotation=45, ha='right', fontsize=11)
ax3.set_ylim(0, max(bar_counts) * 1.15) 

# -------------------------------------------------------------
# Global Legend & Layout
# -------------------------------------------------------------
# Extract handles from ax1, fix the alpha so the legend dots are solid
handles, labels = ax1.get_legend_handles_labels()
leg = fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.93, 0.5), 
                 fontsize=11, markerscale=1.5, frameon=True, title="Gene Body State")

for lh in leg.legend_handles:
    lh.set_alpha(1)

plt.subplots_adjust(left=0.05, right=0.91, wspace=0.25, bottom=0.25)
plt.savefig("AMP_gene_bodies_figure_withCHG.png", dpi=300, bbox_inches='tight')
plt.close()

print("Saved 'AMP_gene_bodies_figure_withCHG.png'")

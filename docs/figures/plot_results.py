#!/usr/bin/env python3
"""Generate publication-quality figures for FUNGUS-SV yeast results."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

strains = ['S288C', 'BJ4', 'IMX2600', 'SX2', 'Makgeolli']
consensus = [277, 165, 314, 290, 250]
high_conf = [133, 55, 153, 152, 109]
perfect = [16, 13, 15, 19, 9]
deletions = [248, 140, 285, 261, 225]
inversions = [11, 12, 9, 11, 8]
duplications = [18, 13, 20, 18, 17]
colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0', '#F44336']

# Figure 1: SVs per strain
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(strains))
width = 0.35
bars1 = ax.bar(x - width/2, consensus, width, label='Consensus SVs', color='#2196F3', edgecolor='navy')
bars2 = ax.bar(x + width/2, high_conf, width, label='HIGH Confidence (T≥0.6)', color='#4CAF50', edgecolor='darkgreen')
for bar in bars1: ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+5, str(int(bar.get_height())), ha='center', fontsize=9)
for bar in bars2: ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+3, str(int(bar.get_height())), ha='center', fontsize=8)
ax.set_ylabel('Number of SVs')
ax.set_title('FUNGUS-SV: CICC-1445 Structural Variants per Reference Genome')
ax.set_xticks(x)
ax.set_xticklabels(strains)
ax.legend(loc='upper right')
ax.set_ylim(0, max(consensus)*1.2)
plt.tight_layout()
plt.savefig('docs/figures/figure1_svs_per_strain.png', dpi=150)
print("Saved figure1_svs_per_strain.png")

# Figure 2: Confidence distribution per strain
fig, ax = plt.subplots(figsize=(10, 6))
tiers = ['TRIPLE\n(T≥0.80)', 'DOUBLE\n(T≥0.60)', 'SINGLE\n(T≥0.40)', 'WEAK\n(T≥0.20)', 'CONTRADICTED\n(<0.20)']
tier_data = {
    'S288C': [34, 99, 16, 113, 15],
    'BJ4': [21, 34, 29, 69, 12],
    'IMX2600': [44, 109, 30, 115, 16],
    'SX2': [38, 114, 19, 102, 17],
    'Makgeolli': [23, 86, 35, 95, 11],
}
x = np.arange(len(tiers))
width = 0.15
for i, (strain, counts) in enumerate(tier_data.items()):
    ax.bar(x + i*width, counts, width, label=strain, color=colors[i], edgecolor='black', alpha=0.85)
ax.set_ylabel('Number of SVs')
ax.set_title('Confidence Tier Distribution by Reference Strain')
ax.set_xticks(x + width*2)
ax.set_xticklabels(tiers)
ax.legend(loc='upper right')
plt.tight_layout()
plt.savefig('docs/figures/figure2_confidence_tiers.png', dpi=150)
print("Saved figure2_confidence_tiers.png")

# Figure 3: SV type breakdown
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(strains))
width = 0.25
ax.bar(x - width, deletions, width, label='DEL', color='#E53935', edgecolor='darkred')
ax.bar(x, inversions, width, label='INV', color='#FFB300', edgecolor='darkorange')
ax.bar(x + width, duplications, width, label='DUP', color='#1E88E5', edgecolor='darkblue')
ax.set_ylabel('Number of SVs')
ax.set_title('SV Type Distribution per Reference Strain')
ax.set_xticks(x)
ax.set_xticklabels(strains)
ax.legend()
plt.tight_layout()
plt.savefig('docs/figures/figure3_sv_types.png', dpi=150)
print("Saved figure3_sv_types.png")

# Figure 4: Phylogenetic distance
fig, ax = plt.subplots(figsize=(8, 5))
strains_sorted = ['BJ4', 'Makgeolli', 'S288C', 'SX2', 'IMX2600']
sv_counts = [165, 250, 277, 290, 314]
bar_colors = ['#4CAF50', '#8BC34A', '#FFC107', '#FF9800', '#F44336']
ax.barh(strains_sorted, sv_counts, color=bar_colors, edgecolor='black')
ax.set_xlabel('Number of Consensus SVs (fewer = closer relative)')
ax.set_title('CICC-1445 Phylogenetic Distance (based on SV count)')
for i, v in enumerate(sv_counts):
    ax.text(v + 3, i, str(v), va='center')
plt.tight_layout()
plt.savefig('docs/figures/figure4_phylogenetic_distance.png', dpi=150)
print("Saved figure4_phylogenetic_distance.png")

print("\nAll figures generated in docs/figures/")

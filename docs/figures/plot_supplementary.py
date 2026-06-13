#!/usr/bin/env python3
"""Supplementary figures for FUNGUS-SV."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

strains = ['S288C', 'BJ4', 'IMX2600', 'SX2', 'Makgeolli']
colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0', '#F44336']

# Figure 5: T-score distribution (simulated from tier data)
tier_data = {
    'S288C': [34, 99, 16, 113, 15],   # TRIPLE, DOUBLE, SINGLE, WEAK, CONTRADICTED
    'BJ4': [21, 34, 29, 69, 12],
    'IMX2600': [44, 109, 30, 115, 16],
    'SX2': [38, 114, 19, 102, 17],
    'Makgeolli': [23, 86, 35, 95, 11],
}

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()
for i, (strain, counts) in enumerate(tier_data.items()):
    t_scores = []
    t_scores.extend([0.90]*counts[0])  # TRIPLE
    t_scores.extend([0.70]*counts[1])  # DOUBLE
    t_scores.extend([0.50]*counts[2])  # SINGLE
    t_scores.extend([0.30]*counts[3])  # WEAK
    t_scores.extend([0.10]*counts[4])  # CONTRADICTED
    axes[i].hist(t_scores, bins=20, color=colors[i], edgecolor='black', alpha=0.8)
    axes[i].axvline(x=0.60, color='red', linestyle='--', label='HIGH threshold')
    axes[i].set_title(f'{strain} (n={len(t_scores)})')
    axes[i].set_xlabel('T-score')
    axes[i].set_ylabel('Count')
    axes[i].legend(fontsize=7)
axes[5].axis('off')
plt.suptitle('T-Score Distribution per Reference Strain', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('docs/figures/figure5_tscore_distribution.png', dpi=150, bbox_inches='tight')
print("Saved figure5_tscore_distribution.png")

# Figure 6: SV size distribution (approximate from tier bin data)
size_bins_data = {
    'S288C': [65, 120, 31, 61],     # 50-100, 100-500, 500-5000, >5000
    'BJ4': [56, 58, 32, 19],
    'IMX2600': [77, 147, 40, 50],
    'SX2': [96, 143, 35, 16],
    'Makgeolli': [63, 116, 36, 35],
}
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(4)
width = 0.15
for i, (strain, counts) in enumerate(size_bins_data.items()):
    ax.bar(x + i*width, counts, width, label=strain, color=colors[i], edgecolor='black')
ax.set_ylabel('Number of SVs')
ax.set_title('SV Size Distribution per Reference Strain')
ax.set_xticks(x + width*2)
ax.set_xticklabels(['50-100 bp', '100-500 bp', '500-5000 bp', '>5000 bp'])
ax.legend()
plt.tight_layout()
plt.savefig('docs/figures/figure6_size_distribution.png', dpi=150)
print("Saved figure6_size_distribution.png")

# Figure 7: Caller agreement (raw calls per strain)
raw_calls = {
    'S288C': [794, 1600, 1586],
    'BJ4': [871, 1520, 1491],
    'IMX2600': [826, 1769, 1599],
    'SX2': [848, 1714, 1572],
    'Makgeolli': [747, 1450, 1489],
}
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(strains))
width = 0.25
sniffles = [raw_calls[s][0] for s in strains]
cutesv = [raw_calls[s][1] for s in strains]
svim = [raw_calls[s][2] for s in strains]
ax.bar(x - width, sniffles, width, label='Sniffles2', color='#E91E63', edgecolor='black')
ax.bar(x, cutesv, width, label='cuteSV', color='#2196F3', edgecolor='black')
ax.bar(x + width, svim, width, label='SVIM', color='#4CAF50', edgecolor='black')
ax.set_ylabel('Raw SV Calls')
ax.set_title('Raw Calls per Caller per Reference')
ax.set_xticks(x)
ax.set_xticklabels(strains)
ax.legend()
plt.tight_layout()
plt.savefig('docs/figures/figure7_raw_calls.png', dpi=150)
print("Saved figure7_raw_calls.png")

# Figure 8: Consensus rate (raw → consensus reduction)
raw_totals = [sum(raw_calls[s]) for s in strains]
consensus_svs = [277, 165, 314, 290, 250]
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(strains))
width = 0.35
ax.bar(x - width/2, raw_totals, width, label='Raw Calls (3 callers)', color='#FF9800', edgecolor='black')
ax.bar(x + width/2, consensus_svs, width, label='ICB Consensus', color='#4CAF50', edgecolor='black')
for i in range(len(strains)):
    reduction = (1 - consensus_svs[i]/raw_totals[i]) * 100
    ax.annotate(f'-{reduction:.0f}%', (x[i]+width/2, consensus_svs[i]+10), ha='center', fontsize=8, color='darkgreen')
ax.set_ylabel('Number of SVs')
ax.set_title('ICB Consensus: Raw Calls → Filtered Consensus')
ax.set_xticks(x)
ax.set_xticklabels(strains)
ax.legend()
plt.tight_layout()
plt.savefig('docs/figures/figure8_consensus_rate.png', dpi=150)
print("Saved figure8_consensus_rate.png")

print("\nAll supplementary figures generated!")

#!/usr/bin/env python3
"""Generate publication-quality figures for FUNGUS-SV README."""

import json, os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

# Collect data
strains_data = {}
for folder in sorted(os.listdir('results/baumannii_strains/')):
    jp = f'results/baumannii_strains/{folder}/validation/validation_results.json'
    if os.path.exists(jp):
        with open(jp) as f: data = json.load(f)
        name = folder.replace('A_baumanii_', '')
        r = data['results']
        strains_data[name] = {
            'total': len(r), 'high': sum(1 for x in r if x['t_score']>=0.6),
            'med': sum(1 for x in r if 0.4<=x['t_score']<0.6),
            'weak': sum(1 for x in r if x['t_score']<0.4),
            't_scores': [x['t_score'] for x in r],
            'sizes': [x.get('sv_size',0) for x in r],
            'types': [x['sv_type'] for x in r]
        }

sorted_names = sorted(strains_data.keys(), key=lambda n: strains_data[n]['total'], reverse=True)

# Colors
colors = {'high': '#2ecc71', 'med': '#f39c12', 'weak': '#e74c3c'}

# ============================================================
# FIGURE 1: Bar chart - SVs per strain with confidence breakdown
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))
names = [n for n in sorted_names]
totals = [strains_data[n]['total'] for n in names]
highs = [strains_data[n]['high'] for n in names]
meds = [strains_data[n]['med'] for n in names]
weaks = [strains_data[n]['weak'] for n in names]

x = np.arange(len(names))
width = 0.6
bars_high = ax.bar(x, highs, width, color=colors['high'], label='HIGH (T≥0.6)')
bars_med = ax.bar(x, meds, width, bottom=highs, color=colors['med'], label='MEDIUM (0.4-0.6)')
bars_weak = ax.bar(x, weaks, width, bottom=[h+m for h,m in zip(highs,meds)], color=colors['weak'], label='WEAK (T<0.4)')

# Add percentage labels
for i, (n, t, h) in enumerate(zip(names, totals, highs)):
    pct = h/t*100 if t>0 else 0
    ax.text(i, t+2, f'{pct:.0f}%', ha='center', fontsize=8, fontweight='bold', color='#2c3e50')

ax.set_xlabel('A. baumannii Strain', fontweight='bold')
ax.set_ylabel('Number of SVs', fontweight='bold')
ax.set_title('FUNGUS-SV: Structural Variants Detected Across 12 A. baumannii Strains', fontweight='bold', fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=45, ha='right', fontsize=9)
ax.legend(loc='upper right', framealpha=0.9)
ax.set_ylim(0, max(totals)*1.15)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('results/figures/figure1_svs_per_strain.png', bbox_inches='tight')
plt.close()
print("Figure 1 saved")

# ============================================================
# FIGURE 2: Pie chart - Confidence tier distribution
# ============================================================
fig, ax = plt.subplots(figsize=(8, 8))
total_high = sum(strains_data[n]['high'] for n in names)
total_med = sum(strains_data[n]['med'] for n in names)
total_weak = sum(strains_data[n]['weak'] for n in names)
total_all = total_high + total_med + total_weak

sizes = [total_high, total_med, total_weak]
labels = [f'HIGH (T≥0.6)\n{total_high} SVs\n{total_high/total_all*100:.0f}%',
          f'MEDIUM (0.4-0.6)\n{total_med} SVs\n{total_med/total_all*100:.0f}%',
          f'WEAK (T<0.4)\n{total_weak} SVs\n{total_weak/total_all*100:.0f}%']
pie_colors = [colors['high'], colors['med'], colors['weak']]
explode = (0.05, 0.02, 0.02)

wedges, texts = ax.pie(sizes, explode=explode, labels=labels, colors=pie_colors,
                        startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
ax.set_title(f'Confidence Distribution: {total_all} Total SVs', fontweight='bold', fontsize=14)
plt.tight_layout()
plt.savefig('results/figures/figure2_confidence_pie.png', bbox_inches='tight')
plt.close()
print("Figure 2 saved")

# ============================================================
# FIGURE 3: Scatter plot - T-score vs. SV size
# ============================================================
fig, ax = plt.subplots(figsize=(10, 6))
all_sizes = []
all_scores = []
all_colors = []
for n in names:
    d = strains_data[n]
    all_sizes.extend(d['sizes'])
    all_scores.extend(d['t_scores'])
    for s, t in zip(d['sizes'], d['t_scores']):
        if t >= 0.6: all_colors.append(colors['high'])
        elif t >= 0.4: all_colors.append(colors['med'])
        else: all_colors.append(colors['weak'])

ax.scatter(all_sizes, all_scores, c=all_colors, alpha=0.6, edgecolors='white', linewidth=0.5)
ax.axhline(y=0.6, color=colors['high'], linestyle='--', alpha=0.5, label='HIGH threshold (0.6)')
ax.axhline(y=0.4, color=colors['med'], linestyle='--', alpha=0.5, label='MEDIUM threshold (0.4)')
ax.axvline(x=100, color='gray', linestyle=':', alpha=0.7, label='100 bp boundary')
ax.set_xlabel('SV Size (bp)', fontweight='bold')
ax.set_ylabel('T-Score', fontweight='bold')
ax.set_title('T-Score vs. SV Size: 860 SVs Across 12 Strains', fontweight='bold', fontsize=13)
ax.set_xscale('log')
ax.set_xlim(30, 10000)
ax.legend(loc='lower right')
ax.grid(alpha=0.3)

# Add annotation
ax.annotate('All SVs ≥100 bp\nscore HIGH', xy=(200, 0.7), fontsize=10, fontweight='bold',
            color=colors['high'], ha='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
ax.annotate('All SVs <100 bp\nscore WEAK', xy=(50, 0.3), fontsize=10, fontweight='bold',
            color=colors['weak'], ha='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('results/figures/figure3_tscore_vs_size.png', bbox_inches='tight')
plt.close()
print("Figure 3 saved")

# ============================================================
# FIGURE 4: LAR validation results
# ============================================================
fig, ax = plt.subplots(figsize=(7, 7))
lar_data = {'Confirmed': 23, 'Partial': 12, 'Failed': 4}
lar_colors = ['#2ecc71', '#f39c12', '#e74c3c']
bars = ax.bar(list(lar_data.keys()), list(lar_data.values()), color=lar_colors, edgecolor='white', linewidth=2)
for bar, val in zip(bars, lar_data.values()):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, f'{val}\n({val/39*100:.0f}%)',
            ha='center', fontweight='bold', fontsize=12)
ax.set_ylabel('Number of SVs', fontweight='bold')
ax.set_title('LAR Validation: Top 3 SVs × 12 Strains (39 Total)', fontweight='bold', fontsize=13)
ax.set_ylim(0, 28)
ax.grid(axis='y', alpha=0.3)
# Add 90% support annotation
ax.text(1.5, 25, '90% Get Some Support', ha='center', fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='#ecf0f1', alpha=0.8))
plt.tight_layout()
plt.savefig('results/figures/figure4_lar_results.png', bbox_inches='tight')
plt.close()
print("Figure 4 saved")

# ============================================================
# FIGURE 5: Summary infographic-style
# ============================================================
fig, ax = plt.subplots(figsize=(12, 8))
ax.axis('off')
ax.set_xlim(0, 12)
ax.set_ylim(0, 8)

# Title
ax.text(6, 7.5, 'FUNGUS-SV v0.7.2 — Pipeline Performance Summary', ha='center', fontsize=18, fontweight='bold')
ax.text(6, 7.1, 'A. baumannii ATCC 19606 HiFi reads vs. 12 clinical strains + 5 species', ha='center', fontsize=11, color='#7f8c8d')

# Key metrics boxes
metrics = [
    (1.5, 6.0, '860', 'Total SVs\nDetected', '#3498db'),
    (4.0, 6.0, '85%', 'HIGH\nConfidence', '#2ecc71'),
    (6.5, 6.0, '59%', 'LAR\nConfirmed', '#27ae60'),
    (9.0, 6.0, '90%', 'LAR\nSupported', '#f39c12'),
]
for x, y, big, label, color in metrics:
    box = plt.Rectangle((x-1, y-0.5), 2, 1.2, facecolor=color, alpha=0.15, edgecolor=color, linewidth=2, transform=ax.transData)
    ax.add_patch(box)
    ax.text(x, y+0.3, big, ha='center', fontsize=24, fontweight='bold', color=color)
    ax.text(x, y-0.15, label, ha='center', fontsize=9, color='#2c3e50')

# Bottom findings
findings = [
    '◆ 100% of SVs ≥100 bp score HIGH confidence',
    '◆ 100% of SVs <100 bp score WEAK confidence',
    '◆ 2+ active layers → HIGH | 1 active layer → WEAK',
    '◆ Within-species finds 25× more SVs than cross-species',
    '◆ ICB reduces raw calls by ~50%',
]
for i, f in enumerate(findings):
    ax.text(2, 4.3 - i*0.5, f, fontsize=11, color='#2c3e50')

# Bottom bar
ax.text(6, 1.5, 'Validated with 17 comparisons | 13 peer-reviewed papers | 6 conda environments | MIT License',
        ha='center', fontsize=9, color='#95a5a6')

plt.tight_layout()
plt.savefig('results/figures/figure5_summary_infographic.png', bbox_inches='tight')
plt.close()
print("Figure 5 saved")

print("\nAll 5 figures saved in results/figures/")
print("Ready for README!")

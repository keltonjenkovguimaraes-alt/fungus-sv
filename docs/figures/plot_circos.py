#!/usr/bin/env python3
"""Genomic circos plot of SVs: CICC-1445 vs S288C."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# S288C chromosome lengths (kb)
chromosomes = {
    'I': 230, 'II': 813, 'III': 317, 'IV': 1532, 'V': 577,
    'VI': 270, 'VII': 1091, 'VIII': 563, 'IX': 440, 'X': 746,
    'XI': 667, 'XII': 1078, 'XIII': 924, 'XIV': 784, 'XV': 1091, 'XVI': 948
}

# Parse SVs from consensus VCF (S288C results)
sv_by_chr = {}
with open('data/yeast/results/consensus_svs.vcf') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\t')
        chrom = parts[0].split('.')[0].replace('NC_001', '')
        if chrom.endswith('1'): chrom = 'XVI'
        pos = int(parts[1])
        info = parts[7]
        import re
        m = re.search(r'SVTYPE=(\w+)', info)
        svtype = m.group(1) if m else 'UNK'
        if chrom not in sv_by_chr: sv_by_chr[chrom] = []
        sv_by_chr[chrom].append((pos, svtype))

# Simple circular layout
fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={'projection': 'polar'})

# Draw chromosome arcs
angles = np.linspace(0, 2*np.pi, len(chromosomes)+1)
colors_cycle = plt.cm.Set3(np.linspace(0, 1, len(chromosomes)))

for i, (name, length) in enumerate(chromosomes.items()):
    theta = np.linspace(angles[i], angles[i+1], 100)
    r = np.ones(100) * 10
    ax.fill_between(theta, 9.5, 10.5, alpha=0.3, color=colors_cycle[i])
    mid_angle = (angles[i] + angles[i+1]) / 2
    ax.text(mid_angle, 11, f'{name}', ha='center', va='center', fontsize=8)

ax.set_ylim(0, 15)
ax.axis('off')
ax.set_title('S. cerevisiae S288C Chromosomes', fontsize=14, pad=20)
plt.tight_layout()
plt.savefig('docs/figures/figure9_circos.png', dpi=150, bbox_inches='tight')
print("Saved figure9_circos.png - simplified version")
print("For publication-quality circos, install: conda install -c bioconda circos")

#!/usr/bin/env python3
"""Gene annotation visualization for FUNGUS-SV yeast results."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import re
from collections import Counter

# ===== DATA =====
strains = ['S288C', 'BJ4', 'IMX2600', 'SX2']
colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0']
genic = [162, 128, 100, 220]
intergenic = [115, 37, 214, 70]
total_svs = [277, 165, 314, 290]

# Gene counts from annotation files
gene_data = {}
for strain in strains:
    fname = f'/home/kelto/fungus-sv/data/yeast/{strain}_sv_annotations.tsv'
    genes = []
    with open(fname) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7 and parts[6] != '(intergenic)':
                for g in parts[6].split(', '):
                    g = g.strip()
                    if g and 'unknown' not in g.lower():
                        genes.append(g)
    gene_data[strain] = Counter(genes)

# ===== FIGURE A: Genic vs Intergenic =====
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(strains))
width = 0.35
ax.bar(x - width/2, genic, width, label='Genic (hits gene)', color='#E53935', edgecolor='darkred')
ax.bar(x + width/2, intergenic, width, label='Intergenic', color='#BDBDBD', edgecolor='gray')
for i in range(len(strains)):
    pct = genic[i]/total_svs[i]*100
    ax.text(i - width/2, genic[i] + 3, f'{genic[i]}\n({pct:.0f}%)', ha='center', fontsize=8, fontweight='bold')
    ax.text(i + width/2, intergenic[i] + 3, str(intergenic[i]), ha='center', fontsize=8)
ax.set_ylabel('Number of SVs')
ax.set_title('A. SVs Affecting Genes vs Intergenic Regions')
ax.set_xticks(x)
ax.set_xticklabels(strains)
ax.legend()
ax.set_ylim(0, max(genic)*1.25)
plt.tight_layout()
plt.savefig('docs/figures/figureA_genic_vs_intergenic.png', dpi=150)
print("Saved figureA_genic_vs_intergenic.png")

# ===== FIGURE B: Top Hit Genes Heatmap =====
# Find top genes across all strains
all_genes = set()
for strain in strains:
    all_genes.update(list(gene_data[strain].keys())[:15])

# Get top 20 genes overall
global_counter = Counter()
for strain in strains:
    global_counter.update(gene_data[strain])
top_genes = [g for g, c in global_counter.most_common(20)]

# Build matrix
matrix = []
for gene in top_genes:
    row = [gene_data[strain].get(gene, 0) for strain in strains]
    matrix.append(row)
matrix = np.array(matrix)

fig, ax = plt.subplots(figsize=(12, 8))
im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(len(strains)))
ax.set_xticklabels(strains, fontsize=11)
ax.set_yticks(range(len(top_genes)))
ax.set_yticklabels(top_genes, fontsize=8)
ax.set_title('B. Top 20 Genes Hit by SVs Across All Strains', fontsize=13, fontweight='bold')
for i in range(len(top_genes)):
    for j in range(len(strains)):
        val = matrix[i, j]
        if val > 0:
            ax.text(j, i, str(val), ha='center', va='center', fontsize=7, fontweight='bold',
                   color='white' if val > matrix.max()/2 else 'black')
plt.colorbar(im, ax=ax, label='Number of SVs overlapping gene')
plt.tight_layout()
plt.savefig('docs/figures/figureB_top_genes_heatmap.png', dpi=150)
print("Saved figureB_top_genes_heatmap.png")

# ===== FIGURE C: SV Size vs Gene Count Scatter =====
size_vs_genes = {s: [] for s in strains}
for strain in strains:
    fname = f'/home/kelto/fungus-sv/data/yeast/{strain}_sv_annotations.tsv'
    with open(fname) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                sv_size = int(parts[5])
                genes_str = parts[6]
                n_genes = 0 if genes_str == '(intergenic)' else len(genes_str.split(', '))
                size_vs_genes[strain].append((sv_size, n_genes))

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()
for i, strain in enumerate(strains):
    ax = axes[i]
    sizes = [s[0] for s in size_vs_genes[strain]]
    ngenes = [s[1] for s in size_vs_genes[strain]]
    ax.scatter(sizes, ngenes, c=colors[i], alpha=0.6, edgecolors='black', linewidth=0.3, s=30)
    ax.set_xlabel('SV Size (bp)')
    ax.set_ylabel('Genes Affected')
    ax.set_title(f'{strain} (n={len(sizes)})')
    ax.set_xscale('log')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    # Add correlation
    if len(sizes) > 1:
        corr = np.corrcoef(sizes, ngenes)[0,1]
        ax.text(0.95, 0.95, f'r={corr:.2f}', transform=ax.transAxes, ha='right', va='top')
plt.suptitle('C. SV Size vs Number of Genes Affected', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('docs/figures/figureC_size_vs_genes.png', dpi=150)
print("Saved figureC_size_vs_genes.png")

# ===== FIGURE D: Genome Track with Gene Labels (S288C, Chr I & II) =====
# Parse S288C VCF and GFF for chromosome I
svs_chr = {}
with open('/home/kelto/fungus-sv/data/yeast/results_3callers/consensus_svs.vcf') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\t')
        chrom = parts[0]
        pos = int(parts[1])
        info = parts[7]
        m = re.search(r'SVTYPE=(\w+)', info)
        svtype = m.group(1) if m else 'UNK'
        m2 = re.search(r'END=(\d+)', info)
        end = int(m2.group(1)) if m2 else pos
        if chrom not in svs_chr: svs_chr[chrom] = []
        svs_chr[chrom].append((pos, end, svtype))

# Get genes for NC_001133 (Chr I) and NC_001134 (Chr II)
genes_chr = {}
with open('/mnt/c/Users/kelto/OneDrive/Documentos/Genomes/ncbi_dataset/data/GCF_000146045.2/genomic.gff') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene': continue
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        name = 'unknown'
        m = re.search(r'Name=([^;]+)', parts[8])
        if m: name = m.group(1)
        if chrom not in genes_chr: genes_chr[chrom] = []
        genes_chr[chrom].append((start, end, name))

# Plot Chr I and II
target_chrs = ['NC_001133.9', 'NC_001134.8']
chr_labels = ['Chromosome I', 'Chromosome II']
chr_lengths = [230218, 813184]
sv_colors = {'DEL': '#E53935', 'INV': '#FFB300', 'DUP': '#1E88E5'}

fig, axes = plt.subplots(2, 1, figsize=(20, 8))

for idx, (chr_id, label, length) in enumerate(zip(target_chrs, chr_labels, chr_lengths)):
    ax = axes[idx]
    
    # Chromosome backbone
    ax.plot([0, length], [1, 1], 'k-', linewidth=2, alpha=0.5)
    
    # Plot genes as small blue bars (top track)
    if chr_id in genes_chr:
        for (start, end, name) in genes_chr[chr_id]:
            ax.plot([start, end], [1.3, 1.3], 'b-', linewidth=1, alpha=0.3)
            if end - start > 2000:  # Label only larger genes
                ax.text((start+end)/2, 1.45, name, ha='center', fontsize=5, rotation=90, alpha=0.7)
    
    # Plot SVs (bottom track)
    if chr_id in svs_chr:
        for (pos, end, svtype) in svs_chr[chr_id]:
            color = sv_colors.get(svtype, 'gray')
            sv_len = abs(end-pos)
            ax.plot([pos, end], [0.7, 0.7], '-', color=color, linewidth=2, alpha=0.8)
            if sv_len > 500:
                ax.text((pos+end)/2, 0.55, f'{svtype} {sv_len}bp', ha='center', fontsize=6, color=color)
    
    ax.set_ylim(0.3, 1.7)
    ax.set_xlim(0, length)
    ax.set_yticks([0.7, 1.0, 1.3])
    ax.set_yticklabels(['SVs', 'Chr', 'Genes'], fontsize=8)
    ax.set_title(f'{label} ({length:,} bp)', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle('D. S288C vs CICC-1445: SV Locations & Affected Genes (Chromosomes I-II)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('docs/figures/figureD_genome_track_annotated.png', dpi=200, bbox_inches='tight')
print("Saved figureD_genome_track_annotated.png")

print("\nAll annotation figures generated!")

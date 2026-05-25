#!/usr/bin/env python3
"""Linear genome tracks: SV positions per chromosome for all 5 strains."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import re
import os

os.chdir('/home/kelto/fungus-sv')

# S288C chromosome lengths (bp)
chrom_lengths = {
    'I': 230218, 'II': 813184, 'III': 316620, 'IV': 1531933, 'V': 576874,
    'VI': 270161, 'VII': 1090940, 'VIII': 562643, 'IX': 439888, 'X': 745751,
    'XI': 666816, 'XII': 1078177, 'XIII': 924431, 'XIV': 784333, 'XV': 1091291, 'XVI': 948066
}

def get_chrom_roman(acc):
    """Map any chromosome naming scheme to roman numeral."""
    # NCBI RefSeq (NC_001133.9 -> I)
    ncbi_map = {
        'NC_001133': 'I', 'NC_001134': 'II', 'NC_001135': 'III', 'NC_001136': 'IV',
        'NC_001137': 'V', 'NC_001138': 'VI', 'NC_001139': 'VII', 'NC_001140': 'VIII',
        'NC_001141': 'IX', 'NC_001142': 'X', 'NC_001143': 'XI', 'NC_001144': 'XII',
        'NC_001145': 'XIII', 'NC_001146': 'XIV', 'NC_001147': 'XV', 'NC_001148': 'XVI',
        'NC_001224': 'MT'
    }
    for key in ncbi_map:
        if key in acc: return ncbi_map[key]
    # ENA/EMBL (LR813517.2 -> chromosome 1 -> I)
    chr_num = re.search(r'chromosome:?\s*(\d+)', acc)
    if chr_num:
        num = int(chr_num.group(1))
        roman = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
        if 1 <= num <= 16: return roman[num-1]
    # GenBank CP (CP127195.1 -> chromosome I)
    cp_num = re.search(r'chromosome\s+([XVI]+)', acc, re.IGNORECASE)
    if cp_num: return cp_num.group(1)
    # Try to find a roman numeral anywhere
    roman_match = re.search(r'(X{0,3})(IX|IV|V?I{0,3})', acc)
    if roman_match: return roman_match.group(0)
    return None

def parse_vcf(vcf_path):
    """Parse consensus VCF, return SVs per chromosome."""
    svs = {}
    if not os.path.exists(vcf_path):
        print(f"    FILE NOT FOUND: {vcf_path}")
        return svs
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 8: continue
            acc = parts[0]
            chrom = get_chrom_roman(acc)
            if chrom is None or chrom == 'MT': continue
            pos = int(parts[1])
            info = parts[7]
            m = re.search(r'SVTYPE=(\w+)', info)
            svtype = m.group(1) if m else 'UNK'
            m2 = re.search(r'END=(\d+)', info)
            end = int(m2.group(1)) if m2 else pos
            if chrom not in svs: svs[chrom] = []
            svs[chrom].append((pos, end, svtype))
    return svs

# Parse all 5 strains
strain_dirs = {
    'S288C': 'data/yeast/results_3callers/consensus_svs.vcf',
    'BJ4': 'data/yeast/results_BJ4/consensus_svs.vcf',
    'IMX2600': 'data/yeast/results_IMX2600/consensus_svs.vcf',
    'SX2': 'data/yeast/results_SX2/consensus_svs.vcf',
    'Makgeolli': 'data/yeast/results_Makgeolli/consensus_svs.vcf',
}

strains_data = {}
for strain, path in strain_dirs.items():
    print(f"Parsing {strain}: {path}")
    strains_data[strain] = parse_vcf(path)
    total = sum(len(v) for v in strains_data[strain].values())
    print(f"  -> {total} SVs in {len(strains_data[strain])} chromosomes")

# Colors
sv_colors = {'DEL': '#E53935', 'INS': '#43A047', 'INV': '#FFB300', 'DUP': '#1E88E5'}

# Create figure
strains_ordered = ['S288C', 'BJ4', 'IMX2600', 'SX2', 'Makgeolli']
chrs_ordered = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']
max_pos = max(chrom_lengths.values())

fig, axes = plt.subplots(5, 1, figsize=(22, 16))
fig.subplots_adjust(hspace=0.45)

for strain_idx, strain in enumerate(strains_ordered):
    ax = axes[strain_idx]
    sv_data = strains_data.get(strain, {})
    total_svs = sum(len(v) for v in sv_data.values())
    
    for chr_idx, chrom in enumerate(chrs_ordered):
        y_pos = len(chrs_ordered) - chr_idx
        length = chrom_lengths.get(chrom, 500000)
        
        # Chromosome line
        ax.plot([0, length], [y_pos, y_pos], 'k-', linewidth=0.4, alpha=0.25)
        
        # SV markers
        if chrom in sv_data:
            for (pos, end, svtype) in sv_data[chrom]:
                color = sv_colors.get(svtype, 'gray')
                sv_len = abs(end - pos)
                size = max(15, min(sv_len/1000, 100))
                ax.scatter(pos, y_pos, s=size, c=color, alpha=0.75, edgecolors='black', linewidth=0.2, zorder=3)
    
    ax.set_ylim(0.5, 17)
    ax.set_yticks(range(1, 17))
    ax.set_yticklabels(chrs_ordered[::-1], fontsize=8)
    ax.set_ylabel(f'{strain}\n({total_svs} SVs)', fontsize=10, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, max_pos * 1.02)
    ax.tick_params(axis='x', labelsize=8)

axes[-1].set_xlabel('Genomic Position (bp) — Chromosome lengths from S288C reference', fontsize=11)

# Legend
legend_patches = [mpatches.Patch(color=sv_colors[t], label=f'{t}') for t in ['DEL', 'INV', 'DUP']]
fig.legend(handles=legend_patches, loc='upper right', fontsize=10, ncol=3, frameon=True)

plt.suptitle('FUNGUS-SV: CICC-1445 SVs Across 5 S. cerevisiae Reference Genomes\nDot size ∝ SV length | Red=DEL, Yellow=INV, Blue=DUP',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('docs/figures/figure_genome_tracks.png', dpi=200, bbox_inches='tight')
print("\nSaved figure_genome_tracks.png")

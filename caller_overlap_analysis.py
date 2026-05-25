#!/usr/bin/env python3
"""Deep comparison: Sniffles2 vs cuteSV vs ICB overlap."""

import os, re
from collections import defaultdict

def parse_vcf_simple(vcf_file):
    """Parse VCF to set of (chrom, pos, end, type) tuples."""
    svs = set()
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            info = parts[7]
            m_type = re.search(r'SVTYPE=(\w+)', info)
            svtype = m_type.group(1) if m_type else 'UNK'
            m_end = re.search(r'END=(\d+)', info)
            end = int(m_end.group(1)) if m_end else pos
            svs.add((chrom, pos, end, svtype))
    return svs

def overlap(sv1, sv2, max_dist=2000):
    """Check if two SVs overlap (same chrom, type, within max_dist)."""
    return (sv1[0] == sv2[0] and sv1[3] == sv2[3] and 
            abs(sv1[1] - sv2[1]) <= max_dist)

def count_overlap(set_a, set_b, max_dist=2000):
    """Count how many SVs in set_a have a match in set_b."""
    count = 0
    for sv_a in set_a:
        for sv_b in set_b:
            if overlap(sv_a, sv_b, max_dist):
                count += 1
                break
    return count

# Collect all data
all_snf = set()
all_csv = set()
all_icb = set()
strain_stats = []

for folder in sorted(os.listdir('results/baumannii_strains/')):
    base = f'results/baumannii_strains/{folder}'
    name = folder.replace('A_baumanii_', '')
    
    snf_file = f'{base}/variants/sniffles2_svs.vcf'
    csv_file = f'{base}/variants/cutesv_svs.vcf'
    icb_file = f'{base}/variants/consensus_svs.vcf'
    
    if not all(os.path.exists(f) for f in [snf_file, csv_file, icb_file]):
        continue
    
    snf_set = parse_vcf_simple(snf_file)
    csv_set = parse_vcf_simple(csv_file)
    icb_set = parse_vcf_simple(icb_file)
    
    all_snf.update(snf_set)
    all_csv.update(csv_set)
    all_icb.update(icb_set)
    
    # Per-strain overlap analysis
    snf_in_csv = count_overlap(snf_set, csv_set)
    csv_in_snf = count_overlap(csv_set, snf_set)
    icb_in_snf = count_overlap(icb_set, snf_set)
    icb_in_csv = count_overlap(icb_set, csv_set)
    icb_in_both = sum(1 for sv in icb_set 
                      if any(overlap(sv, s, 500) for s in snf_set) 
                      and any(overlap(sv, s, 500) for s in csv_set))
    icb_in_one_only = len(icb_set) - icb_in_both
    
    strain_stats.append({
        'name': name,
        'snf_total': len(snf_set), 'csv_total': len(csv_set), 'icb_total': len(icb_set),
        'snf_in_csv': snf_in_csv, 'csv_in_snf': csv_in_snf,
        'icb_in_snf': icb_in_snf, 'icb_in_csv': icb_in_csv,
        'icb_in_both': icb_in_both, 'icb_in_one_only': icb_in_one_only
    })

# Print results
print("\n" + "=" * 100)
print("  CALLER OVERLAP ANALYSIS: Sniffles2 vs. cuteSV vs. ICB Consensus")
print("=" * 100)

# Q1: Sniffles2 ↔ cuteSV overlap
snf_in_csv_total = sum(s['snf_in_csv'] for s in strain_stats)
csv_in_snf_total = sum(s['csv_in_snf'] for s in strain_stats)
snf_total = sum(s['snf_total'] for s in strain_stats)
csv_total = sum(s['csv_total'] for s in strain_stats)

print(f"\n  Q1: How many single-caller SVs overlap between callers?")
print(f"  {'─'*70}")
print(f"  Sniffles2 total: {snf_total}")
print(f"  cuteSV total: {csv_total}")
print(f"  Sniffles2 calls found in cuteSV: {snf_in_csv_total}/{snf_total} = {snf_in_csv_total/snf_total*100:.0f}%")
print(f"  cuteSV calls found in Sniffles2: {csv_in_snf_total}/{csv_total} = {csv_in_snf_total/csv_total*100:.0f}%")
print(f"  → Only {snf_in_csv_total/snf_total*100:.0f}% of Sniffles2 calls are confirmed by cuteSV")
print(f"  → {(1-csv_in_snf_total/csv_total)*100:.0f}% of cuteSV calls are NOT in Sniffles2 (unique calls)")

# Q2: ICB vs. single callers
icb_total = sum(s['icb_total'] for s in strain_stats)
icb_in_snf_total = sum(s['icb_in_snf'] for s in strain_stats)
icb_in_csv_total = sum(s['icb_in_csv'] for s in strain_stats)
icb_in_both_total = sum(s['icb_in_both'] for s in strain_stats)
icb_in_one_total = sum(s['icb_in_one_only'] for s in strain_stats)

print(f"\n  Q2: How many ICB SVs are found by both callers vs. only one?")
print(f"  {'─'*70}")
print(f"  ICB total: {icb_total}")
print(f"  ICB found in Sniffles2: {icb_in_snf_total}/{icb_total} = {icb_in_snf_total/icb_total*100:.0f}%")
print(f"  ICB found in cuteSV: {icb_in_csv_total}/{icb_total} = {icb_in_csv_total/icb_total*100:.0f}%")
print(f"  ICB found in BOTH: {icb_in_both_total}/{icb_total} = {icb_in_both_total/icb_total*100:.0f}%")
print(f"  ICB found in ONE only: {icb_in_one_total}/{icb_total} = {icb_in_one_total/icb_total*100:.0f}%")
print(f"  → {icb_in_both_total/icb_total*100:.0f}% of ICB SVs are confirmed by both callers independently")
print(f"  → {icb_in_one_total/icb_total*100:.0f}% are rescued from a single caller by ICB matching")

# Q3: Are HIGH-confidence SVs enriched for dual-caller support?
print(f"\n  Q3: Are HIGH-confidence SVs enriched for dual-caller support?")
print(f"  {'─'*70}")
# We need to load validation data for this
import json
dual_high = 0
dual_total = 0
single_high = 0
single_total = 0

for folder in sorted(os.listdir('results/baumannii_strains/')):
    base = f'results/baumannii_strains/{folder}'
    val_file = f'{base}/validation/validation_results.json'
    icb_file = f'{base}/variants/consensus_svs.vcf'
    snf_file = f'{base}/variants/sniffles2_svs.vcf'
    csv_file = f'{base}/variants/cutesv_svs.vcf'
    
    if not all(os.path.exists(f) for f in [val_file, icb_file, snf_file, csv_file]):
        continue
    
    icb_set = parse_vcf_simple(icb_file)
    snf_set = parse_vcf_simple(snf_file)
    csv_set = parse_vcf_simple(csv_file)
    
    with open(val_file) as f:
        val_data = json.load(f)
    
    for r in val_data['results']:
        # Find matching ICB SV
        sv_key = (r['sv_chrom'], r['sv_start'], r['sv_end'], r['sv_type'])
        in_snf = any(overlap(sv_key, s, 500) for s in snf_set)
        in_csv = any(overlap(sv_key, s, 500) for s in csv_set)
        in_both = in_snf and in_csv
        
        if in_both:
            dual_total += 1
            if r['t_score'] >= 0.6: dual_high += 1
        else:
            single_total += 1
            if r['t_score'] >= 0.6: single_high += 1

print(f"  Dual-caller SVs that are HIGH: {dual_high}/{dual_total} = {dual_high/dual_total*100:.0f}%")
print(f"  Single-caller SVs that are HIGH: {single_high}/{single_total} = {single_high/single_total*100:.0f}%" if single_total > 0 else "  N/A")
print(f"  → Dual-caller SVs are {dual_high/dual_total*100 - single_high/single_total*100:.0f}% more likely to score HIGH" if single_total > 0 else "")

# Per-strain detail
print(f"\n  PER-STRAIN DETAIL:")
print(f"  {'Strain':<12} {'Snf→cSV':<10} {'cSV→Snf':<10} {'ICB in BOTH':<12} {'ICB in ONE':<12}")
for s in strain_stats:
    print(f"  {s['name']:<12} {s['snf_in_csv']/s['snf_total']*100:>6.0f}%{'':>3} {s['csv_in_snf']/s['csv_total']*100:>6.0f}%{'':>3} {s['icb_in_both']:>6}/{s['icb_total']:<4} {s['icb_in_one_only']:>6}/{s['icb_total']:<4}")

print("\n" + "=" * 100 + "\n")

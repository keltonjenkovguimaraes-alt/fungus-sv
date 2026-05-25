#!/usr/bin/env python3
"""Compare FUNGUS-SV ICB consensus vs. individual callers."""

import json, os, re
from collections import defaultdict

# Strain data from our 12 runs
strains_data = {}
for folder in sorted(os.listdir('results/baumannii_strains/')):
    base = f'results/baumannii_strains/{folder}'
    name = folder.replace('A_baumanii_', '')
    
    icb_file = f'{base}/variants/consensus_svs.vcf'
    snf_file = f'{base}/variants/sniffles2_svs.vcf'
    csv_file = f'{base}/variants/cutesv_svs.vcf'
    val_file = f'{base}/validation/validation_results.json'
    
    if not all(os.path.exists(f) for f in [icb_file, snf_file, csv_file, val_file]):
        continue
    
    # Count raw calls
    def count_vcf(f):
        return sum(1 for l in open(f) if not l.startswith('#'))
    
    snf_count = count_vcf(snf_file)
    csv_count = count_vcf(csv_file)
    icb_count = count_vcf(icb_file)
    
    # Load validation (T-scores)
    with open(val_file) as f:
        val_data = json.load(f)
    
    icb_high = sum(1 for r in val_data['results'] if r['t_score'] >= 0.6)
    icb_total = len(val_data['results'])
    
    strains_data[name] = {
        'sniffles2': snf_count,
        'cutesv': csv_count,
        'icb': icb_count,
        'icb_high': icb_high,
        'icb_total': icb_total
    }

# Print comparison table
print("\n" + "=" * 95)
print("  SINGLE CALLER vs. ICB CONSENSUS COMPARISON")
print("=" * 95)
print(f"  {'Strain':<14} {'Sniffles2':<12} {'cuteSV':<12} {'ICB':<10} {'ICB HIGH':<10} {'Filter %':<10}")
print(f"  {'':14} {'(raw)':<12} {'(raw)':<12} {'(≥2 callers)':<10} {'(T≥0.6)':<10}")

totals = {'snf': 0, 'csv': 0, 'icb': 0, 'icb_high': 0}
for name in sorted(strains_data.keys()):
    d = strains_data[name]
    totals['snf'] += d['sniffles2']
    totals['csv'] += d['cutesv']
    totals['icb'] += d['icb']
    totals['icb_high'] += d['icb_high']
    filter_pct = (1 - d['icb']/max(d['sniffles2'], d['cutesv'])) * 100
    print(f"  {name:<14} {d['sniffles2']:<12} {d['cutesv']:<12} {d['icb']:<10} {d['icb_high']:<10} {filter_pct:.0f}%")

t = totals
raw_avg = (t['snf'] + t['csv']) / 2
print(f"  {'─'*90}")
print(f"  {'TOTAL':<14} {t['snf']:<12} {t['csv']:<12} {t['icb']:<10} {t['icb_high']:<10}")
print(f"\n  COMPARISON METRICS:")
print(f"  ─────────────────")
print(f"  Average raw calls per caller: {raw_avg/len(strains_data):.0f}")
print(f"  Average ICB consensus calls: {t['icb']/len(strains_data):.0f}")
print(f"  ICB filters: {t['icb']}/{raw_avg:.0f} = {t['icb']/raw_avg*100:.0f}% retention")
print(f"  ICB HIGH confidence: {t['icb_high']}/{t['icb']} = {t['icb_high']/t['icb']*100:.0f}% of consensus")
print(f"")
print(f"  KEY FINDING:")
print(f"  Raw calls (Sniffles2 + cuteSV): {t['snf'] + t['csv']} total")
print(f"  ICB consensus (≥2 callers): {t['icb']} — removes {(1-t['icb']/((t['snf']+t['csv'])/2))*100:.0f}% of single-caller calls")
print(f"  ICB HIGH confidence: {t['icb_high']} — {(t['icb_high']/(t['snf']+t['csv']))*100:.0f}% of original raw calls survive as HIGH")
print(f"")
print(f"  If you used Sniffles2 alone: {t['snf']} SVs to validate")
print(f"  If you used cuteSV alone: {t['csv']} SVs to validate")
print(f"  If you use FUNGUS-SV ICB + validation: {t['icb_high']} HIGH-confidence SVs")
print(f"  Reduction in validation burden: {(1 - t['icb_high']/max(t['snf'], t['csv']))*100:.0f}%")
print("=" * 95 + "\n")

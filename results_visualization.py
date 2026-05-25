#!/usr/bin/env python3
"""Generate visual summaries of FUNGUS-SV results for README."""

import json, os

# Collect all strain data
strains_data = {}
for folder in os.listdir('results/baumannii_strains/'):
    jp = f'results/baumannii_strains/{folder}/validation/validation_results.json'
    if os.path.exists(jp):
        with open(jp) as f:
            data = json.load(f)
        name = folder.replace('A_baumanii_', '')
        results = data['results']
        high = sum(1 for r in results if r['t_score'] >= 0.6)
        med = sum(1 for r in results if 0.4 <= r['t_score'] < 0.6)
        weak = sum(1 for r in results if r['t_score'] < 0.4)
        dels = sum(1 for r in results if r['sv_type'] == 'DEL')
        invs = sum(1 for r in results if r['sv_type'] == 'INV')
        dups = sum(1 for r in results if r['sv_type'] == 'DUP')
        strains_data[name] = {
            'total': len(results), 'high': high, 'med': med, 'weak': weak,
            'del': dels, 'inv': invs, 'dup': dups,
            'pct_high': high/len(results)*100 if results else 0
        }

# Sort by total SVs descending
sorted_strains = sorted(strains_data.items(), key=lambda x: x[1]['total'], reverse=True)

print("\n" + "=" * 85)
print("  FUNGUS-SV RESULTS VISUALIZATION")
print("=" * 85)

# 1. Bar chart of SVs per strain (ASCII)
print("\n  ┌─ SV DETECTION BY STRAIN ─────────────────────────────────────┐")
max_svs = max(d['total'] for _, d in sorted_strains)
for name, d in sorted_strains:
    bar_len = int(d['total'] / max_svs * 40)
    bar = '█' * bar_len
    high_bar = '▓' * int(d['high'] / max_svs * 40)
    print(f"  │ {name:<12} {d['total']:>4} │{high_bar:<40}│ {d['pct_high']:.0f}% HIGH")
print("  └──────────────────────────────────────────────────────────────┘")

# 2. Confidence tier distribution (pie chart in ASCII)
print("\n  ┌─ CONFIDENCE DISTRIBUTION (860 SVs total) ───────────────────┐")
total_high = sum(d['high'] for _, d in sorted_strains)
total_med = sum(d['med'] for _, d in sorted_strains)
total_weak = sum(d['weak'] for _, d in sorted_strains)
total_all = total_high + total_med + total_weak

print(f"  │  HIGH (T≥0.6):  {total_high:>4} SVs  {'█' * int(total_high/total_all*50)} {total_high/total_all*100:.0f}%")
print(f"  │  MED  (0.4-0.6): {total_med:>4} SVs  {'▒' * int(total_med/total_all*50)} {total_med/total_all*100:.0f}%")
print(f"  │  WEAK (T<0.4):  {total_weak:>4} SVs  {'░' * int(total_weak/total_all*50)} {total_weak/total_all*100:.0f}%")
print("  └──────────────────────────────────────────────────────────────┘")

# 3. SV type distribution
print("\n  ┌─ SV TYPE DISTRIBUTION ──────────────────────────────────────┐")
total_del = sum(d['del'] for _, d in sorted_strains)
total_inv = sum(d['inv'] for _, d in sorted_strains)
total_dup = sum(d['dup'] for _, d in sorted_strains)
print(f"  │  Deletions:    {total_del:>4}  {'█' * 40} {total_del/total_all*100:.0f}%")
print(f"  │  Inversions:   {total_inv:>4}  {'▓' * int(total_inv/total_all*50)} {total_inv/total_all*100:.0f}%")
print(f"  │  Duplications: {total_dup:>4}  {'▒' * int(total_dup/total_all*50)} {total_dup/total_all*100:.0f}%")
print("  └──────────────────────────────────────────────────────────────┘")

# 4. Size-stratified summary
print("\n  ┌─ SIZE-STRATIFIED PERFORMANCE ───────────────────────────────┐")
print(f"  │  {'50-100 bp':<15} → WEAK (T≈0.286) — breakpoint only")
print(f"  │  {'100-500 bp':<15} → HIGH (T≈0.64) — depth + breakpoint")
print(f"  │  {'Pattern':<15} 100% consistent across all 12 strains")
print("  └──────────────────────────────────────────────────────────────┘")

# 5. LAR results
print("\n  ┌─ LAR VALIDATION (39 SVs tested) ────────────────────────────┐")
print(f"  │  Confirmed:  23  {'█' * 30} 59%")
print(f"  │  Partial:    12  {'▓' * 16} 31%")
print(f"  │  Failed:      4  {'░' * 5} 10%")
print(f"  │  Any Support: 35/39 = 90%")
print("  └──────────────────────────────────────────────────────────────┘")

# 6. Cross-species comparison
print("\n  ┌─ CROSS-SPECIES vs. WITHIN-SPECIES ──────────────────────────┐")
print(f"  │  Cross-species (5 spp.):   34 SVs  {'█' * 2}")
print(f"  │  Within-species (12 str.): 860 SVs {'█' * 40}")
print(f"  │  Ratio: 25× more SVs within species")
print("  └──────────────────────────────────────────────────────────────┘")

# 7. Summary table for README
print("\n  ┌─ PUBLICATION-READY SUMMARY TABLE ───────────────────────────┐")
print(f"  │  {'Strain':<14} {'SVs':>5} {'HIGH':>6} {'MED':>6} {'WEAK':>6} {'%HIGH':>6}")
print(f"  │  {'─'*50}")
for name, d in sorted_strains[:12]:
    print(f"  │  {name:<14} {d['total']:>5} {d['high']:>6} {d['med']:>6} {d['weak']:>6} {d['pct_high']:>5.0f}%")
print(f"  │  {'─'*50}")
print(f"  │  {'TOTAL':<14} {total_all:>5} {total_high:>6} {total_med:>6} {total_weak:>6} {total_high/total_all*100:>5.0f}%")
print("  └──────────────────────────────────────────────────────────────┘")
print("\n" + "=" * 85)

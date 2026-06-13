# FUNGUS-SV v0.9.4 — Final Fixes & Cleanup

**Date:** 2026-06-09  
**Status:** All fixes applied, tested, committed to GitHub

---

## Summary of Changes

This session addressed fundamental code quality issues, removed inappropriate tools, eliminated dead code, and fixed environment reproducibility.

---

## 1. Removed Longshot — Replaced Ploidy Layer

### Problem
Longshot is a diploid SNV caller designed for human genomes. It models two haplotypes, calls heterozygous variants (0/1), and phases SNPs with HapCUT2. None of this applies to haploid fungi.

Using it just to count heterozygous positions was like hiring a structural engineer to check if a door is locked — massive overkill. It also took 30-60 minutes per run.

### Fix
Replaced with a fast BCFtools-based ploidy check:
- Samples 10 random 50 kb genomic windows via `bcftools mpileup`
- Counts heterozygous calls (0/1)
- Het rate <3% → strongly haploid (score 1.0)
- Het rate <7% → mostly haploid (score 0.8)
- Het rate ≥7% → likely diploid (score 0.3)
- Runs in ~1-2 minutes total, **once before the SV loop**, not per SV

### Files Changed
| File | Change |
|------|--------|
| `valid_sv/evidence/layer_ploidy.py` | Complete rewrite — removed `run_longshot()`, replaced with `analyze_ploidy()` using BCFtools |
| `valid_sv/run_validation.py` | Removed `run_longshot` import, ploidy runs once before SV loop, result cached for all SVs |

---

## 2. Fixed Environment Reproducibility

### Problem
`workflow/envs/sv_calling.yaml` had no Python version pin, used `cutesv=1.0.8`, and included `pbsv=2.11.0` (unused dead weight). This caused Python 3.13 to be installed, which is incompatible with cuteSV.

### Fix
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.7
  - sniffles=2.2
  - cutesv=1.0.11
  - svim=1.4.2
Pinned Python 3.7 (compatible with all three callers)

Updated cuteSV to 1.0.11 (tested working)

Removed pbsv (never implemented, dead weight)

Removed bcftools (not needed in sv_call; already in sv_valid)

3. Fixed sniffles2 Binary Name
Problem
The conda package sniffles=2.2 installs the binary as sniffles (not sniffles2). The ICB script called sniffles2, which didn't exist.

Fix
Changed icb.py to call sniffles instead of sniffles2 in the subprocess command.

File Changed
File	Change
fungus_sv/core/icb.py	'sniffles2': ['sniffles2', ...] → 'sniffles2': ['sniffles', ...]
4. Fixed VCF Parsing for NUMCALLERS
Problem
run_validation.py parsed SUPPORT= from VCF INFO field, but the ICB consensus VCF uses NUMCALLERS=. All SVs were being assigned icb_support=1.

Fix
Parser now checks for SUPPORT= first, then falls back to NUMCALLERS=.

File Changed
File	Change
valid_sv/run_validation.py	Added fallback: if not support_match: support_match = re.search(r'NUMCALLERS=(\d+)', info)
5. Deleted Dead Code
Removed Files
File	Why Deleted
fungus_sv/core/build_consensus.py	Duplicate consensus logic — icb.py has its own build_consensus()
fungus_sv/modules/local_assembly.py	Old assembly module superseded by valid_sv/evidence/layer_lar.py
Removed Code
Location	What
icb.py	pbsv error block ("pbsv is not yet implemented")
layer_ploidy.py	run_longshot() function (replaced by BCFtools method)
run_validation.py	run_longshot import, longshot VCF path logic
6. Ploidy Now Runs Once (Not Per SV)
Problem
The ploidy check was inside the SV loop, running analyze_ploidy() 285 times. With BCFtools, each call takes ~5 seconds → 285 × 5 = ~24 minutes.

Fix
Ploidy runs once before the loop, result is cached in ploidy_result, and reused for every SV.
Before
for i, sv in enumerate(svs_to_validate):
    ...
    ploidy_result = analyze_ploidy(bam_path, reference_path)  # 285 times!
After
# Once before the loop
ploidy_result = analyze_ploidy(bam_path, reference_path)

for i, sv in enumerate(svs_to_validate):
    ...
    # Reuse cached result
    layer_results.append(LayerResult("ploidy_confirmation", ploidy_result.evidence_score, ...))
Result: Pipeline went from ~30+ minutes to 2 minutes 41 seconds.

7. Cleaned Conda Environments
Before (6 environments, 2 dead)
base, sv_align, sv_call, sv_cutesv, sv_longshot, sv_lar, sv_valid
After (5 environments, all active)
base, sv_align, sv_call, sv_lar, sv_valid
Removed	Why
sv_cutesv	No longer needed — cuteSV 1.0.11 works in sv_call with Python 3.7
sv_longshot	Longshot removed from pipeline entirely
Verification Results
Pipeline Test (2026-06-09)

ICB Consensus: 285 SVs (240 three-caller, 45 two-caller)
  DEL: 254, INV: 12, DUP: 19

Validation: 2m41s
  Ploidy: het_rate=0.000, haploid=True
  Mean T-scores: 0.35 (50-100bp) → 0.77 (>5kb)

Tier Distribution:
  TRIPLE_TRIANGULATED: ~94
  DOUBLE_CONFIRMED: ~28
  SINGLE_LINE: ~74
  DUP_SPLIT_READ_CONFIRMED: ~19
  INV_SPLIT_READ_CONFIRMED: ~12
  CONTRADICTED: ~16
LAR Spot Check
SV	Type	Size	Verdict	Reads	Contigs
ICB_NC_001146.8_546729_223	DEL	15.3 kb	✅ CONFIRMED	366	3
ICB_NC_001224.1_11503_264	DUP	36.8 kb	Running...	—	—
What Was NOT Changed
Scorer logic (DUP/INV overrides intact)

Depth layer (DHBFC, size stratification intact)

Breakpoint layer

LAR layer (Flye + Miniasm+Racon)

Genomic context layer

k-mer layer (still present, still optional)

What Still Needs Work (Future)
Task	Priority
Run on second fungal species	High — generalizability
Complete synthetic FDR benchmark	High — quantitative validation
Add BND/translocation support	Medium
Integrate hifiasm as Tier 3 assembler	Medium (needs Colab/bigger machine)
Replace or remove k-mer layer	Low
Docker container	Low (after stabilization)
Commit
text
v0.9.4: Cleanup - removed longshot, deleted duplicate code, 
fixed env YAMLs, fast ploidy, NUMCALLERS fix
All changes pushed to GitHub: github.com/keltonjenkovguimaraes-alt/fungus-sv

*Generated: 2026-06-09*

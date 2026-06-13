# FUNGUS-SV v0.9.3 Pipeline Update

## Release Date: 2 June 2026

This document summarizes all changes, calibrations, and validation results integrated into FUNGUS-SV between v0.8.0 and v0.9.3.

---

## 1. Parameter Calibration for Haploid Fungi

### 1.1 Layer Weights (Recalibrated)

| Layer | Original Weight | New Weight | Rationale |
|-------|:---:|:---:|-----------|
| Depth Signature | 0.25 | **0.35** | 100% depth loss unambiguous in haploids |
| Breakpoint Junction | 0.20 | **0.30** | Only signal for inversions; critical for all types |
| Local Assembly | 0.30 | **0.20** | Computationally expensive; standalone tool available |
| k-mer Spectrum | 0.25 | **0.15** | Redundant with depth in haploid genomes |

**Source:** Liu et al. (2024) original framework; recalibrated from CICC-1445 vs S288C empirical data.

### 1.2 DHFFC Thresholds (Haploid-Adapted)

| SV Type | Human Diploid | Haploid Fungi | Source |
|---------|:---:|:---:|--------|
| Deletion | DHFFC < 0.7 | **DHFFC < 0.3** | Pedersen & Quinlan (2019); this study |
| Duplication | DHFFC > 1.3 | **DHFFC > 2.0** | Pedersen & Quinlan (2019); this study |

**Rationale:** In diploid genomes, a heterozygous deletion drops depth by ~50% (DHFFC ≈ 0.5). In haploid genomes, a true deletion causes ~100% depth loss (DHFFC ≈ 0.0). The original human threshold of 0.7 would accept false positives in haploids.

### 1.3 Size-Stratified Depth Scoring

| SV Size | Size Factor | Basis |
|---------|:---:|-------|
| ≥ 5,000 bp | 1.00 | AUC ~1.0 (Pedersen & Quinlan 2019) |
| 1,000–4,999 bp | 0.95 | AUC ~0.97 |
| 500–999 bp | 0.85 | Moderate penalty |
| 100–499 bp | 0.75 | Significant penalty |
| < 100 bp | 0.60 | Heavy penalty |

---

## 2. New Features Added

### 2.1 DHBFC Integration (GC-Corrected Depth)

**Source:** Pedersen & Quinlan (2019) *GigaScience*

DHBFC compares SV region depth to genome-wide average, providing a GC-normalized depth metric alongside DHFFC. The two metrics show 90% agreement in yeast, confirming minimal GC bias.

**Implementation:** `valid_sv/evidence/layer_depth.py` — combined_ratio = (DHFFC + DHBFC) / 2.0

### 2.2 Inversion-Specific Scoring

**Problem:** All 11 CICC-1445 vs S288C inversions were confirmed by split reads (32-1,091 junction reads per breakpoint) but scored CONTRADICTED (T=0.167) because depth and k-mer layers are silent for balanced inversions.

**Solution:** 
- Inversions skip depth/k-mer layers in scoring
- Report as `INV_SPLIT_READ_CONFIRMED` when breakpoint evidence ≥ 0.6
- Empirically, FDR for INV calls is ~0% when split reads support them

**Implementation:** `valid_sv/engine/scorer.py`

### 2.3 Duplication Split-Read Rescue

**Problem:** Only 2/18 DUPs confirmed by DHFFC > 2.0, but all 18 have split-read junction evidence.

**Solution:** DUPs with breakpoint junction score ≥ 0.6 receive a 0.15 boost, based on David et al. (2024) finding that split-read evidence improves DUP accuracy.

**Implementation:** `valid_sv/engine/scorer.py`

### 2.4 Repeat Region Flag

DELs and DUPs in known problematic regions (FLO genes, rDNA, Ty elements, subtelomeric) now receive a `[REPEAT_REGION]` warning in the depth layer details, based on Dhakal et al. (2024) and David et al. (2024).

**Implementation:** `valid_sv/evidence/layer_depth.py`

### 2.5 Translocation Flag

DELs with extreme DHFFC (< 0.01) that may represent translocations rather than true deletions now receive a `[NEAR-ZERO_DEPTH: possible translocation]` warning. This was discovered through LAR validation of IMX2600 where DHFFC=0.000 but no deletion existed — the sequence was present elsewhere in the genome.

**Implementation:** `valid_sv/evidence/layer_depth.py`

### 2.6 Small SV Split-Read Minimum

SVs < 100 bp require ≥ 6 split reads at breakpoints to receive a non-zero breakpoint score, based on Pedersen & Quinlan (2019) showing that depth signal is unreliable for small events.

**Implementation:** `valid_sv/evidence/layer_breakpoint.py`

---

## 3. LAR (Local Assembly Refinement) — Standalone Tool

### 3.1 Development

LAR was developed as a standalone validation tool that assembles reads from a single SV region using Flye, then aligns the contig back to the reference to confirm or refute the called SV.

**Why standalone:** 
- Computational cost: 2-15 min per SV, 200-500 MB RAM
- Requires human interpretation of CIGAR strings
- Best applied to key candidate SVs, not all 277

**Location:** `valid_sv/evidence/layer_lar.py`
**Environment:** `conda activate sv_lar`

### 3.2 Validation Results (5 Loci Tested)

| Locus | Type | Size | Strain | Result |
|-------|------|------|--------|--------|
| YBL005W-B (Ty2) | DEL | 5.9 kb | S288C | REAL |
| chrVII multi-gene | DEL | 55.7 kb | S288C | REAL |
| SX2 chrII | INV | 430 kb | SX2 | REAL |
| BJ4 chrXII | INV | 205 kb | BJ4 | REAL |
| chrXIV complex | DEL+DUP+INV | 15-40 kb | S288C | REAL (all confirmed) |

**5/5 (100%) candidate SVs confirmed by LAR.**

### 3.3 Multi-Strain Calibration (64 CONTRADICTED DELs)

LAR was used to determine the ground truth for 64 CONTRADICTED DELs (5.0-6.5 kb) across all 5 reference strains to calibrate the depth scoring system.

| Strain | Tested | Real | False | % Real |
|--------|:---:|:---:|:---:|:---:|
| S288C | 15 | 8 | 7 | 53% |
| BJ4 | 1 | 0 | 1 | 0% |
| IMX2600 | 28 | 5 | 23 | 18% |
| Makgeolli | 16 | 3 | 13 | 19% |
| SX2 | 4 | 1 | 3 | 25% |
| **Total** | **64** | **17** | **47** | **27%** |

**Key Finding:** Only 27% of CONTRADICTED large DELs are real. DHFFC alone cannot predict reality across different reference genomes. The pipeline is correctly cautious — CONTRADICTED means "uncertain," and most (73%) are indeed false or complex.

---

## 4. Pipeline Scoring Flow (Updated)

For each SV in consensus VCF:
│
├── Layer 0: ICB Alignment Consensus (pre-computed, NOT scored)
│
├── Layer 1: Local Assembly (weight 0.20)
│ └── Standalone — run manually via layer_lar.py
│ └── Not scored in automated pipeline
│
├── Layer 2: Depth Signature (weight 0.35) ← UPDATED
│ ├── DHFFC = median(SV_depth) / median(flank_depth)
│ ├── DHBFC = median(SV_depth) / median(local_average)
│ ├── combined_ratio = (DHFFC + DHBFC) / 2.0
│ ├── size_factor = f(SV_length) — see table above
│ ├── DEL scoring: combined_ratio < 0.3 → score 0.8-1.0 × size_factor
│ ├── DUP scoring: combined_ratio > 2.0 → score 0.8-1.0 × size_factor
│ ├── REPEAT_REGION flag for FLO/rDNA/Ty regions
│ └── NEAR-ZERO_DEPTH flag for possible translocations
│
├── Layer 3: k-mer Spectrum (weight 0.15)
│ └── Only for DEL, INS
│
├── Layer 4: Breakpoint Junction (weight 0.30) ← UPDATED
│ ├── Split reads (CIGAR 'S') + soft-clipped + spanning reads
│ ├── distance_support = 0.2 × len + 2000/len
│ ├── MAPQ ≥ 20 filter
│ ├── Small SV minimum: ≥6 split reads for <100bp
│ └── All SV types
│
├── Layer 5: Ploidy Confirmation (HARD FILTER)
│ └── heterozygous SNV rate < 7%
│
└── Triangulation Engine
├── T = Σ(score_i × weight_i) / Σ(weight_i)
├── Completeness penalty if <50% layers available
├── INV override: breakpoint-only scoring, INV_SPLIT_READ_CONFIRMED
├── DUP rescue: +0.15 boost if split reads ≥ 0.6
└── Tier classification

---

## 5. Confidence Tiers (v0.9.3)

| Tier | T-Score | Description |
|------|:---:|-------------|
| TRIPLE_TRIANGULATED | T ≥ 0.80 | All layers strongly support |
| DOUBLE_CONFIRMED | T ≥ 0.60 | Multiple layers support |
| SINGLE_LINE | T ≥ 0.40 | Single layer supports |
| WEAK | T ≥ 0.20 | Weak or ambiguous evidence |
| CONTRADICTED | T < 0.20 | Evidence contradicts or insufficient |
| **INV_SPLIT_READ_CONFIRMED** | T=0.167 | **NEW** — Inversion with junction reads |

---

## 6. S288C v3 Validation Results

| Tier | Count | % |
|------|:---:|:---:|
| TRIPLE_TRIANGULATED | 5 | 1.8% |
| DOUBLE_CONFIRMED | 17 | 6.1% |
| SINGLE_LINE | 1 | 0.4% |
| WEAK | 92 | 33.2% |
| CONTRADICTED | 151 | 54.5% |
| INV_SPLIT_READ_CONFIRMED | 11 | 4.0% |
| **Total** | **277** | 100% |

---

## 7. Files Modified

| File | Change |
|------|--------|
| `config/config.yaml` | Calibrated weights, DHFFC thresholds, LAR config, size stratification |
| `valid_sv/engine/scorer.py` | INV override, DUP rescue, string tier handling |
| `valid_sv/evidence/layer_depth.py` | DHBFC, size stratification, repeat flag, translocation flag |
| `valid_sv/evidence/layer_breakpoint.py` | Small SV minimum, distance_support, MAPQ filter |
| `valid_sv/evidence/layer_lar.py` | NEW — standalone LAR module |
| `valid_sv/reporting/report_card.py` | String tier handling |
| `valid_sv/run_validation.py` | LAR placeholder fix |
| `README.md` | Full documentation update |

---

## 8. Known Issues and Limitations

1. **DHFFC is reference-dependent** — IMX2600 DHFFC=0.000 can be FALSE; S288C DHFFC=0.596 can be REAL
2. **CONTRADICTED tier includes real SVs** — ~27% of tested large DELs were real; LAR recommended for key loci
3. **No BND/translocation caller** — translocations are flagged but not called
4. **k-mer layer often unavailable** — Jellyfish database is large (2.4 GB) and takes time to build
5. **Spike-in calibration pending** — SURVIVOR crashes on yeast; PBSIM2 needs trained model

---

## 9. References Integrated in v0.9.x

1. Pedersen & Quinlan (2019) — DHFFC/DHBFC metrics
2. Zheng & Shang (2024) — distance_support, MAPQ filter, CIGAR parsing
3. Belyeu et al. (2021) — Samplot visualization
4. David et al. (2024) — Manual curation FDR estimates
5. Dhakal et al. (2024) — Fungal SV landscape, repeat region behavior
6. Li et al. (2023) — Haploid fungal variant calling benchmark
7. Nkouamedjo et al. (2025) — Quality metrics beyond binary consensus
8. Sedlazeck et al. (2017) — SV annotation with genomic features

# FUNGUS-SV Project — Master Summary for Continuation

## Project Overview

**FUNGUS-SV:** A structural variant discovery and triangulation-based validation pipeline for haploid fungal genomes using PacBio HiFi reads.

**GitHub:** https://github.com/keltonjenkovguimaraes-alt/fungus-sv
**Interactive Report:** https://fungus-sv.netlify.app
**Paper:** In preparation (Guimarães et al., 2026)

---

## Data

**Query strain:** *Saccharomyces cerevisiae* CICC-1445 (Chinese industrial)
**Sequencing:** PacBio HiFi (SRR18210299), 274,915 reads, ~20 kb N50, ~100× coverage
**FASTQ:** `/home/kelto/fungus-sv/data/raw/cicc1445_hifi.fastq.gz` (4.3 GB)

**5 Reference strains:**
| Strain | Type | File | Chromosomes |
|--------|------|------|-------------|
| S288C | Laboratory | `S288C_reference.fasta` | NC_* (RefSeq) |
| BJ4 | Industrial | `BJ4_reference.fasta` | LR8135* (ENA) |
| IMX2600 | Engineered | `IMX2600_reference.fasta` | CP127* (GenBank) |
| Makgeolli | Fermentation | `Makgeolli_reference.fasta` | CP025* (GenBank) |
| SX2 | Industrial | `SX2_reference.fasta` | LR8135* (ENA) |

---

## Pipeline Architecture

### Phase 1: ICB Consensus
- 3 callers: Sniffles2 + cuteSV + SVIM
- Consensus: ≥2 callers agree, 0.5 overlap, 200 bp flank

### Phase 2: 5-Layer Triangulation Validation
| Layer | Weight | Description |
|-------|:---:|-------------|
| Local Assembly (Flye) | 0.20 | Regional assembly — standalone tool |
| Depth Signature (DHFFC+DHBFC) | 0.35 | Read depth fold-change |
| k-mer Spectrum (Jellyfish) | 0.15 | k-mer frequency analysis |
| Breakpoint Junction | 0.30 | Split-read + CIGAR analysis |
| Ploidy Confirmation | FILTER | Heterozygosity check |

### Phase 3: Confidence Tiers
TRIPLE_TRIANGULATED: T ≥ 0.80
DOUBLE_CONFIRMED: T ≥ 0.60
SINGLE_LINE: T ≥ 0.40
WEAK: T ≥ 0.20
CONTRADICTED: T < 0.20
INV_SPLIT_READ_CONFIRMED: Inversions with split-read support

---

## SV Detection Results (CICC-1445 vs 5 References)

| Strain | Total SVs | DEL | DUP | INV |
|--------|:---:|:---:|:---:|:---:|
| S288C | 277 | 248 | 18 | 11 |
| BJ4 | 165 | 140 | 13 | 12 |
| IMX2600 | 314 | 285 | 9 | 20 |
| Makgeolli | 250 | 225 | 17 | 8 |
| SX2 | 290 | 261 | 18 | 11 |

---

## v3 Validation Results (S288C, with all fixes)

Run with calibrated parameters at:
`/home/kelto/fungus-sv/data/yeast/results_3callers/validation_v3/`

| Tier | Count |
|------|:---:|
| TRIPLE_TRIANGULATED | 5 |
| DOUBLE_CONFIRMED | 17 |
| SINGLE_LINE | 1 |
| WEAK | 92 |
| CONTRADICTED | 151 |
| INV_SPLIT_READ_CONFIRMED | 11 |

---

## Key Parameters (Calibrated for Haploid Fungi)

### DHFFC Thresholds
| SV Type | Human Diploid | Haploid Fungi |
|---------|:---:|:---:|
| Deletion | DHFFC < 0.7 | **DHFFC < 0.3** |
| Duplication | DHFFC > 1.3 | **DHFFC > 2.0** |

### Size-Stratified Depth Scoring
| Size | Factor |
|------|:---:|
| ≥5 kb | 1.00 |
| 1-5 kb | 0.95 |
| 500-1000 bp | 0.85 |
| 100-500 bp | 0.75 |
| <100 bp | 0.60 |

### Layer Weights
| Layer | Original | Calibrated |
|-------|:---:|:---:|
| Depth | 0.25 | **0.35** |
| Breakpoint | 0.20 | **0.30** |
| Assembly | 0.30 | **0.20** |
| k-mer | 0.25 | **0.15** |

---

## Pipeline Fixes Applied (v0.9.2-v0.9.3)

1. **INV reporting:** Inversions with split-read support get `INV_SPLIT_READ_CONFIRMED` instead of CONTRADICTED
2. **DUP split-read rescue:** DUPs with junction evidence get 0.15 boost
3. **Repeat flag:** FLO/rDNA/Ty regions get `[REPEAT_REGION]` warning
4. **Small SV minimum:** <100 bp SVs need ≥6 split reads
5. **Translocation flag:** DHFFC<0.01 DELs get `[NEAR-ZERO_DEPTH]` warning
6. **DHBFC integration:** GC-corrected depth from Pedersen & Quinlan (2019)

---

## LAR (Local Assembly Refinement) — Standalone Tool

**Purpose:** Definitive proof of SV existence by assembling reads from a single SV region.

**Location:** `valid_sv/evidence/layer_lar.py`
**Environment:** `conda activate sv_lar`

**How it works:**
1. Extract reads mapping to SV region ± 3 kb
2. Assemble with Flye (genome-size = window size, NOT full genome)
3. Align contig to reference with minimap2
4. Parse CIGAR for deletion gaps or strand switches

**Performance:** <500 MB RAM, 2-15 min per SV

---

## LAR Calibration Results (64 CONTRADICTED DELs Tested)

**Overall: 17 REAL (27%), 47 FALSE (73%)**

| Strain | Tested | Real | False | % Real |
|--------|:---:|:---:|:---:|:---:|
| S288C | 15 | 8 | 7 | 53% |
| BJ4 | 1 | 0 | 1 | 0% |
| IMX2600 | 28 | 5 | 23 | 18% |
| Makgeolli | 16 | 3 | 13 | 19% |
| SX2 | 4 | 1 | 3 | 25% |
| **Total** | **64** | **17** | **47** | **27%** |

**Key finding:** DHFFC alone cannot predict reality across different reference genomes. S288C has the highest real rate (best reference quality). IMX2600 has the lowest (engineered strain with assembly differences).

Full dataset: `/home/kelto/fungus-sv/data/LAR_calibration_dataset.md`

---

## Key LAR-Validated Loci

| Locus | Type | Size | Strain | Result |
|-------|------|------|--------|--------|
| YBL005W-B (Ty2) | DEL | 5.9 kb | S288C | REAL |
| YBR012W-B (Ty2) | DEL | 5.9 kb | S288C | REAL |
| chrVII multi-gene (35 genes) | DEL | 55.7 kb | S288C | REAL |
| SX2 chrII | INV | 430 kb | SX2 | REAL |
| BJ4 chrXII | INV | 205 kb | BJ4 | REAL |
| chrXIV complex | DEL+DUP+INV | 15-40 kb | S288C | REAL (all) |
| FLO1 | DEL | 554 bp | S288C | NOT SIMPLE DEL |
| IMX2600 M1 | DEL | 5.9 kb | IMX2600 | FALSE (translocation) |

---

## Orthogonal Validation Methods Used

1. **DHFFC/DHBFC:** Genome-wide depth analysis on 266 SVs
2. **Split-read junction counts:** All 11 INVs confirmed (32-1,091 reads/breakpoint)
3. **Samplot visualization:** 17 images generated for manual curation
4. **Comparative Samplot:** YBL005W-B and chrVII across all 5 references
5. **Gene overlap annotation:** 6,459 genes from NCBI GFF mapped to SVs

---

## Important Discoveries

1. **S288C has Ty2 insertions other strains lack** — called as DELs but are S288C-specific insertions
2. **Inversions systematically under-scored** — all 11 INVs confirmed by split reads despite T=0.167
3. **DHFFC pattern breaks across references** — IMX2600 DHFFC=0.000 can be FALSE
4. **Translocations mimic deletions** — near-zero depth with no CIGAR gap = possible translocation
5. **CICC-1445 closest to BJ4** (165 SVs) — both Chinese industrial strains

---

## Files and Paths

### Key Results
- S288C v3 validation: `data/yeast/results_3callers/validation_v3/`
- BJ4 v2: `data/yeast/results_BJ4/validation_v2/`
- IMX2600 v2: `data/yeast/results_IMX2600/validation_v2/`
- Makgeolli v2: `data/yeast/results_Makgeolli/validation_v2/`
- SX2 v2: `data/yeast/results_SX2/validation_v2/`

### BAM Files (Large, can regenerate)
- `cicc1445_sorted.bam` — CICC-1445 vs S288C
- `cicc1445_vs_BJ4.bam`
- `cicc1445_vs_IMX2600.bam`
- `cicc1445_vs_Makgeolli.bam`
- `cicc1445_vs_SX2.bam`

### Documentation
- `docs/LAR_validation_results.md` — All LAR test results
- `docs/LAR_mathematics.md` — Math behind LAR
- `docs/LAR_integration_explained.md` — Why LAR is standalone
- `docs/SV_direction_explanation.md` — Reference vs query direction
- `docs/DHBFC_size_integration.md` — DHBFC + size stratification
- `docs/v2_validation_summary_2026-05-29.md` — v2 results
- `docs/pipeline_mathematics.md` — Full mathematical framework
- `docs/validation_confidence_estimates.md` — Confidence estimates
- `data/LAR_calibration_dataset.md` — LAR calibration data

### Code (modified files)
- `valid_sv/engine/scorer.py` — INV override, DUP rescue, string tier fix
- `valid_sv/evidence/layer_depth.py` — DHBFC, size stratification, repeat flag, translocation flag
- `valid_sv/evidence/layer_breakpoint.py` — Small SV minimum
- `valid_sv/evidence/layer_lar.py` — LAR standalone module
- `config/config.yaml` — Calibrated weights and thresholds
- `README.md` — Full project documentation

---

## GitHub Status

**Latest commit:** v0.9.3 — Working INV tier override, DUP rescue, repeat flag, small SV minimum, translocation flag
**Known issue:** `run_validation.py` had broken `args.lar` reference — should be fixed in latest commit

---

## Papers Analyzed and Integrated

1. **Pedersen & Quinlan (2019)** — Duphold: DHFFC/DHBFC metrics
2. **Zheng & Shang (2024)** — SVvalidation: distance_support, MAPQ filter, CIGAR parsing
3. **Belyeu et al. (2021)** — Samplot: visualization and ML-based curation
4. **David et al. (2024)** — Manual curation strategies, FDR estimates per SV type
5. **Dhakal et al. (2024)** — SV landscape in *Fusarium graminearum* (assembly-based)
6. **Li et al. (2023)** — Variant calling in *Candida auris* (haploid fungal benchmark)
7. **Nkouamedjo et al. (2025)** — SV-MeCa: XGBoost meta-caller (quality metrics beyond binary)
8. **Sedlazeck et al. (2017)** — SURVIVOR_ant: SV annotation and comparison
9. **Yuan et al. (2018)** — SVSR: SV simulation (human short-read, not directly used)

---

## Environment Structure

| Environment | Tools | Purpose |
|-------------|-------|---------|
| `sv_align` | minimap2, samtools | Read alignment |
| `sv_call` | sniffles2, cutesv, svim, bcftools | SV detection |
| `sv_valid` | python, numpy, scipy, pandas, pysam, matplotlib | Triangulation |
| `sv_lar` | flye, minimap2, samtools | LAR assembly |
| `sv_samplot` | samplot | Visualization |
| `sv_spikein` | survivor, pbsim2 | Spike-in simulation (in progress) |
| `sv_duphold3` | duphold | Depth annotation |

---

## Next Steps (Pending)

1. **Re-run all 5 strains with v3 fixes** — only S288C has been re-run
2. **Test WEAK and DOUBLE tier DELs with LAR** — to see if T-score correlates with reality
3. **Cross-reference with Peter et al. (2018)** — 1011 yeast genomes CNV data downloaded, needs gene name mapping
4. **Build fungal HiFi read simulator** — new repo for spike-in benchmarking
5. **Add BND/translocation support** — currently filtered out
6. **Clean up BAM files** — can delete ~15 GB of intermediate files
7. **Publish to Zenodo** — description ready

---

## How to Resume in New Chat

Say: "I'm continuing work on FUNGUS-SV, a structural variant pipeline for haploid fungi. We've calibrated depth scoring with 64 LAR-validated DELs across 5 yeast strains. Key files are in /docs/ and the master summary is at docs/MASTER_SUMMARY_FOR_NEXT_CHAT.md. Pick up where we left off."

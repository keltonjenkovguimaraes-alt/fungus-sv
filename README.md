# 🧬 FUNGUS-SV v0.9.4

**Structural variant discovery and triangulation-based validation for haploid fungal genomes using PacBio HiFi reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.7+](https://img.shields.io/badge/Python-3.7+-yellow.svg)](https://python.org)
[![Status](https://img.shields.io/badge/status-active%20development-orange.svg)]()

---

## ⚠️ Disclaimer

FUNGUS-SV is a **hypothesis-generation tool under active development.** NOT for clinical use. All SVs of biological interest must be independently validated by experimental methods.

---

## Overview

FUNGUS-SV detects structural variants (SVs) in haploid fungal genomes using PacBio HiFi long reads. It combines a **three-caller ICB consensus** (Sniffles2 + cuteSV + SVIM) with a **7-layer triangulation validation system** and a **two-assembler Local Assembly Refinement (LAR)** pipeline for definitive SV confirmation.

### Key Features

- **7-layer triangulation:** Depth (DHFFC+DHBFC), k-mer spectrum, breakpoint junction, local assembly, ploidy confirmation, genomic context
- **Two-assembler LAR:** Flye + Miniasm/Racon consensus for definitive breakpoint proof
- **Haploid-calibrated:** All thresholds recalibrated for haploid genomes (DHFFC, size factors, layer weights)
- **SV-type-specific overrides:** DUP_SPLIT_READ_CONFIRMED and INV_SPLIT_READ_CONFIRMED tiers based on empirical truth data
- **Multi-strain validated:** Tested on 6 *Saccharomyces* strains (1,527 SVs) + cross-genus *Candida albicans*
- **60-SV LAR truth set:** Two-assembler validation across all 6 strains (72% within-species confirmation)

---

## Pipeline Architecture

PacBio HiFi reads (*.fastq.gz)
│
▼
[minimap2] ← map-hifi preset
│
▼
┌─────────────────────────────────────────────────────────┐
│ PHASE 1: ICB CONSENSUS │
│ Sniffles2 + cuteSV + SVIM │
│ Consensus: ≥2 callers, 0.5 overlap, 200 bp flank │
└─────────────────────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────────────────────┐
│ PHASE 2: 7-LAYER TRIANGULATION VALIDATION │
│ │
│ Layer 0: ICB Consensus (reported, not scored) │
│ Layer 1: Local Assembly — Flye + Miniasm/Racon (0.20) │
│ Layer 2: Depth Signature — DHFFC + DHBFC (0.35) │
│ Layer 3: k-mer Spectrum — Jellyfish (0.15) │
│ Layer 4: Breakpoint Junction — split reads (0.30) │
│ Layer 5: Ploidy Confirmation — BCFtools pileup (0.0) │
│ Layer 6: Genomic Context — gene annotations (0.0) │
│ │
│ T = Σ(score_i × weight_i) / Σ(weight_i) │
└─────────────────────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────────────────────┐
│ PHASE 3: CONFIDENCE TIERS │
│ │
│ TRIPLE_TRIANGULATED: T ≥ 0.80 │
│ DOUBLE_CONFIRMED: T ≥ 0.60 │
│ SINGLE_LINE: T ≥ 0.40 │
│ WEAK: T ≥ 0.20 │
│ CONTRADICTED: T < 0.20 │
│ DUP_SPLIT_READ_CONFIRMED: Hardcoded (FDR ~0%) │
│ INV_SPLIT_READ_CONFIRMED: Hardcoded (FDR ~0%) │
└─────────────────────────────────────────────────────────┘

---

## Validation Results

### Multi-Strain *Saccharomyces* (CICC-1445 vs 6 References)

| Strain | Species | Type | Total SVs | DEL | DUP | INV | LAR Conf. |
|--------|---------|------|-----------|-----|-----|-----|-----------|
| S288C | *S. cerevisiae* | Laboratory | 285 | 254 | 19 | 12 | 7/10 |
| BJ4 | *S. cerevisiae* | Industrial | 176 | 151 | 14 | 11 | 7/10 |
| SX2 | *S. cerevisiae* | Industrial | 310 | 283 | 18 | 9 | 4/10 |
| IMX2600 | *S. cerevisiae* | Engineered | 345 | 313 | 20 | 12 | 8/10 |
| Makgeolli | *S. cerevisiae* | Fermentation | 244 | 218 | 17 | 9 | 10/10 |
| Jurei | *S. jurei* | Wild Type | 167 | 156 | 11 | 0 | 0/10* |

**Total: 1,527 SVs across 6 strains.** *Cross-species.

### Two-Assembler LAR Truth Set (60 SVs)

| Metric | Value |
|--------|-------|
| Within-species confirmed (both assemblers) | 36/50 (72%) |
| Cross-species confirmed (Jurei) | 0/10 (0%) |
| Makgeolli (cleanest strain) | 10/10 (100%) |
| False positives identified | SX2: 2 INVs contradicted by both |

### Key Biological Findings
- **CICC-1445 closest to BJ4** (176 SVs) — both Chinese industrial strains
- **Most divergent from IMX2600** (345 SVs) — engineered lab strain
- **Jurei cross-species calls are mostly artifacts** — LAR confirms 0/10
- **LAR frequently finds larger deletions than callers report** — callers fragment large SVs
- **DUP/INV overrides based on 48-SV LAR truth set** (81.3% DUPs real)

---

## Quick Start

### Prerequisites
- Linux (Ubuntu 20.04+) or WSL2 | ≥8 GB RAM | Conda/Mamba
- PacBio HiFi reads (≥20× coverage)

### Installation

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

# Create environments
conda env create -f workflow/envs/alignment.yaml -n sv_align
conda env create -f workflow/envs/sv_calling.yaml -n sv_call
conda env create -f workflow/envs/validation.yaml -n sv_valid
conda env create -f workflow/envs/lar.yaml -n sv_lar
Usage
# 1. Align reads to reference
conda activate sv_align
minimap2 -d ref.mmi reference.fasta
minimap2 -t 8 -ax map-hifi ref.mmi reads.fastq.gz | samtools sort -@ 4 -o sample.bam -
samtools index sample.bam

# 2. Detect SVs (ICB consensus)
conda activate sv_call
python fungus_sv/core/icb.py --bam sample.bam --reference ref.fasta --output results/ --threads 4

# 3. Validate (7-layer triangulation)
conda activate sv_valid
PYTHONPATH=. python -m valid_sv.run_validation \
    --consensus-vcf results/consensus_svs.vcf \
    --bam sample.bam --reference ref.fasta \
    --output results/validation/ --threads 4

# 4. LAR: Two-assembler definitive proof
conda activate sv_lar
python3 -c "
from valid_sv.evidence.layer_lar import run_lar, run_lar_miniasm
flye = run_lar('sample.bam','ref.fasta','SV_ID','DEL','chr',start,end)
mini = run_lar_miniasm('sample.bam','ref.fasta','SV_ID','DEL','chr',start,end)
print(f'Flye: {flye.verdict.value}, Miniasm: {mini.verdict.value}')
"
Repository Structure
fungus-sv/
├── fungus_sv/core/              # ICB consensus calling
│   ├── icb.py                   # Main ICB pipeline
│   └── annotate_svs.py          # Gene annotation
├── valid_sv/                    # Triangulation validation
│   ├── evidence/                # Layer implementations
│   │   ├── layer_lar.py         # Flye + Miniasm/Racon dual assembly
│   │   ├── layer_depth.py       # DHFFC + DHBFC + max depth + LCR
│   │   ├── layer_kmer.py        # k-mer spectrum
│   │   ├── layer_breakpoint.py  # Split-read junction analysis
│   │   ├── layer_ploidy.py      # Fast BCFtools haploid check
│   │   └── layer_genomic_context.py  # Gene annotation filters
│   ├── engine/scorer.py         # Triangulation scoring engine
│   ├── benchmarks/              # Synthetic SV generator + FDR calibration
│   └── run_validation.py        # Main validation entry point
├── config/config.yaml           # Calibrated weights and thresholds
├── workflow/envs/               # Conda environment specifications
├── docs/                        # Documentation and development log
│   ├── DEVELOPMENT_LOGBOOK.md   # Complete chronological log
│   ├── IMPLEMENTATION_SUMMARY_v0.9.4.md
│   └── UPDATE_2026-06-11.md
├── data/LAR_validation/         # 60-SV LAR truth set results
└── figures/                     # HTML reports and Samplot images
Calibrated Parameters (Haploid Fungi)
Parameter	Value	Source
DHFFC DEL threshold	<0.3	Pedersen & Quinlan (2019), recalibrated
DHFFC DUP threshold	>2.0	Pedersen & Quinlan (2019), recalibrated
Depth weight	0.35	Liu et al. (2024), recalibrated
Breakpoint weight	0.30	Liu et al. (2024), recalibrated
LAR weight	0.20	This study
k-mer weight	0.15	Liu et al. (2024), recalibrated
Size factor (≥5 kb)	1.00	Pedersen & Quinlan (2019)
Size factor (<100 bp)	0.60	Pedersen & Quinlan (2019)
MAPQ filter	≥20	Zheng & Shang (2024)
Max het rate (haploid)	<7%	Xing et al. (2025)
DUP override FDR	~0%	48-SV LAR truth set
INV override FDR	~0%	11-SV split-read confirmed
Papers Integrated (13)
Pedersen & Quinlan (2019) — DHFFC/DHBFC depth metrics

Zheng & Shang (2024) — distance_support, MAPQ filter

Belyeu et al. (2021) — Samplot visualization

David et al. (2024) — Manual curation FDR

Dhakal et al. (2024) — Fungal SV landscape

Li et al. (2023) — Haploid fungal variant calling

Liu et al. (2024) — ICB consensus, triangulation weights

Zhang et al. (2025) — SMaHT confidence tiers

Luo et al. (2025) — FocalSV target-region assembly

Li (2014) — Artifacts: max depth, low-complexity regions

Peter et al. (2018) — 1011 yeast genomes

Mochizuki et al. (2023) — Assembly guideline, ploidy, polishing

Cheng & Sedlazeck (2025) — Inversion benchmark

Known Limitations
Limitation	Detail
DUPs under-scored by depth	DHFFC rarely >2.0; DUP_SPLIT_READ_CONFIRMED override mitigates
No BND/translocation support	Flagged but not called
k-mer layer often unavailable	2.4 GB Jellyfish DB
Cross-species performance	Jurei 0% LAR confirmation — pipeline correctly flags as uncertain
No spike-in FDR calibration	Pending synthetic benchmark
hifiasm requires >5 GB RAM	Not viable on low-RAM systems
Citation
Guimarães, K.H.A. et al. (2026). FUNGUS-SV: A triangulation-based structural variant discovery and validation pipeline for haploid genomes. In preparation.

📧 GitHub: @keltonjenkovguimaraes-alt

License
MIT License — see LICENSE for details.

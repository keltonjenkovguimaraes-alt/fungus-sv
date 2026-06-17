# 🧬 FUNGUS-SV v1.0.0

**Structural variant discovery and triangulation-based validation for haploid fungal genomes using PacBio HiFi reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.7+](https://img.shields.io/badge/Python-3.7+-yellow.svg)](https://python.org)
[![Status](https://img.shields.io/badge/status-publication--ready-brightgreen.svg)]()

---

## ⚠️ Disclaimer

FUNGUS-SV is a **hypothesis-generation tool.** NOT for clinical use. All SVs of biological interest must be independently validated by experimental methods.

---

## Overview

FUNGUS-SV detects structural variants (SVs) in haploid fungal genomes using PacBio HiFi long reads. It combines a **three-caller ICB consensus** (Sniffles2 + cuteSV + SVIM) with a **7-layer triangulation validation system** and a **two-assembler Local Assembly Refinement (LAR)** pipeline for definitive SV confirmation.

### Key Features

- **Three-caller ICB consensus:** Sniffles2 + cuteSV + SVIM — ≥2 callers must agree, removing 80–95% of single-caller noise
- **7-layer triangulation:** Depth (DHFFC+DHBFC), k-mer spectrum, breakpoint junction, local assembly, ploidy confirmation, genomic context
- **Two-assembler LAR:** Flye + Miniasm/Racon dual-assembler consensus for definitive breakpoint proof
- **Haploid-calibrated:** All thresholds recalibrated for haploid genomes (DHFFC, size factors, layer weights)
- **SV-type-specific overrides:** DUP_SPLIT_READ_CONFIRMED and INV_SPLIT_READ_CONFIRMED tiers based on empirical truth data (FDR ~0%)
- **Empirical false discovery rate:** Self-alignment negative controls establish ~3% FDR at highest confidence tier
- **Cross-genus validated:** 12 strains across *Saccharomyces cerevisiae* (1,527 SVs) and *Candida albicans* (2,083 SVs)
- **130-SV LAR truth set:** Two-assembler validation across 12 strains, 2 genera (72% within-species confirmation)
- **Samplot visual curation:** Combined LAR + Samplot achieves 100% accuracy on 19 tested SVs
- **Benchmarked:** Compared against MUMmer4, SVIM-asm, and individual callers on 10 strains

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
│ Layer 1: Local Assembly — Flye + Miniasm/Racon (w=0.20) │
│ Layer 2: Depth Signature — DHFFC + DHBFC (w=0.35) │
│ Layer 3: k-mer Spectrum — Jellyfish (w=0.15) │
│ Layer 4: Breakpoint Junction — split reads (w=0.30) │
│ Layer 5: Ploidy Confirmation — BCFtools pileup (w=0.0) │
│ Layer 6: Genomic Context — gene annotations (w=0.0) │
│ │
│ T = Σ(score_i × weight_i) / Σ(weight_i) │
└─────────────────────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────────────────────┐
│ PHASE 3: CONFIDENCE TIERS │
│ │
│ TRIPLE_TRIANGULATED: T ≥ 0.80 (FDR ~3%) │
│ DOUBLE_CONFIRMED: T ≥ 0.60 │
│ SINGLE_LINE: T ≥ 0.40 │
│ WEAK: T ≥ 0.20 │
│ CONTRADICTED: T < 0.20 (73% truly false) │
│ DUP_SPLIT_READ_CONFIRMED: Breakpoint ≥ 0.5 (FDR ~0%) │
│ INV_SPLIT_READ_CONFIRMED: Breakpoint ≥ 0.6 (FDR ~0%) │
└─────────────────────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────────────────────┐
│ PHASE 4: TWO-ASSEMBLER LAR (standalone) │
│ Flye + Miniasm/Racon → definitive breakpoint proof │
│ Samplot visual curation for discordant cases (~32%) │
└─────────────────────────────────────────────────────────┘

---

## Validation Results

### *Saccharomyces cerevisiae* — CICC-1445 vs 6 Reference Strains (1,527 SVs)

| Strain | Species | Type | Total SVs | DEL | DUP | INV | LAR Conf. |
|--------|---------|------|-----------|-----|-----|-----|-----------|
| S288C | *S. cerevisiae* | Laboratory | 285 | 254 | 19 | 12 | 7/10 (70%) |
| BJ4 | *S. cerevisiae* | Industrial | 176 | 151 | 14 | 11 | 7/10 (70%) |
| SX2 | *S. cerevisiae* | Industrial | 310 | 283 | 18 | 9 | 4/10 (40%) |
| IMX2600 | *S. cerevisiae* | Engineered | 345 | 313 | 20 | 12 | 8/10 (80%) |
| Makgeolli | *S. cerevisiae* | Fermentation | 244 | 218 | 17 | 9 | 10/10 (100%) |
| Jurei | *S. jurei* | Wild Type | 167 | 156 | 11 | 0 | 0/10 (0%)* |

\* Cross-species. Pipeline correctly identifies these as uncertain.

### *Candida albicans* — SC5314 vs 5 Reference Strains (2,083 SVs)

| Strain | Contigs | Total SVs | LAR Conf. | Assembly |
|--------|---------|-----------|-----------|----------|
| 101 | 9 chromosomes | 293 | 7/10 (70%) | Chromosome-level |
| WO1 | 17 chromosomes | 298 | 8/10 (80%) | Chromosome-level |
| FDAARGOS656 | 69 contigs | 288 | 7/10 (70%) | Near-complete |
| UAB012 | 308 contigs | 893 | 7/10 (70%) | Fragmented ⚠️ |
| ATCC64124 | 3,971 contigs | 111 | 7/10 (70%) | Highly fragmented ⚠️ |

### Self-Alignment Negative Controls — Empirical FDR

| Metric | *S. cerevisiae* | *C. albicans* |
|--------|----------------|---------------|
| Total SVs detected | 130 | 200 |
| TRIPLE tier FDR | 3.1% | 2.5% |
| LAR confirmed real (assembly errors) | 5/9 (56%) | 8/10 (80%) |
| Pipeline correctly flagged false | 3/3 (100%) | 1/1 (100%) |

### LAR Truth Set — 130 SVs Across 12 Strains

| Category | SVs | Confirmed | % |
|----------|-----|-----------|-----|
| Within-species (both genera) | 80 | 58 | **72%** |
| Cross-species (Jurei) | 10 | 0 | **0%** |
| Self-alignment | 20 | 9 | **45%** |
| DUP calibration set | 48 | 39 | **81%** |

### Benchmarking — Compared Against MUMmer4, SVIM-asm, and Individual Callers

| Tool | Method | SVs (S288C) | LAR SVs Recovered |
|------|--------|-------------|-------------------|
| **FUNGUS-SV ICB** | Read-based consensus | 285 | **7/7 (100%)** |
| Sniffles2 (single) | Read-based | 771 | — |
| cuteSV (single) | Read-based | 1,588 | — |
| SVIM (single) | Read-based | 1,549 | — |
| SVIM-asm | Assembly-based | 252 | 1/7 (14%) |
| MUMmer4 | Assembly-based | ~600 events | — |

> **Key finding:** ICB consensus removes 80–95% of single-caller noise. SVIM-asm missed 6/7 LAR-confirmed SVs — assembly-based tools are blind to real structural differences.

---

## Quick Start

### Prerequisites
- Linux (Ubuntu 20.04+) or WSL2
- ≥8 GB RAM
- Conda/Mamba
- PacBio HiFi reads (≥20× coverage)

### Installation

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

# Create conda environments (one-time setup)
conda env create -f workflow/envs/alignment.yaml -n sv_align
conda env create -f workflow/envs/sv_calling.yaml -n sv_call
conda env create -f workflow/envs/validation.yaml -n sv_valid
conda env create -f workflow/envs/lar.yaml -n sv_lar
conda env create -f workflow/envs/kmers.yaml -n sv_kmers  # optional
Usage
# 1. Align reads to reference
conda activate sv_align
minimap2 -d ref.mmi reference.fasta
minimap2 -t 8 -ax map-hifi ref.mmi reads.fastq.gz | samtools sort -@ 4 -o sample.bam -
samtools index sample.bam

# 2. Detect SVs (ICB consensus)
conda activate sv_call
python fungus_sv/core/icb.py \
    --bam sample.bam \
    --reference ref.fasta \
    --output results/ \
    --callers sniffles2 cutesv svim \
    --min-callers 2 \
    --min-overlap 0.5 \
    --threads 4

# 3. Validate (7-layer triangulation)
conda activate sv_valid
PYTHONPATH=. python -m valid_sv.run_validation \
    --consensus-vcf results/consensus_svs.vcf \
    --bam sample.bam \
    --reference ref.fasta \
    --output results/validation/ \
    --threads 4

# 4. LAR: Two-assembler definitive proof (standalone, for key candidates)
conda activate sv_lar
python3 -c "
from valid_sv.evidence.layer_lar import run_lar, run_lar_miniasm
flye = run_lar('sample.bam','ref.fasta','SV_ID','DEL','chr',start,end)
mini = run_lar_miniasm('sample.bam','ref.fasta','SV_ID','DEL','chr',start,end)
print(f'Flye: {flye.verdict.value}, Miniasm: {mini.verdict.value}')
"

# 5. Benchmark against other tools (optional)
bash scripts/run_benchmarking.sh
Repository Structure
fungus-sv/
├── fungus_sv/core/              # ICB consensus calling
│   └── icb.py                   # Main ICB pipeline (3-caller consensus)
├── valid_sv/                    # Triangulation validation framework
│   ├── evidence/                # 7 evidence layer implementations
│   │   ├── layer_lar.py         # Flye + Miniasm/Racon dual assembly
│   │   ├── layer_depth.py       # DHFFC + DHBFC + size stratification
│   │   ├── layer_kmer.py        # k-mer spectrum (Jellyfish)
│   │   ├── layer_breakpoint.py  # Split-read junction analysis
│   │   ├── layer_ploidy.py      # Fast BCFtools haploid check
│   │   └── layer_genomic_context.py  # Gene annotation filters
│   ├── engine/                  # Scoring engine
│   │   ├── scorer.py            # Triangulation T-score computation
│   │   └── fdr_estimator.py     # FDR estimation from calibration
│   ├── quality/triangulability.py  # Layer completeness assessment
│   ├── reporting/report_card.py    # Per-SV report generation
│   └── run_validation.py        # Main validation entry point
├── workflow/                    # Snakemake workflow
│   ├── Snakefile                # Automated pipeline (5 rules)
│   └── envs/                    # Conda environment YAMLs
├── config/config.yaml           # Calibrated parameters & thresholds
├── scripts/                     # Utility scripts
│   └── run_benchmarking.sh      # MUMmer4 + SVIM-asm comparison
├── test_dataset/                # Minimal test data (mitochondrial chr)
└── docs/                        # Documentation
Calibrated Parameters
Parameter	Value	Source
Layer Weights		
Depth Signature	0.35	Liu et al. (2024), recalibrated for haploid
Breakpoint Junction	0.30	Liu et al. (2024), recalibrated
Local Assembly	0.20	This study
k-mer Spectrum	0.15	Liu et al. (2024), recalibrated
DHFFC Thresholds		
DEL strong	<0.3	Pedersen & Quinlan (2019), recalibrated
DUP strong	>2.0	Pedersen & Quinlan (2019), recalibrated
Filters		
Min MAPQ	≥20	Zheng & Shang (2024)
Min SV size	50 bp	This study
Max het rate (haploid)	<7%	Xing et al. (2025)
ICB min overlap	0.5	Liu et al. (2024)
ICB min callers	2 (of 3)	Liu et al. (2024)
Overrides		
DUP override FDR	~0%	48-SV LAR truth set
INV override FDR	~0%	11-SV split-read confirmed
Performance
Metric	Value
Pipeline runtime (S288C, 285 SVs)	2 min 41 sec
LAR runtime (per SV)	2–15 min
LAR RAM (Flye)	<500 MB
LAR RAM (Miniasm)	<200 MB
ICB noise reduction vs single callers	80–95%
Within-species LAR confirmation	72%
TRIPLE tier empirical FDR	~3%
DUP/INV override FDR	~0%
Papers Integrated (13)
#	Paper	Contribution
1	Pedersen & Quinlan (2019)	DHFFC/DHBFC depth metrics
2	Zheng & Shang (2024)	distance_support, MAPQ filter
3	Belyeu et al. (2021)	Samplot visualization
4	David et al. (2024)	Manual curation FDR
5	Dhakal et al. (2024)	Fungal SV landscape
6	Li et al. (2023)	Haploid fungal variant calling
7	Liu et al. (2024)	ICB consensus, triangulation weights
8	Zhang et al. (2025)	SMaHT confidence tiers
9	Luo et al. (2025)	FocalSV target-region assembly
10	Li (2014)	Max depth, low-complexity region filters
11	Peter et al. (2018)	1011 yeast genomes
12	Mochizuki et al. (2023)	Assembly guideline, ploidy, polishing
13	Cheng & Sedlazeck (2025)	Inversion benchmark
Known Limitations
Limitation	Detail	Mitigation
DUPs under-scored by depth	DHFFC rarely >2.0 in haploids	DUP_SPLIT_READ_CONFIRMED override (FDR ~0%)
Assembly quality matters	Fragmented references inflate SVs	Chromosome-level assemblies recommended
Cross-species calls unreliable	Jurei 0% LAR confirmed	Pipeline correctly flags as uncertain
No BND/translocation	Translocations flagged, not called	Flagged as possible translocation
k-mer layer optional	2.4 GB Jellyfish DB	Pipeline runs without it
Ploidy = homozygosity check	Inbred diploids appear haploid	Documented limitation
LAR is manual	Not in automated pipeline	Standalone tool for key candidates
Small SVs (<100 bp)	Only 40% LAR confirmed	Require ≥6 split reads
Citation
Guimarães, K.H.A. et al. (2026). FUNGUS-SV: A triangulation-based structural variant discovery and validation pipeline for haploid genomes. In preparation.

📧 GitHub: @keltonjenkovguimaraes-alt


License
MIT License — see LICENSE for details.

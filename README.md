# 🧬 FUNGUS-SV v0.9.4

**Structural variant discovery and triangulation-based validation for haploid fungal genomes using PacBio HiFi reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)
[![Status](https://img.shields.io/badge/status-active%20development-orange.svg)]()

---

## ⚠️ Disclaimer

FUNGUS-SV is a **hypothesis-generation tool under active development.** NOT for clinical use. All SVs of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing).

---

## Overview

FUNGUS-SV detects structural variants (SVs) in haploid fungal genomes using PacBio HiFi long reads. It combines a **three-caller ICB consensus** (Sniffles2 + cuteSV + SVIM) with a **multi-layer validation system** that scores each SV using orthogonal evidence and SV-type-specific filters.

The pipeline was calibrated and validated on *Saccharomyces cerevisiae* CICC-1445 against five reference strains (S288C, BJ4, IMX2600, Makgeolli, SX2).

### Pipeline Architecture

PacBio HiFi reads (*.fastq.gz)
│
▼
[minimap2] ← map-hifi preset
│
▼
┌─────────────────────────────────────────┐
│ PHASE 1: ICB CONSENSUS │
│ Sniffles2 │ cuteSV │ SVIM │
│ Consensus: ≥2 callers, 0.5 overlap │
└─────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────┐
│ PHASE 2: TRIANGULATION VALIDATION │
│ │
│ Layer 1: Local Assembly (Flye/Miniasm) │
│ Layer 2: Depth (DHFFC + DHBFC) │
│ Layer 3: k-mer Spectrum (Jellyfish) │
│ Layer 4: Breakpoint Junction │
│ Layer 5: Ploidy Confirmation (HARD FILTER)│
│ Layer 6: Genomic Context (HARD FILTER) │
│ │
│ T = Σ(score × weight) / Σ(weights) │
└─────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────┐
│ PHASE 3: CONFIDENCE TIERS │
│ │
│ TRIPLE_TRIANGULATED: T ≥ 0.80 │
│ DOUBLE_CONFIRMED: T ≥ 0.60 │
│ SINGLE_LINE: T ≥ 0.40 │
│ WEAK: T ≥ 0.20 │
│ CONTRADICTED: T < 0.20 │
│ INV_SPLIT_READ_CONFIRMED │
└─────────────────────────────────────────┘

### Calibrated Weights (Haploid Fungi)

| Layer | Weight | Rationale |
|-------|--------|-----------|
| Depth (DHFFC + DHBFC) | **0.35** | 100% depth loss unambiguous in haploids |
| Breakpoint Junction | **0.30** | Only signal for inversions; critical for all types |
| Local Assembly (LAR) | **0.20** | Computationally expensive; standalone tool available |
| k-mer Spectrum | **0.15** | Redundant with depth in haploids |
| Ploidy Confirmation | **0.00** | Hard filter, not scored |
| Genomic Context | **0.00** | Hard filter, not scored |

---

## Key Features

### 1. Multi-Assembler LAR (Local Assembly Refinement)
Definitive proof of SV existence using regional assembly:
- **Tier 1:** Flye (<500 MB RAM, 2-15 min) — default
- **Tier 2:** Miniasm + Racon (<200 MB RAM, <2 min) — orthogonal confirmation
- **Tier 3:** hifiasm — high-accuracy for unresolved calls (requires >5 GB RAM)

### 2. Haploid-Calibrated Depth Analysis
- DHFFC + DHBFC combined ratio
- Size-stratified scoring factors
- Max depth filter (Li 2014)
- Low-complexity region flag

### 3. Genomic Context Filters (Layer 6)
SV-type-specific hard filters using gene annotation:
- S288C-specific insertion detection (Ty2 elements)
- Essential gene dosage warnings
- Mitochondrial genome exception
- Size-stratified reliability for small DELs
- Repeat-flanked inversion escalation

### 4. Multi-Strain Validation
Tested against 5 reference strains with LAR-validated truth sets:
- S288C (laboratory) — best reference quality, 53% LAR confirmation
- BJ4, SX2 (industrial) — Chinese strain relatives
- IMX2600 (engineered) — assembly differences cause false signals
- Makgeolli (fermentation)

---

## Validation Results

### CICC-1445 vs S288C (277 SVs)

| Tier | Count | % |
|------|-------|---|
| TRIPLE_TRIANGULATED | 5 | 1.8% |
| DOUBLE_CONFIRMED | 17 | 6.1% |
| SINGLE_LINE | 1 | 0.4% |
| WEAK | 92 | 33.2% |
| CONTRADICTED | 151 | 54.5% |
| INV_SPLIT_READ_CONFIRMED | 11 | 4.0% |

### LAR-Calibrated Truth Data

| Strain | SV Type | Tested | LAR Real | % Real |
|--------|---------|--------|----------|--------|
| S288C | DEL (5-6.5 kb) | 15 | 8 | 53% |
| S288C | DUP (largest 5) | 5 | 4 | 80% |
| SX2 | INV | 3 | 2 | 67% |
| IMX2600 | DEL | 28 | 5 | 18% |
| Makgeolli | DEL | 16 | 3 | 19% |

---

## Quick Start

### Prerequisites
- Linux (Ubuntu 20.04+) | ≥16 GB RAM | Conda/Mamba
- PacBio HiFi reads (≥20× coverage)

### Installation

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

# Create environments
conda create -n sv_align -c bioconda -c conda-forge minimap2 samtools -y
conda create -n sv_call -c bioconda -c conda-forge sniffles=2.2 cutesv svim bcftools -y
conda create -n sv_valid -c conda-forge -c bioconda python=3.11 numpy scipy pandas pysam pyyaml matplotlib -y
conda create -n sv_lar -c bioconda flye minimap2 samtools miniasm racon -y
conda create -n sv_samplot -c bioconda samplot -y
Usage
# 1. Align reads to reference
conda activate sv_align
minimap2 -t 8 -ax map-hifi -R '@RG\tID:sample\tSM:sample' ref.fasta reads.fastq.gz \
    | samtools sort -@ 4 -o sample.sorted.bam -
samtools index sample.sorted.bam

# 2. Detect SVs (ICB consensus)
conda activate sv_call
python fungus_sv/core/icb.py --bam sample.sorted.bam --reference ref.fasta \
    --output results/ --threads 4

# 3. Validate (triangulation + genomic context)
conda activate sv_valid
PYTHONPATH=. python -m valid_sv.run_validation \
    --consensus-vcf results/consensus_svs.vcf \
    --bam sample.sorted.bam --reference ref.fasta --fastq reads.fastq.gz \
    --output results/validation/ --threads 4

# 4. LAR: Definitive proof for key SVs (standalone)
conda activate sv_lar
python valid_sv/evidence/layer_lar.py
Repository Structure
fungus-sv/
├── config/                      # Pipeline configuration
├── fungus_sv/core/              # ICB consensus calling
├── valid_sv/                    # Triangulation validation
│   ├── evidence/                # Layer implementations
│   │   ├── layer_lar.py         # Flye + Miniasm+Racon assembly
│   │   ├── layer_depth.py       # DHFFC + DHBFC + max depth + LCR
│   │   ├── layer_kmer.py        # k-mer spectrum
│   │   ├── layer_breakpoint.py  # Split-read junction analysis
│   │   ├── layer_ploidy.py      # Haploid confirmation
│   │   └── layer_genomic_context.py  # Gene annotation filters
│   ├── engine/                  # Scoring engine
│   └── reporting/               # Report generation
├── data/yeast/                  # Validation data (CICC-1445 vs 5 refs)
├── docs/                        # Documentation
│   ├── DEVELOPMENT_LOGBOOK.md   # Chronological development log
│   └── IMPLEMENTATION_SUMMARY_v0.9.4.md  # Implementation details
└── figures/                     # Samplot images and plots
Known Limitations
Limitation	Detail
DUPs under-scored	Depth layer cannot confirm DUPs (DHFFC rarely >2.0); LAR required
Inversions	Breakpoint-only scoring; depth/k-mer silent for balanced events
Small SVs (<100 bp)	Unreliable depth/breakpoint signals
Repetitive regions	FLO, rDNA, Ty elements produce complex signals
k-mer layer	Jellyfish DB is 2.4 GB; often unavailable
No spike-in calibration	FDR estimates approximate; truth set from LAR
No BND/translocation support	Flagged but not called
hifiasm requires >5 GB RAM	Not viable on low-RAM systems for Tier 3 assembly
Papers Integrated
Pedersen & Quinlan (2019) — DHFFC/DHBFC depth metrics

Zheng & Shang (2024) — distance_support, MAPQ filter

Belyeu et al. (2021) — Samplot visualization, ML curation

David et al. (2024) — Manual curation FDR estimates

Dhakal et al. (2024) — Fungal SV landscape

Li et al. (2023) — Haploid fungal variant calling

Liu et al. (2024) — ICB consensus, triangulation weights

Zhang et al. (2025) — SMaHT benchmark, confidence tiers

Luo et al. (2025) — FocalSV: target region assembly SV detection

Li (2014) — Artifacts: max depth filter, low-complexity regions

Peter et al. (2018) — 1011 yeast genomes

Mochizuki et al. (2023) — Assembly guideline, ploidy, polishing

Cheng & Sedlazeck (2025) — Inversion benchmark

Citation
Guimarães, K.H.A. et al. (2026). FUNGUS-SV: A triangulation-based structural variant discovery and validation pipeline for haploid genomes. In preparation.

📧 GitHub: @keltonjenkovguimaraes-alt

License
MIT License — see LICENSE for details.

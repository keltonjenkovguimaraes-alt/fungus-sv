# 🧬 FUNGUS-SV

**A benchmark-free structural variant discovery and validation pipeline for non-model haploid fungi using PacBio HiFi long reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)
[![Snakemake](https://img.shields.io/badge/Snakemake-8.0-blue.svg)](https://snakemake.github.io)

---

## 📖 Overview

Structural variant (SV) calling in non-model organisms faces a fundamental challenge: **how do you validate variants when no benchmark truth set exists?** For humans, the Genome in a Bottle (GIAB) consortium provides gold-standard calls. For *Sporothrix schenckii*, *Candida auris*, *Aspergillus fumigatus*, and countless other medically important fungi — nothing.

FUNGUS-SV solves this through **three methodological innovations**:

| Innovation | Description |
|-----------|-------------|
| 🔵 **ICB** (Iterative Consensus Builder) | Runs 3 orthogonal SV callers (pbsv, Sniffles2, cuteSV), clusters overlapping SVs by reciprocal overlap, scores by multi-caller agreement, and iteratively refines using high-confidence calls as pseudo-truth |
| 🟢 **LAR** (Local Assembly Refinement) | Extracts reads spanning SV breakpoints, performs local de novo assembly with Flye, realigns to reference, and corrects breakpoint coordinates to base-pair precision |
| 🟣 **Pan-Atlas** | Maps SVs across species via ortholog matching, classifies variants as conserved or species-specific, annotates functional impact on virulence-associated genes |

---

## 🎯 Why Haploid Fungi?

FUNGUS-SV is explicitly designed for **haploid organisms**. This is not a limitation — it is a deliberate design choice that simplifies variant calling (no heterozygosity, no allelic phasing) and makes the pipeline faster and more interpretable. Compatible organisms include:

- All haploid fungi (*Sporothrix*, *Candida*, *Aspergillus*, *Histoplasma*, *Coccidioides*, *Cryptococcus*, *Fusarium*)
- Bacteria and archaea
- Haploid plants (moss, algal gametophytes)
- Any haploid eukaryote with PacBio HiFi data at ≥30× coverage

---

## 📊 Pipeline Architecture

```
PacBio HiFi reads (*.fastq.gz)
        │
        ▼
   [minimap2]  ← map-hifi preset, 99.9% mapping rate
        │
        ▼
   ┌────────────────────────────┐
   │  ICB: 3-Way SV Calling     │
   │  pbsv + Sniffles2 + cuteSV │
   │         ↓                  │
   │  Consensus Scoring         │
   │         ↓                  │
   │  Iterative Refinement      │
   └────────────────────────────┘
        │
        ▼
   ┌────────────────────────────┐
   │  LAR: Local Assembly       │
   │  samtools → Flye → minimap2│
   │         ↓                  │
   │  Breakpoint Correction     │
   └────────────────────────────┘
        │
        ▼
   ┌────────────────────────────┐
   │  Pan-Atlas: Cross-Species  │
   │  Ortholog Matching         │
   │         ↓                  │
   │  Conservation Analysis     │
   └────────────────────────────┘
```

---

## 🧬 Validation: *Sporothrix schenckii* NBRC32961

The pipeline was developed and validated on a medically important dimorphic fungal pathogen:

### Sequencing

| Metric | Value |
|--------|-------|
| Platform | PacBio Revio |
| Chemistry | HiFi (CCS) |
| Total reads | 197,830 |
| Read N50 | ~18 kb |
| Mean coverage | 58× |
| SRA accession | DRR631664 |

### Reference Genome

| Metric | Value |
|--------|-------|
| Species | *Sporothrix schenckii* 1099-18 |
| Assembly | GCF_000961545.1 |
| Size | 32.4 Mb |
| Contigs | 16 |
| N50 | 4.3 Mb |
| Protein-coding genes | 10,389 |

### Results

| Metric | Value |
|--------|-------|
| Mapping rate | **99.91%** |
| Raw SV calls (3 callers) | 1,915 |
| ICB consensus SVs | **275** |
| High-confidence (3-caller) | **174** |
| SV size range | 29 bp – 139 kb |
| Genes affected (HIGH impact) | **83** |
| Conserved in *S. brasiliensis* | 72 (86.7%) |
| Species-specific genes | 11 (13.3%) |

### Key Virulence Genes Affected by SVs

| Gene | Product | Relevance |
|------|---------|-----------|
| SPSK_02606 | Superoxide dismutase (Fe-Mn) | Oxidative stress resistance |
| SPSK_04255 | GPI-anchored protein | Cell wall — host adhesion |
| SPSK_02728 | Secretory component Shr3 | ER protein secretion |
| SPSK_02129 | COPII vesicle coat protein | Virulence factor delivery |
| SPSK_02379 | WASP-interacting protein | Actin cytoskeleton |
| SPSK_02154 | Phospholipid-translocating ATPase | **S. schenckii-specific** |

### LAR Validation: 139 kb → 652 bp Inversion

| Metric | ICB (pre-LAR) | LAR (post-refinement) |
|--------|--------------|----------------------|
| Inversion size | 139,348 bp | **652 bp** |
| Assembly identity | N/A | **99.85%** |
| Reads used | N/A | 850 (70× local coverage) |

This demonstrates why LAR is essential: alignment-based callers systematically overestimate large SV sizes.

---

## ⚡ Quick Start

### Prerequisites

- Linux (Ubuntu 20.04+ recommended)
- ≥32 GB RAM, ≥8 CPU cores
- Conda or Mamba
- PacBio HiFi reads at ≥30× coverage

### Installation

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv
conda env create -f workflow/envs/environment.yaml
conda activate snp_svant_pacbio
```

### Prepare Your Data

```
data/raw/your_sample.fastq.gz    # PacBio HiFi reads
data/reference/reference.fasta    # Reference genome
data/reference/reference.gff      # Gene annotations
```

### Configure

Edit `config/config.yaml` with your sample name and file paths.

### Run

```bash
snakemake -s workflow/Snakefile --cores 8
```

---

## 📁 Repository Structure

```
fungus-sv/
├── README.md                      # This file
├── LICENSE                        # MIT License
├── workflow/
│   ├── Snakefile                  # Main Snakemake pipeline
│   └── envs/environment.yaml      # Conda environment
├── fungus_sv/
│   ├── core/
│   │   ├── icb.py                 # Iterative Consensus Builder
│   │   ├── build_consensus.py     # Consensus scoring algorithm
│   │   └── annotate_svs.py        # Custom SV annotator
│   └── modules/
│       └── local_assembly.py      # LAR module
├── config/config.yaml             # Pipeline configuration
├── docs/methods.md                # Detailed methods
└── paper/                         # Manuscript directory
```

---

## 🔬 Methods Summary

### ICB: Iterative Consensus Builder

- **Step 1**: Run pbsv, Sniffles2, and cuteSV on aligned BAM
- **Step 2**: Parse all VCFs, cluster SVs by reciprocal overlap (≥0.5)
- **Step 3**: Score clusters by caller agreement (1-3 callers)
- **Step 4**: Output consensus VCF with SUPPORT and CALLERS tags
- **Step 5**: Iteratively upgrade 2-caller SVs matching truth profile

### LAR: Local Assembly Refinement

- **Step 1**: Extract reads spanning SV ±5 kb (samtools)
- **Step 2**: Local de novo assembly (Flye, `--pacbio-hifi`)
- **Step 3**: Align assembly to reference (minimap2, `-cx asm5`)
- **Step 4**: Extract refined breakpoints from PAF

### Pan-Atlas

- Ortholog matching via locus tag conversion (SPSK_ → SPBR_)
- Classification: Conserved (ortholog present) vs. Species-specific
- Functional annotation from GFF mRNA product fields

---

## 📚 References

This pipeline builds upon:

1. **Gunasekaran et al. (2024)** — SNP-SVant: computational workflow for organisms lacking benchmarked variants. *Current Protocols*, 4, e1046.
2. **Holt et al. (2024)** — HiPhase: jointly phasing small, structural, and tandem repeat variants from HiFi sequencing. *Bioinformatics*, 40, btae042.
3. **Liu et al. (2024)** — Tradeoffs in alignment and assembly-based methods for structural variant detection with long-read sequencing data. *Nature Communications*, 15, 2447.
4. **Hops et al. (2025)** — HiFi long-read genomes for difficult-to-detect, clinically relevant variants. *American Journal of Human Genetics*, 112, 450-456.
5. **Teixeira et al. (2014)** — Comparative genomics of *Sporothrix schenckii* and *Sporothrix brasiliensis*. *BMC Genomics*, 15, 943.
6. **Hartmann et al. (2022)** — vembrane: filtering and transforming VCF/BCF files. *Bioinformatics*, 38, 5300-5302.

---

## 📄 Citation

> Guimarães, K.H.A; Philippsen H.K., et al. (2026). FUNGUS-SV: A benchmark-free structural variant discovery pipeline for non-model haploid fungi using PacBio HiFi sequencing. *In preparation*.

---

## 📧 Contact

**Kelton Jenkov Guimarães**
GitHub: [@keltonjenkovguimaraes-alt](https://github.com/keltonjenkovguimaraes-alt)

---

## 🤝 Acknowledgments

Built for the *Sporothrix* research community and all scientists working on non-model fungal pathogens.

*If this pipeline helps your research, please star ⭐ the repository and cite the paper.*

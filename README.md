# 🧬 FUNGUS-SV

**A structural variant discovery and triangulation-based prioritization pipeline for non-model haploid fungi using PacBio HiFi long reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)
[![Snakemake](https://img.shields.io/badge/Snakemake-8.0-blue.svg)](https://snakemake.github.io)

---

## ⚠️ CRITICAL: Read Before Using

**FUNGUS-SV is a hypothesis-generation tool, not a truth machine.**

There is no benchmark truth set for *Sporothrix schenckii* (or most non-model fungi). No external gold standard exists to validate structural variant calls against.

### What This Pipeline Actually Provides

| Pipeline Output | Honest Interpretation |
|----------------|----------------------|
| "275 ICB consensus SVs" | 275 SVs where ≥2 of 3 alignment-based callers agree |
| "174 three-caller SVs" | 174 SVs where all 3 callers agree on presence of an SV |
| "50 triple-triangulated SVs" | 50 SVs supported by ≥2 orthogonal evidence layers (not including caller agreement) |
| "T-score" | A weighted evidence score (0-1), NOT a probability. Higher = more orthogonal support |
| "Confidence estimate" | Approximate, based on mixture modeling of T-score distribution. **NOT calibrated** |

### Known Limitations (Listed First For Transparency)

1. **LAR is NOT run automatically.** The local assembly layer is a placeholder. You must run `fungus_sv/modules/local_assembly.py` separately. Currently all SVs show LAR as "not_run".

2. **Caller agreement is NOT used in T-score.** We report ICB support but exclude it from triangulation scoring to avoid circular validation (using the prediction as its own validation).

3. **No external truth set exists.** T-scores and confidence estimates are internally derived and APPROXIMATE.

4. **Systematic false positives are possible.** All SV callers share alignment assumptions. A reference assembly error could produce confident but false consensus calls.

5. **Small SVs (<100 bp) have limited orthogonal validation.** Depth signature requires ≥100 bp. k-mer resolution depends on k-mer size (31 bp).

6. **False negative rate is unknown.** SVs missed by all three callers are invisible.

7. **These results have not been peer-reviewed.** Use with appropriate skepticism.

**Any SV of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing) before functional interpretation.**

---

## 📖 Overview

Structural variant (SV) calling in non-model organisms faces a fundamental challenge: **how do you assess confidence when no benchmark truth set exists?**

FUNGUS-SV addresses this through a **two-phase architecture**:

| Phase | Component | Description |
|-------|-----------|-------------|
| 🔵 **Prediction** | ICB (Iterative Consensus Builder) | Runs 3 SV callers (pbsv, Sniffles2, cuteSV), clusters overlapping SVs, scores by multi-caller agreement |
| 🟢 **Validation** | VALID-SV Triangulation Core | Combines 4-5 independent evidence layers with uncorrelated failure modes. **Caller agreement is reported but excluded from T-score** to avoid circular validation |
| 🟣 **Annotation** | Pan-Atlas | Maps SVs across species via ortholog matching and annotates functional impact |

---

## 🎯 Why Haploid Fungi?

FUNGUS-SV is explicitly designed for **haploid organisms**. This simplifies variant calling (no heterozygosity, no allelic phasing). Compatible organisms include:

- All haploid fungi (*Sporothrix*, *Candida*, *Aspergillus*, *Histoplasma*, *Coccidioides*, *Cryptococcus*, *Fusarium*)
- Bacteria and archaea
- Haploid plants (moss, algal gametophytes)
- Any haploid eukaryote with PacBio HiFi data at ≥30× coverage

---

## 📊 Pipeline Architecture

PacBio HiFi reads (*.fastq.gz)
│
▼
[minimap2] ← map-hifi preset
│
▼
┌────────────────────────────────────────┐
│ PHASE 1: PREDICTION (ICB) │
│ pbsv + Sniffles2 + cuteSV │
│ ↓ │
│ Consensus clustering + scoring │
│ Output: 275 candidate SVs │
└────────────────────────────────────────┘
│
▼
┌────────────────────────────────────────┐
│ PHASE 2: VALIDATION (VALID-SV) │
│ │
│ Layer 1: ICB agreement (REPORTED ONLY)│
│ Layer 2: Local Assembly (MANUAL) │
│ Layer 3: Read-Depth Signature │
│ Layer 4: k-mer Spectrum │
│ Layer 5: Breakpoint Junctions │
│ Layer 6: Ploidy Confirmation (SNV) │
│ ↓ │
│ Triangulation Engine → T-score │
└────────────────────────────────────────┘
│
▼
┌────────────────────────────────────────┐
│ PHASE 3: ANNOTATION + REPORTING │
│ Gene overlap + Pan-Atlas + Report cards│
└────────────────────────────────────────┘

---

## 🧬 Results: *Sporothrix schenckii* NBRC32961

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

> **Note:** NBRC32961 reads are mapped to strain 1099-18 reference. SVs may include strain-specific structural differences.

### ICB Prediction Results

| Metric | Value |
|--------|-------|
| Mapping rate | 99.91% |
| Raw SV calls (3 callers) | 1,915 |
| ICB consensus SVs (≥2 callers) | **275** |
| Three-caller agreement | **174** |
| Two-caller agreement | **101** |
| SV size range | 29 bp – 139 kb |
| SV types | 269 DEL, 5 INV, 1 DUP |

### VALID-SV Triangulation Results

**These results are from computational triangulation only. They are not experimentally validated.**

| Triangulation Tier | T-score | 3-caller SVs (n=174) | 2-caller SVs (n=101) | All SVs (n=275) |
|-------------------|---------|---------------------|---------------------|-----------------|
| Triple-triangulated | ≥0.80 | 45 (25.9%) | 5 (5.0%) | 50 (18.2%) |
| Double-confirmed | 0.60–0.79 | 51 (29.3%) | 5 (5.0%) | 56 (20.4%) |
| Single-line evidence | 0.40–0.59 | 32 (18.4%) | 58 (57.4%) | 90 (32.7%) |
| Weak/contradicted | <0.40 | 46 (26.4%) | 33 (32.7%) | 79 (28.7%) |

**Key finding:** Of the 174 SVs where all three callers agree, 46 (26.4%) show weak or contradictory evidence from orthogonal layers. This demonstrates why consensus-alone approaches are insufficient and orthogonal validation is essential.

### Evidence Layer Performance

| Layer | Independence | Used in T-score | Works best for | Limitation |
|-------|-------------|-----------------|----------------|------------|
| ICB consensus | Circular — excluded | **No** (reported only) | All SV types | Shares alignment assumptions |
| Local assembly (LAR) | Medium (different paradigm) | **Yes** — but must run manually | ≥500 bp | Currently placeholder |
| Read-depth signature | High (counts reads) | **Yes** | DEL, DUP ≥100 bp | Not for INV |
| k-mer spectrum | Highest (no reference) | **Yes** | DEL, INS | Requires jellyfish DB |
| Breakpoint junctions | Medium (split-read) | **Yes** | All types | Low sensitivity on HiFi |
| Ploidy confirmation | Independent (SNV het rate) | **Yes** | All (validates haploid assumption) | Requires Longshot |

---

## ⚡ Quick Start

### Prerequisites

- Linux (Ubuntu 20.04+ recommended)
- ≥32 GB RAM, ≥8 CPU cores
- Conda or Mamba
- PacBio HiFi reads at ≥30× coverage

### Installation

```bash
# Clone the repository
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

# Create and activate conda environment
conda env create -f workflow/envs/environment.yaml
conda activate snp_svant_pacbio

# Verify installation
python --version  # Should show Python 3.11
snakemake --version  # Should show 8.0+
minimap2 --version
samtools --version
Prepare Data
# Create data directories
mkdir -p data/raw data/reference

# Copy your data files
cp /path/to/your_sample.fastq.gz data/raw/
cp /path/to/reference.fasta data/reference/
cp /path/to/reference.gff data/reference/  # Optional: for annotation

# Index reference
samtools faidx data/reference/reference.fasta
minimap2 -d data/reference/reference.mmi data/reference/reference.fasta
Run Pipeline
# Phase 1: SV Prediction
snakemake -s workflow/Snakefile --cores 8

# Phase 2: SV Validation
python -m valid_sv.run_validation \
    --consensus-vcf results/variants/consensus/sample.consensus_svs.vcf \
    --bam results/alignment/sample.sorted.bam \
    --reference data/reference/reference.fasta \
    --fastq data/raw/sample.fastq.gz \
    --output results/validation/
# LAR must be run manually for each SV of interest
python fungus_sv/modules/local_assembly.py \
    --consensus results/variants/consensus/sample.consensus_svs.vcf \
    --bam results/alignment/sample.sorted.bam \
    --reference data/reference/reference.fasta \
    --output results/variants/refined/sample.refined_svs.vcf.gz
📁 Repository Structure
fungus-sv/
├── README.md                         # This file
├── LICENSE                           # MIT License
├── workflow/
│   ├── Snakefile                     # Main Snakemake pipeline
│   └── envs/environment.yaml         # Conda environment
├── fungus_sv/                        # Prediction phase
│   ├── core/
│   │   ├── icb.py                    # Iterative Consensus Builder
│   │   ├── build_consensus.py        # Consensus scoring algorithm
│   │   └── annotate_svs.py           # Custom SV annotator
│   └── modules/
│       └── local_assembly.py         # LAR module (run separately)
├── valid_sv/                         # Validation phase
│   ├── evidence/
│   │   ├── layer_depth.py            # Read-depth signature
│   │   ├── layer_kmer.py             # k-mer spectrum analysis
│   │   ├── layer_breakpoint.py       # Breakpoint junction analysis
│   │   └── layer_ploidy.py           # Ploidy confirmation via SNV het rate
│   ├── quality/
│   │   └── triangulability.py        # Assess which layers can validate each SV
│   ├── engine/
│   │   ├── scorer.py                 # Weighted T-score computation
│   │   └── fdr_estimator.py          # Mixture model confidence estimation
│   ├── benchmarks/
│   │   └── spike_in.py               # Synthetic benchmark for calibration
│   ├── reporting/
│   │   └── report_card.py            # Per-SV evidence summaries
│   └── run_validation.py             # Main validation entry point
├── config/config.yaml                # Pipeline configuration
└── docs/methods.md                   # Detailed methods
 References
Ko & Brandizzi (2025) — NEEDLE: Network-enabled gene discovery pipeline for non-model plants. Cell Reports Methods, 5, 100963.

Turner et al. (2026) — OrthoGarden: Phylogenomics for non-model organisms without reference orthologs. Molecular Biology and Evolution, 43, msag053.

Salam et al. (2025) — FAIRification of DMRichR for non-model epigenetics. Bioinformatics Advances, vbaf024.

Gunasekaran et al. (2024) — SNP-SVant: computational workflow for organisms lacking benchmarked variants. Current Protocols, 4, e1046.

Holt et al. (2024) — HiPhase: jointly phasing small, structural, and tandem repeat variants. Bioinformatics, 40, btae042.

Liu et al. (2024) — Tradeoffs in alignment and assembly-based methods for SV detection. Nature Communications, 15, 2447.

Teixeira et al. (2014) — Comparative genomics of Sporothrix schenckii and Sporothrix brasiliensis. BMC Genomics, 15, 943.
📄 Citation
Guimarães, K.H.A; Philippsen H.K., et al. (2026). FUNGUS-SV: A structural variant discovery and triangulation-based prioritization pipeline for non-model haploid fungi using PacBio HiFi sequencing. In preparation.

📧 Contact
Kelton Jenkov Guimarães
GitHub: @keltonjenkovguimaraes-alt

🤝 Acknowledgments
Built for the Sporothrix research community and all scientists working on non-model fungal pathogens. The VALID-SV triangulation approach was inspired by NEEDLE's prediction-validation split, OrthoGarden's all-vs-all evidence inference, and DMRichR's FAIRification for non-model organisms.

If this pipeline helps your research, please star ⭐ the repository and cite the paper.

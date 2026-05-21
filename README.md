# 🧬 FUNGUS-SV

**A structural variant discovery and triangulation-based prioritization pipeline for haploid genomes using PacBio HiFi reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)

---

## ⚠️ HONEST ASSESSMENT

FUNGUS-SV is a **hypothesis-generation tool under active development.** It is NOT production-ready.

### What This Pipeline Does

Detects structural variants (SVs ≥50 bp) in haploid genomes by combining multiple SV callers through an Intersection-Consensus-Builder (ICB), then validating each consensus SV using five orthogonal evidence layers:

| Layer | Method | What It Measures |
|-------|--------|-----------------|
| Local Assembly Refinement (LAR) | Flye assembly of reads at breakpoint | Confirms SV by assembling the variant allele |
| Read-Depth Signature | Coverage drop/increase at SV region | Detects deletions/duplications via copy number |
| k-mer Spectrum | Jellyfish k-mer presence/absence | Confirms sequence gain/loss independent of alignment |
| Breakpoint Junction | Split-read and soft-clip analysis | Confirms precise breakpoint locations |
| Ploidy Confirmation | SNV heterozygosity rate | Verifies haploid assumption |

### First Real Results (Cross-Species Validation)

| Metric | Value |
|--------|-------|
| Reference | *Acinetobacter bouvetii* JCM 18991 (3.4 Mb) |
| Reads | *Acinetobacter baumannii* ATCC 19606 PacBio HiFi (19,568 reads, ~82×) |
| Sniffles2 calls | 6 |
| cuteSV calls | 9 |
| **ICB Consensus (≥2 callers)** | **5 deletions** |
| DOUBLE_CONFIRMED (T≥0.6) | 2 (297 bp, 205 bp) |
| WEAK (T<0.4) | 3 (60-76 bp) |

**Key finding:** Larger deletions (>200 bp) received stronger orthogonal validation than small deletions (<100 bp), consistent with expectations from the literature.

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
│ Sniffles2 + cuteSV + SVIM + pbsv │
│ ↓ │
│ Consensus clustering (≥2 callers) │
│ Output: candidate SVs │
└────────────────────────────────────────┘
│
▼
┌────────────────────────────────────────┐
│ PHASE 2: VALIDATION (VALID-SV) │
│ │
│ Layer 1: Local Assembly (Flye) │
│ Layer 2: Read-Depth Signature │
│ Layer 3: k-mer Spectrum (Jellyfish) │
│ Layer 4: Breakpoint Junctions │
│ Layer 5: Ploidy Confirmation │
│ ↓ │
│ Triangulation Engine → T-score (0-1) │
└────────────────────────────────────────┘
│
▼
┌────────────────────────────────────────┐
│ PHASE 3: REPORTING │
│ Per-SV report cards + size-stratified │
│ summary + SMaHT confidence tiers │
└────────────────────────────────────────┘

---

## 🚀 Quick Start

### Prerequisites

- Linux (Ubuntu 20.04+)
- ≥16 GB RAM (32 GB recommended)
- Conda or Mamba
- PacBio HiFi reads (≥20× coverage)

### Installation

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

# Create isolated conda environments
conda create -n sv_align -c bioconda -c conda-forge minimap2 samtools -y
conda create -n sv_call -c bioconda -c conda-forge sniffles=2.2 cutesv svim bcftools -y
conda create -n sv_valid -c conda-forge -c bioconda python=3.11 numpy scipy pandas pysam pyyaml samtools minimap2 -y
conda create -n sv_lar -c bioconda -c conda-forge python=3.10 flye samtools minimap2 -y
conda create -n sv_kmers -c bioconda -c conda-forge python=3.10 jellyfish -y
Usage
# 1. Align reads to reference
conda activate sv_align
minimap2 -t 8 -ax map-hifi -R '@RG\tID:sample\tSM:sample' ref.fasta reads.fastq.gz \
    | samtools sort -@ 4 -o sample.sorted.bam -
samtools index sample.sorted.bam
conda deactivate

# 2. Run ICB consensus (multi-caller SV detection)
conda activate sv_call
python fungus_sv/core/icb.py \
    --bam sample.sorted.bam \
    --reference ref.fasta \
    --output results/variants/ \
    --callers sniffles2 cutesv svim \
    --min-callers 2
conda deactivate

# 3. Run orthogonal validation
conda activate sv_valid
python -m valid_sv.run_validation \
    --consensus-vcf results/variants/consensus_svs.vcf \
    --bam sample.sorted.bam \
    --reference ref.fasta \
    --fastq reads.fastq.gz \
    --output results/validation/
conda deactivate
📁 Repository Structure
fungus-sv/
├── workflow/Snakefile                 # Snakemake workflow
├── workflow/envs/                     # Conda environment YAMLs
├── fungus_sv/
│   ├── core/icb.py                    # ICB consensus builder
│   └── modules/local_assembly.py      # LAR (Flye-based)
├── valid_sv/
│   ├── evidence/
│   │   ├── layer_depth.py             # Read-depth signature
│   │   ├── layer_kmer.py              # k-mer spectrum
│   │   ├── layer_breakpoint.py        # Split-read analysis
│   │   └── layer_ploidy.py            # SNV het rate
│   ├── engine/
│   │   ├── scorer.py                  # T-score calculation
│   │   └── fdr_estimator.py           # Mixture model FDR
│   ├── benchmarks/
│   │   ├── spike_in.py                # Truth set generator
│   │   ├── run_calibration.py         # Calibration runner
│   │   └── ablation.py                # Layer contribution study
│   ├── reporting/report_card.py       # Per-SV reports
│   └── run_validation.py              # Main entry point
├── config/config.yaml                 # Configuration
└── data/reference/                    # Reference genome + index
📚 Key References
Parameters and methods informed by:

Liu et al. (2024) Nature Communications — SV caller benchmarking, ICB overlap thresholds

Liu et al. (2024) Genome Biology — Multi-pipeline evaluation, merging strategies

Dunn et al. (2024) Genome Biology — Joint small+structural variant evaluation

Kronenberg et al. (2025) Nature Methods — Platinum Pedigree, SV merging

Hammond et al. (2025) Genome Research — HiFi validation

Chen et al. (2023) Nature Communications — DeBreak: local assembly for SVs

Helal et al. (2024) Scientific Reports — SV caller evaluation across aligners

Zhang et al. (2025) bioRxiv — SMaHT mosaic SV benchmark

Zheng & Shang (2024) PLOS ONE — SVvalidation: breakpoint validation method

Todd et al. (2025) Methods — SV-JIM: multi-caller consensus pipeline

Nkouamedjo et al. (2025) BMC Bioinformatics — SV-MeCa: ML-based meta-caller
Layer Weights: Uniform Priors
All evidence layers weighted equally (0.25). No published study provides empirical weights for combining orthogonal SV evidence types. Weights require spike-in calibration.
📄 Citation
Guimarães, K.H.A. et al. (2026). FUNGUS-SV: A triangulation-based structural variant discovery and validation pipeline for haploid genomes. In preparation.

📧 Contact
GitHub: @keltonjenkovguimaraes-alt

This pipeline is under active development. All parameters are documented with their sources or explicitly marked as uncalibrated. Any SV of biological interest must be independently validated by experimental methods.

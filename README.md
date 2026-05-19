# 🧬 FUNGUS-SV

**A structural variant discovery and triangulation-based prioritization pipeline for non-model haploid fungi using PacBio HiFi long reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)
[![Snakemake](https://img.shields.io/badge/Snakemake-6.15-blue.svg)](https://snakemake.github.io)

---

## ⚠️ HONEST ASSESSMENT — READ BEFORE USING

FUNGUS-SV is a **hypothesis-generation tool under active development.** It is NOT production-ready.

### What This Pipeline Actually Does

FUNGUS-SV detects structural variants (SVs ≥50 bp) in haploid fungal genomes using PacBio HiFi reads. It combines three SV callers (pbsv, Sniffles2, cuteSV) through an Intersection-Consensus-Builder (ICB), then validates each consensus SV using five orthogonal evidence layers:

| Layer | Method | What It Measures |
|-------|--------|-----------------|
| Local Assembly Refinement (LAR) | Flye assembly of reads at breakpoint | Confirms SV by assembling the variant allele |
| Read-Depth Signature | Coverage drop/increase at SV region | Detects deletions and duplications via copy number change |
| k-mer Spectrum | Jellyfish k-mer presence/absence | Confirms sequence gain/loss independent of alignment |
| Breakpoint Junction | Split-read and soft-clip analysis | Confirms precise breakpoint locations |
| Ploidy Confirmation | SNV heterozygosity rate | Verifies haploid assumption holds |

These layers are combined into a **T-score (0-1)** that ranks SVs by confidence.

### What We Know vs. What We Don't

| We Know | We Don't Know |
|---------|---------------|
| The pipeline compiles and runs end-to-end | Whether triangulation actually separates true from false SVs |
| ICB consensus (≥2 of 3 callers) reduces false positives | The optimal number of callers for fungal genomes |
| LAR via Flye can refine breakpoints | The false positive rate of LAR alone |
| Depth and k-mer are alignment-independent | Whether these layers are truly independent or correlated |
| The code is modular with 6 isolated conda environments | How the pipeline performs on real fungal sequencing data |

### First Calibration Results (Synthetic *C. auris* Data)

| Metric | Value |
|--------|-------|
| Genome | *Candidozyma auris* B11220 (RefSeq GCF_003013715.1, 12.25 Mb, 7 chromosomes) |
| Truth SVs | 50 spike-in SVs (20 DEL, 10 INS, 10 INV, 10 DUP) |
| Simulated reads | 47,365 HiFi-like reads, 58× coverage |
| Consensus SVs detected | 9 |
| True positives | 7 (all INV) |
| False positives | 2 |
| **Precision** | **77.8%** |
| **Recall** | **14.0%** |
| **F1** | **23.7%** |

**Caveat:** These results use a simple read simulator. Real HiFi data with proper SV breakpoints would likely yield different performance. This is a floor, not a ceiling.

### Known Limitations

1. **Layer weights are uniform (0.25 each).** No published study provides empirical weights for combining SV evidence layers. These are uninformative priors pending calibration.

2. **T-score thresholds are arbitrary.** We use 0.80/0.60/0.40 based on convention, not empirical FDR. The spike-in calibration has not been run at scale.

3. **63% of SVs are <100 bp** (from *S. schenckii* test data). These have reduced orthogonal validation — only k-mer and breakpoint layers apply.

4. **No fungal benchmark exists.** There is no GIAB-equivalent truth set for any fungus. All validation must use synthetic or orthogonal approaches.

5. **Not peer-reviewed.** Manuscript in preparation.

6. **Any SV of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing).**

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
│ Consensus clustering (≥2 callers) │
│ Output: candidate SVs │
└────────────────────────────────────────┘
│
▼
┌────────────────────────────────────────┐
│ PHASE 2: VALIDATION (VALID-SV) │
│ │
│ Layer 1: ICB agreement (reported only) │
│ Layer 2: Local Assembly (Flye) │
│ Layer 3: Read-Depth Signature │
│ Layer 4: k-mer Spectrum (Jellyfish) │
│ Layer 5: Breakpoint Junctions │
│ Layer 6: Ploidy Confirmation (SNV) │
│ ↓ │
│ Triangulation Engine → T-score │
└────────────────────────────────────────┘
│
▼
┌────────────────────────────────────────┐
│ PHASE 3: REPORTING │
│ Per-SV report cards + summary table │
└────────────────────────────────────────┘

---

## 🚀 Quick Start

### Prerequisites

- Linux (Ubuntu 20.04+)
- ≥32 GB RAM, ≥8 CPU cores (for real data; development possible with less)
- Conda or Mamba
- PacBio HiFi reads at ≥20× coverage (≥30× recommended per Hammond et al. 2025)

### Installation

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

# Create conda environments (one per tool group to avoid conflicts)
conda create -n sv_align -c bioconda -c conda-forge minimap2 samtools -y
conda create -n sv_call -c bioconda -c conda-forge sniffles=2.2 cutesv pbsv -y
conda create -n sv_valid -c conda-forge python=3.11 numpy scipy pandas pysam pyyaml biopython samtools -y
conda create -n sv_lar -c bioconda -c conda-forge python=3.10 flye samtools minimap2 -y
conda create -n sv_kmers -c bioconda -c conda-forge python=3.10 jellyfish -y
Usage
# 1. Align reads
conda activate sv_align
minimap2 -t 8 -ax map-hifi -R '@RG\tID:sample\tSM:sample' ref.fasta reads.fastq.gz \
    | samtools sort -@ 4 -o sample.sorted.bam -
samtools index sample.sorted.bam
conda deactivate

# 2. Run ICB consensus
conda activate sv_call
python fungus_sv/core/icb.py \
    --bam sample.sorted.bam \
    --reference ref.fasta \
    --output results/variants/ \
    --callers pbsv sniffles2 cutesv \
    --min-callers 2
conda deactivate

# 3. Run validation
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
│   │   └── run_calibration.py         # Calibration runner
│   ├── reporting/report_card.py       # Per-SV reports
│   └── run_validation.py              # Main entry point
├── config/config.yaml                 # Configuration
└── tests/                             # Test data
📚 Key References
Parameters and methods are informed by:

Liu et al. (2024) Nature Communications — SV caller benchmarking, ICB overlap thresholds

Liu et al. (2024) Genome Biology — Multi-pipeline evaluation, merging strategies

Dunn et al. (2024) Genome Biology — Joint small+structural variant evaluation (vcfdist)

Kronenberg et al. (2025) Nature Methods — Platinum Pedigree truth set, SV merging

Hammond et al. (2025) Genome Research — HiFi small variant validation

Chen et al. (2023) Nature Communications — DeBreak: local assembly for SV breakpoints

Helal et al. (2024) Scientific Reports — SV caller evaluation across aligners

Zhang et al. (2025) bioRxiv — SMaHT mosaic SV benchmark, binomial error model
Layer Weights: Uniform Priors
All evidence layers are weighted equally (0.25) because no published study provides empirical weights for combining orthogonal SV evidence types. This is a known gap in the literature. Weights will be calibrated via spike-in benchmarks in a future release.
📄 Citation
Guimarães, K.H.A; Philippsen H.K., et al. (2026). FUNGUS-SV. In preparation.

📧 Contact
Kelton Jenkov Guimarães — GitHub: @keltonjenkovguimaraes-alt
This README was written with the principle that honesty about limitations is more valuable than impressive-sounding claims. All parameters are documented with their sources or explicitly marked as uncalibrated.

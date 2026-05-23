# 🧬 FUNGUS-SV

**A structural variant discovery and triangulation-based prioritization pipeline for haploid genomes using PacBio HiFi reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)

---

## ⚠️ HONEST ASSESSMENT

FUNGUS-SV is a **hypothesis-generation tool under active development.** It is NOT production-ready for clinical use. All SVs of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing).

### What This Pipeline Actually Does

FUNGUS-SV detects structural variants (SVs ≥50 bp) in haploid genomes using PacBio HiFi reads. It combines multiple SV callers through an Intersection-Consensus-Builder (ICB), then validates each consensus SV using five orthogonal evidence layers. The output is a ranked list of SVs with confidence scores (T-scores), telling researchers **which SVs to prioritize for experimental validation.**

---

## 📊 Proven Performance

### Within-Species: A. baumannii ATCC 19606 vs. 6 Clinical Strains

| Strain | ICB SVs | HIGH (T≥0.6) | MED | WEAK | % HIGH |
|--------|---------|-------------|-----|------|--------|
| AB30 | 107 | 93 | 2 | 12 | 87% |
| MRSN15313 | 74 | 56 | 4 | 14 | 76% |
| DETAB-E51 | 68 | 59 | 0 | 9 | 87% |
| XH1056 | 71 | 65 | 0 | 6 | 92% |
| UC23022 | 73 | 61 | 0 | 12 | 84% |
| 6080 | 69 | 64 | 1 | 4 | **93%** |
| **TOTAL** | **462** | **398** | **7** | **57** | **86%** |

### Cross-Species: A. baumannii vs. 5 Acinetobacter spp.

| Species | ICB SVs |
|---------|---------|
| A. bouvetii | 5 |
| A. lwoffii | 8 |
| A. cumulans | 9 |
| A. lanii | 11 |
| A. larvae | 1 |
| **TOTAL** | **34** |

### Key Findings

1. **Within-species comparisons find 13.6× more SVs** than cross-species (462 vs. 34)
2. **86% of all consensus SVs score HIGH confidence** (T ≥ 0.6)
3. **100% of SVs ≥100 bp score HIGH; 100% of SVs <100 bp score WEAK**
4. **ICB consensus reduces raw calls by ~50%** — removes caller-specific false positives
5. **Strain 6080: 26 SVs with T=1.000** — all orthogonal layers active
6. **INV and DUP detected** where they exist (strain-dependent)
7. **Depth + breakpoint layers drive confidence** for large SVs

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
│ Parameters from Liu et al. 2024 │
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
│ Per-SV report cards │
│ Size-stratified summary │
│ SMaHT confidence tiers │
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

# Create isolated conda environments (no dependency conflicts)
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
    --callers sniffles2 cutesv \
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
├── workflow/                          # Snakemake workflow + env YAMLs
├── data/reference/                    # Reference genome + index
├── INTERNAL_BENCHMARK_METHODOLOGY.md  # Benchmark construction guide
└── RESULTS_SUMMARY.md                 # Detailed results
⚠️ Known Limitations
Limitation	Impact
Layer weights are uniform priors (0.25)	No published study provides empirical weights for combining SV evidence layers
T-score thresholds uncalibrated	T≥0.6 is "HIGH" by convention, not by empirical FDR
k-mer layer limited to same-species	Control k-mers absent in cross-species comparisons
Ploidy layer uses mpileup fallback	Longshot integration is unstable
LAR not run by default	Local assembly must be executed separately
No fungal benchmark exists	All validation uses bacterial data; fungal genomes pending
INS detection is limited	Insertions are the hardest SV type for all callers
63% of SVs <100 bp have reduced validation	Only breakpoint layer applies to small SVs
📚 Key References
Parameters and methods informed by 15+ peer-reviewed publications:

Liu et al. (2024) Nature Communications — SV caller benchmarking, ICB overlap thresholds

Liu et al. (2024) Genome Biology — Multi-pipeline evaluation, merging strategies

Dunn et al. (2024) Genome Biology — Joint SV evaluation (vcfdist)

Kronenberg et al. (2025) Nature Methods — Platinum Pedigree, SV merging

Hammond et al. (2025) Genome Research — HiFi validation benchmarks

Chen et al. (2023) Nature Communications — DeBreak: local assembly for SVs

Helal et al. (2024) Scientific Reports — SV caller evaluation across aligners

Zhang et al. (2025) bioRxiv — SMaHT mosaic SV benchmark

Zheng & Shang (2024) PLOS ONE — SVvalidation: breakpoint validation

Todd et al. (2025) Methods — SV-JIM: multi-caller consensus pipeline

Nkouamedjo et al. (2025) BMC Bioinformatics — SV-MeCa: ML meta-caller

📄 Citation
Guimarães, K.H.A. et al. (2026). FUNGUS-SV: A triangulation-based structural variant discovery and validation pipeline for haploid genomes. In preparation.

📧 Contact
GitHub: @keltonjenkovguimaraes-alt

This pipeline is under active development. All parameters are documented with their sources or explicitly marked as uncalibrated. Any SV of biological interest must be independently validated by experimental methods.

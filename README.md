# 🧬 FUNGUS-SV

**A structural variant discovery and triangulation-based prioritization pipeline for haploid genomes using PacBio HiFi reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)

---

## ⚠️ HONEST ASSESSMENT

FUNGUS-SV is a **hypothesis-generation tool under active development.** It is NOT production-ready for clinical use. All SVs of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing).

---

## 📊 Proven Performance

### Within-Species Validation: *A. baumannii* ATCC 19606 vs. 12 Clinical Strains

| Strain | ICB SVs | HIGH | MED | WEAK | % HIGH |
|--------|---------|------|-----|------|--------|
| AB30 | 107 | 93 | 2 | 12 | 87% |
| SRM25 | 85 | 77 | 1 | 7 | 91% |
| Aci4735 | 75 | 60 | 4 | 11 | 80% |
| MRSN15313 | 74 | 56 | 4 | 14 | 76% |
| XH1056 | 71 | 65 | 0 | 6 | 92% |
| UC23022 | 73 | 61 | 0 | 12 | 84% |
| 6080 | 69 | 64 | 1 | 4 | 93% |
| DETAB-E51 | 68 | 59 | 0 | 9 | 87% |
| AR_0083 | 68 | 51 | 0 | 17 | 75% |
| XH1037 | 64 | 57 | 1 | 6 | 89% |
| 966CSF | 60 | 47 | 1 | 12 | 78% |
| 280820 | 46 | 37 | 0 | 9 | 80% |
| **TOTAL** | **860** | **727** | **14** | **119** | **85%** |

### Cross-Species Validation: *A. baumannii* vs. 5 *Acinetobacter* spp.

| Species | ICB SVs |
|---------|---------|
| *A. lanii* | 11 |
| *A. cumulans* | 9 |
| *A. lwoffii* | 8 |
| *A. bouvetii* | 5 |
| *A. larvae* | 1 |
| **TOTAL** | **34** |

### Key Findings

| # | Finding |
|---|---------|
| 1 | **860 SVs detected across 12 strains** (mean: 72 per strain) |
| 2 | **85% of all SVs score HIGH confidence** (T ≥ 0.6) |
| 3 | **100% of SVs ≥100 bp score HIGH** — depth + breakpoint layers confirm |
| 4 | **100% of SVs <100 bp score WEAK** — limited orthogonal validation |
| 5 | **Within-species finds 25× more SVs** than cross-species (860 vs. 34) |
| 6 | **ICB consensus reduces raw calls by ~50%** — removes false positives |
| 7 | **INV and DUP detected** where they exist (strain-dependent) |
| 8 | **Depth + breakpoint layers drive confidence** for deletions ≥100 bp |

---

## 🏗️ Pipeline Architecture

PacBio HiFi reads (*.fastq.gz)
│
▼
[minimap2] ← map-hifi preset
│
▼
┌──────────────────────────────────────────┐
│ PHASE 1: PREDICTION (ICB) │
│ Sniffles2 + cuteSV + SVIM + pbsv │
│ ↓ │
│ Consensus: ≥2 callers, 0.5 overlap │
│ Parameters from Liu et al. (2024) │
└──────────────────────────────────────────┘
│
▼
┌──────────────────────────────────────────┐
│ PHASE 2: VALIDATION (VALID-SV) │
│ │
│ Layer 1: Local Assembly (Flye) │
│ Layer 2: Read-Depth Signature │
│ Layer 3: k-mer Spectrum (Jellyfish) │
│ Layer 4: Breakpoint Junctions │
│ Layer 5: Ploidy Confirmation │
│ ↓ │
│ Triangulation Engine → T-score (0-1) │
└──────────────────────────────────────────┘
│
▼
┌──────────────────────────────────────────┐
│ PHASE 3: REPORTING │
│ Per-SV report cards │
│ Size-stratified summary │
│ SMaHT confidence tiers │
└──────────────────────────────────────────┘

### Evidence Layers

| Layer | Method | What It Measures | Status |
|-------|--------|-----------------|--------|
| **Local Assembly** | Flye assembly at breakpoint | Confirms SV by assembling variant allele | ⚠️ Manual |
| **Read-Depth** | Coverage drop/increase ratio | Copy number change for DEL/DUP | ✅ Active |
| **k-mer Spectrum** | Jellyfish presence/absence | Sequence gain/loss independent of alignment | ⚠️ Same-species |
| **Breakpoint Junction** | Split-read + soft-clip analysis | Precise breakpoint locations | ✅ Active |
| **Ploidy Confirmation** | SNV heterozygosity rate | Verifies haploid assumption | ⚠️ mpileup |

---

## 🚀 Quick Start

### Prerequisites
- Linux (Ubuntu 20.04+)  |  ≥16 GB RAM  |  Conda/Mamba  |  PacBio HiFi reads (≥20×)

### Installation
```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

conda create -n sv_align -c bioconda -c conda-forge minimap2 samtools -y
conda create -n sv_call -c bioconda -c conda-forge sniffles=2.2 cutesv svim bcftools -y
conda create -n sv_valid -c conda-forge -c bioconda python=3.11 numpy scipy pandas pysam pyyaml samtools minimap2 -y
conda create -n sv_lar -c bioconda -c conda-forge python=3.10 flye samtools minimap2 -y
conda create -n sv_kmers -c bioconda -c conda-forge python=3.10 jellyfish -y
Usage
# 1. Align
conda activate sv_align
minimap2 -t 8 -ax map-hifi -R '@RG\tID:sample\tSM:sample' ref.fasta reads.fastq.gz \
    | samtools sort -@ 4 -o sample.sorted.bam -
samtools index sample.sorted.bam
conda deactivate

# 2. Detect SVs (ICB consensus)
conda activate sv_call
python fungus_sv/core/icb.py --bam sample.sorted.bam --reference ref.fasta \
    --output results/variants/ --callers sniffles2 cutesv --min-callers 2
conda deactivate

# 3. Validate (orthogonal evidence)
conda activate sv_valid
python -m valid_sv.run_validation --consensus-vcf results/variants/consensus_svs.vcf \
    --bam sample.sorted.bam --reference ref.fasta --fastq reads.fastq.gz \
    --output results/validation/
conda deactivate
⚠️ Known Limitations
Limitation	Detail
Layer weights are uniform	0.25 priors — no empirical calibration exists
T-score thresholds	Arbitrary convention, not FDR-calibrated
k-mer layer	Control k-mers absent in cross-species comparisons
Ploidy layer	Uses mpileup fallback; Longshot unstable
LAR	Must be run manually; not wired into pipeline
Small SVs (<100 bp)	Only breakpoint layer applies
Fungal data	Validated on bacteria; fungal testing pending
📚 References
Parameters informed by 13 peer-reviewed publications:

Liu et al. (2024) Nature Communications — ICB thresholds, breakpoint data

Liu et al. (2024) Genome Biology — Multi-pipeline evaluation

Dunn et al. (2024) Genome Biology — Joint SV evaluation (vcfdist)

Kronenberg et al. (2025) Nature Methods — SV merging strategy

Hammond et al. (2025) Genome Research — HiFi validation

Chen et al. (2023) Nature Communications — DeBreak LAR

Helal et al. (2024) Scientific Reports — Caller performance

Zhang et al. (2025) bioRxiv — SMaHT benchmark

Zheng & Shang (2024) PLOS ONE — SVvalidation method

Todd et al. (2025) Methods — SV-JIM pipeline

Nkouamedjo et al. (2025) BMC Bioinformatics — SV-MeCa meta-caller

Joe et al. (2024) BMC Genomics — Short-read caller data

Xing et al. (2025) BMC Genomics — k-mer ploidy method

📄 Citation
Guimarães, K.H.A. et al. (2026). FUNGUS-SV: A triangulation-based structural variant discovery and validation pipeline for haploid genomes. In preparation.

📧 Contact
GitHub: @keltonjenkovguimaraes-alt
All parameters documented with sources or marked uncalibrated. SVs require experimental validation.

# 🧬 FUNGUS-SV

**A structural variant discovery and triangulation-based prioritization pipeline for non-model haploid fungi using PacBio HiFi long reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)
[![Snakemake](https://img.shields.io/badge/Snakemake-8.0-blue.svg)](https://snakemake.github.io)

---

## ⚠️ CRITICAL: Read Before Using

**FUNGUS-SV is a hypothesis-generation tool under active development. It is NOT production-ready.**

### Honest Assessment of Current State

| Component | Status | Details |
|-----------|--------|---------|
| ICB consensus calling | ✅ Working | 3 callers, 275 SVs from *S. schenckii* |
| Depth signature layer | ✅ Working | For SVs ≥50 bp |
| k-mer spectrum layer | ✅ Working | Requires pre-built jellyfish DB |
| Breakpoint junction layer | ⚠️ Fixed, untested | Rewritten — needs validation run |
| LAR (local assembly) | ⚠️ Fixed, untested | Now runs Flye automatically |
| Ploidy confirmation | ✅ Working | Longshot SNV het rate analysis |
| T-score calculation | ⚠️ Revised | Excludes circular validation |
| FDR calibration | ❌ Not calibrated | Spike-in script exists, not run |
| Size-stratified reporting | ⚠️ Added | Not yet populated with data |

### Key Limitations

1. **No ground truth calibration exists.** The `spike_in` benchmark has not been run. T-score thresholds (0.80, 0.60, 0.40) are arbitrary until calibrated.

2. **63% of SVs are <100 bp** (173/275). These have limited orthogonal validation — only k-mer and breakpoint layers apply.

3. **Breakpoint layer was returning 0.0 for all SVs.** This has been fixed but NOT re-run on the data yet.

4. **LAR was a placeholder (score=0.5).** Now rewritten to run Flye assembly automatically, but UNTESTED on the full dataset.

5. **Strain mismatch:** NBRC32961 reads mapped to 1099-18 reference. SVs may include strain differences.

6. **Not peer-reviewed.** Manuscript in preparation.

**Any SV of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing).**

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
│ Layer 2: Local Assembly (Flye) │
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

### Current Data (from previous run — BEFORE fixes)

| Metric | Value |
|--------|-------|
| Platform | PacBio Revio |
| Coverage | 58× |
| Total SVs (≥2 callers) | 275 |
| 3-caller agreement | 174 |
| 2-caller agreement | 101 |
| SV size range | 29 bp – 139 kb |
| **SVs <100 bp** | **173 (62.9%)** |
| SVs 100-1000 bp | 83 (30.2%) |
| SVs >1000 bp | 19 (6.9%) |

### Validation Status by Size (BEFORE fixes)

| Size Range | Depth | k-mer | Breakpoint | LAR | Actually Validated? |
|------------|-------|-------|-------------|-----|---------------------|
| <100 bp | ❌ | ⚠️ ~0.6 | ❌ 0.0 | ❌ placeholder | k-mer only |
| 100-1000 bp | ✅ | ⚠️ ~0.6 | ❌ 0.0 | ❌ placeholder | depth + k-mer |
| >1000 bp | ⚠️ 79% | ⚠️ 74% | ❌ 0.0 | ❌ placeholder | depth + k-mer |

**After fixes applied in this commit:**
- Breakpoint layer rewritten (needs re-run)
- LAR integrated (needs re-run)
- Depth layer extended to ≥50 bp

---

## ⚡ Quick Start

### Prerequisites

- Linux (Ubuntu 20.04+)
- ≥32 GB RAM, ≥8 CPU cores
- Conda or Mamba
- PacBio HiFi reads at ≥30× coverage

### Installation

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv
conda env create -f workflow/envs/environment.yaml
conda activate snp_svant_pacbio
Prepare Data
mkdir -p data/raw data/reference
cp /path/to/your_sample.fastq.gz data/raw/
cp /path/to/reference.fasta data/reference/
samtools faidx data/reference/reference.fasta
minimap2 -d data/reference/reference.mmi data/reference/reference.fasta
# Phase 1: SV Prediction
snakemake -s workflow/Snakefile --cores 8

# Phase 2: SV Validation (now includes LAR)
python -m valid_sv.run_validation \
    --consensus-vcf results/variants/consensus/sample.consensus_svs.vcf \
    --bam results/alignment/sample.sorted.bam \
    --reference data/reference/reference.fasta \
    --fastq data/raw/sample.fastq.gz \
    --output results/validation/
python valid_sv/benchmarks/run_calibration.py \
    --reference data/reference/reference.fasta \
    --bam results/alignment/sample.sorted.bam \
    --fastq data/raw/sample.fastq.gz \
    --output results/calibration/
📁 Repository Structure
fungus-sv/
├── workflow/Snakefile                 # Main pipeline
├── workflow/envs/environment.yaml     # Conda environment
├── fungus_sv/
│   ├── core/icb.py                    # Consensus builder
│   └── modules/local_assembly.py      # LAR (Flye-based)
├── valid_sv/
│   ├── evidence/
│   │   ├── layer_depth.py             # Read-depth (≥50 bp)
│   │   ├── layer_kmer.py              # k-mer spectrum
│   │   ├── layer_breakpoint.py        # Split-read analysis
│   │   └── layer_ploidy.py            # SNV het rate
│   ├── engine/
│   │   ├── scorer.py                  # T-score (excludes consensus)
│   │   └── fdr_estimator.py           # Mixture model
│   ├── benchmarks/
│   │   ├── spike_in.py                # Truth set generator
│   │   └── run_calibration.py         # Calibration runner
│   ├── reporting/report_card.py       # Per-SV reports
│   └── run_validation.py              # Main entry point
├── config/config.yaml
└── results/validation_final/          # Previous run output
📚 References
Ko & Brandizzi (2025) — NEEDLE. Cell Reports Methods, 5, 100963.

Turner et al. (2026) — OrthoGarden. Mol Biol Evol, 43, msag053.

Gunasekaran et al. (2024) — SNP-SVant. Current Protocols, 4, e1046.

Liu et al. (2024) — SV detection tradeoffs. Nature Communications, 15, 2447.
📄 Citation
Guimarães, K.H.A; Philippsen H.K., et al. (2026). FUNGUS-SV. In preparation.
📧 Contact
Kelton Jenkov Guimarães — GitHub: @keltonjenkovguimaraes-alt

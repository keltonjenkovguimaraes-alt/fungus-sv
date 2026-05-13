# 🧬 FUNGUS-SV

**A structural variant discovery and triangulation-based prioritization pipeline for non-model haploid fungi using PacBio HiFi long reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)
[![Snakemake](https://img.shields.io/badge/Snakemake-8.0-blue.svg)](https://snakemake.github.io)

---

## ⚠️ Important: What This Pipeline Is and Is Not

**FUNGUS-SV is a hypothesis-generation tool, not a truth machine.**

There is no benchmark truth set for *Sporothrix schenckii* (or most non-model fungi). No external gold standard exists to validate structural variant calls against. What this pipeline provides:

| Claim | Honest Interpretation |
|-------|----------------------|
| "275 ICB consensus SVs" | 275 SVs where ≥2 of 3 alignment-based callers agree |
| "174 three-caller SVs" | 174 SVs where all 3 callers agree on the presence of an SV |
| "50 triple-triangulated SVs" | 50 SVs supported by ≥3 orthogonal evidence layers with uncorrelated failure modes |
| "T-score" | A weighted evidence score (0-1), NOT a probability. Higher = more independent support |
| "Estimated FDR" | Approximate, based on mixture modeling of T-score distribution. NOT calibrated |

**Any SV of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing) before functional interpretation.**

---

## 📖 Overview

Structural variant (SV) calling in non-model organisms faces a fundamental challenge: **how do you assess confidence when no benchmark truth set exists?** For humans, the Genome in a Bottle (GIAB) consortium provides gold-standard calls. For *Sporothrix schenckii*, *Candida auris*, *Aspergillus fumigatus*, and countless other medically important fungi — nothing.

FUNGUS-SV addresses this through a **two-phase architecture** inspired by NEEDLE (Ko & Brandizzi, 2025) and OrthoGarden (Turner et al., 2026):

| Phase | Component | Description |
|-------|-----------|-------------|
| 🔵 **Prediction** | ICB (Iterative Consensus Builder) | Runs 3 orthogonal SV callers (pbsv, Sniffles2, cuteSV), clusters overlapping SVs by reciprocal overlap, and scores by multi-caller agreement |
| 🟢 **Validation** | VALID-SV Triangulation Core | Combines 5 independent evidence layers with uncorrelated failure modes to estimate confidence without a truth set |
| 🟣 **Annotation** | Pan-Atlas | Maps SVs across species via ortholog matching and annotates functional impact |

---

## 🎯 Why Haploid Fungi?

FUNGUS-SV is explicitly designed for **haploid organisms**. This simplifies variant calling (no heterozygosity, no allelic phasing) and makes the pipeline faster and more interpretable. Compatible organisms include:

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
│ Layer 1: ICB caller agreement │
│ Layer 2: Local Assembly (LAR) │
│ Layer 3: Read-Depth Signature ← NEW │
│ Layer 4: k-mer Spectrum ← NEW │
│ Layer 5: Breakpoint Junctions ← NEW │
│ ↓ │
│ Triangulation Engine → T-score + FDR │
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

These results are from computational triangulation only. **They are not experimentally validated.**

| Triangulation Tier | T-score | 3-caller SVs (n=174) | 2-caller SVs (n=101) | All SVs (n=275) |
|-------------------|---------|---------------------|---------------------|-----------------|
| Triple-triangulated | ≥0.80 | 45 (25.9%) | 5 (5.0%) | 50 (18.2%) |
| Double-confirmed | 0.60–0.79 | 51 (29.3%) | 5 (5.0%) | 56 (20.4%) |
| Single-line evidence | 0.40–0.59 | 32 (18.4%) | 58 (57.4%) | 90 (32.7%) |
| Weak/contradicted | <0.40 | 46 (26.4%) | 33 (32.7%) | 79 (28.7%) |

**Key finding:** Of the 174 SVs where all three callers agree, 46 (26.4%) show weak or contradictory evidence from orthogonal layers. These are candidates for false positives or reference assembly errors that consensus-alone approaches would miss.

### Evidence Layer Performance

| Layer | Independence | Works best for | Limitation |
|-------|-------------|----------------|------------|
| Alignment consensus | Low (3 callers share alignment assumptions) | All SV types | Circular if used alone |
| Local assembly (LAR) | Medium (different algorithmic paradigm) | All types, ≥500 bp | Requires ≥10× local coverage |
| Read-depth signature | High (counts reads, not alignments) | DEL, DUP ≥100 bp | Not applicable for INV |
| k-mer spectrum | Highest (no reference, no alignment) | DEL, INS | Requires pre-built jellyfish DB |
| Breakpoint junctions | Medium (split-read analysis) | All types | Low sensitivity on HiFi data |

---

## ⚠️ Limitations (Please Read Before Using)

1. **No external truth set exists.** T-scores and FDR estimates are internally derived and APPROXIMATE. They have not been calibrated against experimental validation.

2. **Systematic false positives are possible.** All three SV callers share underlying assumptions (linear reference alignment). A reference assembly error could produce confident consensus calls that are all false. The triangulation layers catch some but not all of these.

3. **False negative rate is unknown.** We can only assess SVs that were called. SVs missed by all three callers are invisible to this pipeline.

4. **LAR was run on only 3 of 275 SVs in the original study.** The LAR layer in VALID-SV is currently a placeholder (score=0.5) until LAR is run on each SV. The T-scores reported here do NOT include LAR confirmation.

5. **Small SVs (<100 bp) have limited triangulation.** Depth signature requires ≥100 bp. k-mer resolution depends on k-mer size (31 bp). The smallest SVs rely primarily on caller agreement alone.

6. **k-mer analysis requires a pre-built jellyfish database** (~1.6 GB for 58× coverage of a 32 Mb genome). This is built once and cached.

7. **These results have not been peer-reviewed.** The manuscript is in preparation. Use with appropriate skepticism.

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

# Install additional dependency for k-mer validation
conda install -c bioconda jellyfish
data/raw/your_sample.fastq.gz    # PacBio HiFi reads
data/reference/reference.fasta   # Reference genome
data/reference/reference.gff     # Gene annotations
# Index reference
samtools faidx data/reference/reference.fasta
minimap2 -d data/reference/reference.mmi data/reference/reference.fasta

# Run full pipeline (prediction + validation)
snakemake -s workflow/Snakefile --cores 8
python -m valid_sv.run_validation \
    --consensus-vcf results/variants/consensus/sample.consensus_svs.vcf \
    --bam results/alignment/sample.sorted.bam \
    --reference data/reference/reference.fasta \
    --fastq data/raw/sample.fastq.gz \
    --output results/validation/
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
│       └── local_assembly.py         # LAR module
├── valid_sv/                         # Validation phase (NEW)
│   ├── evidence/
│   │   ├── layer_depth.py            # Read-depth signature
│   │   ├── layer_kmer.py             # k-mer spectrum analysis
│   │   └── layer_breakpoint.py       # Breakpoint junction analysis
│   ├── quality/
│   │   └── triangulability.py        # Assess which layers can validate each SV
│   ├── engine/
│   │   ├── scorer.py                 # Weighted T-score computation
│   │   └── fdr_estimator.py          # Mixture model FDR estimation
│   ├── reporting/
│   │   └── report_card.py            # Per-SV evidence summaries
│   └── run_validation.py             # Main validation entry point
├── config/config.yaml                # Pipeline configuration
└── docs/methods.md                   # Detailed methods
📚 References
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

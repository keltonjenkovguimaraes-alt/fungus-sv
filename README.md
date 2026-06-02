# 🧬 FUNGUS-SV

**Structural variant discovery and triangulation-based validation for haploid fungal genomes using PacBio HiFi reads.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)
[![Status](https://img.shields.io/badge/status-v0.9.3%20calibrated-blue.svg)]()

🌐 **Interactive Report:** [fungus-sv.netlify.app](https://fungus-sv.netlify.app)  
📦 **Zenodo:** [10.5281/zenodo.XXXXXXX](https://zenodo.org/) *(pending)*  
📄 **Paper:** In preparation (Guimarães et al., 2026)

---

## ⚠️ Disclaimer

FUNGUS-SV is a **hypothesis-generation tool under active development.** NOT for clinical use. All SVs of biological interest must be independently validated by experimental methods (PCR, Sanger sequencing).

---

## 📖 Overview

FUNGUS-SV detects structural variants (SVs) in haploid fungal genomes using PacBio HiFi long reads. It combines a **three-caller consensus** (Sniffles2 + cuteSV + SVIM) with a **five-layer triangulation validation system** that scores each SV based on orthogonal evidence: read depth, k-mer spectrum, breakpoint junctions, local assembly, and ploidy.

The pipeline was calibrated and validated on *Saccharomyces cerevisiae* CICC-1445 against five reference strains (S288C, BJ4, IMX2600, Makgeolli, SX2) using genome-wide depth analysis (DHFFC/DHBFC), split-read junction detection, manual curation with Samplot, and local assembly refinement (LAR).

---

## 🍺 Validation Results

### CICC-1445 vs 5 Reference Strains

**Data:** PacBio HiFi reads (SRR18210299, 274,915 reads, ~20 kb N50)  
**Callers:** Sniffles2 + cuteSV + SVIM (ICB consensus, ≥2 agreement)  
**Validation:** 5-layer triangulation + orthogonal depth/split-read analysis + LAR

#### SV Detection Summary

| Strain | Type | Total SVs | DEL | DUP | INV |
|--------|------|-----------|-----|-----|-----|
| BJ4 | Industrial | **165** | 140 | 13 | 12 |
| Makgeolli | Fermentation | 250 | 225 | 17 | 8 |
| S288C | Laboratory | 277 | 248 | 18 | 11 |
| SX2 | Industrial | 290 | 261 | 18 | 11 |
| IMX2600 | Engineered | **314** | 285 | 9 | 20 |

- **CICC-1445 is closest to BJ4 (165 SVs)** — both Chinese industrial strains
- **Most divergent from IMX2600 (314 SVs)** — engineered laboratory strain

#### v2.0 Calibrated Validation (All 5 Strains)

After integrating DHBFC (GC-corrected depth), size-stratified scoring, and haploid-specific thresholds:

| Strain | Total | HIGH (T≥0.6) | WEAK | CONTRADICTED |
|--------|-------|-------------|------|-------------|
| S288C | 277 | 35 (12.6%) | 92 | 149 |
| BJ4 | 165 | 30 (18.1%) | 61 | 71 |
| IMX2600 | 314 | 43 (13.6%) | 113 | 157 |
| Makgeolli | 250 | 31 (12.4%) | 91 | 124 |
| SX2 | 290 | 49 (16.8%) | 105 | 135 |

#### Orthogonal Validation Highlights

#### LAR Calibration (64 DELs Tested Across 5 Strains)

| Strain | Tested | Real | False | % Real |
|--------|:---:|:---:|:---:|:---:|
| S288C | 15 | 8 | 7 | 53% |
| BJ4 | 1 | 0 | 1 | 0% |
| IMX2600 | 28 | 5 | 23 | 18% |
| Makgeolli | 16 | 3 | 13 | 19% |
| SX2 | 4 | 1 | 3 | 25% |
| **Total** | **64** | **17** | **47** | **27%** |

**Key finding:** Only 27% of CONTRADICTED large DELs are real. DHFFC alone cannot predict reality across different reference genomes. LAR is recommended for any CONTRADICTED call of biological interest.

| Finding | Evidence |
|---------|----------|
| **11/11 inversions confirmed** by split reads | 32–1,091 junction reads per breakpoint |
| **115/248 deletions confirmed** by depth (DHFFC < 0.3) | 97% depth drop for rDNA deletion |
| **DHFFC + DHBFC show 90% agreement** | GC bias minimal in yeast genome |
| **Larger DELs more reliable** | 75.6% confirmed for 5–10 kb vs 33.3% for 50–100 bp |

---

## 🔬 LAR: Local Assembly Refinement (Standalone)

For **definitive proof** of an SV, LAR assembles reads from a single SV region using Flye. Unlike whole-genome assembly, LAR uses <500 MB RAM and completes in 2–15 minutes.

```bash
conda activate sv_lar
python valid_sv/evidence/layer_lar.py
LAR Validation Results (30 May 2026)
Tested on CICC-1445 against S288C, BJ4, and SX2 references:

Locus	Type	Size	Strain	Verdict	Reads	RAM	Time
YBL005W-B (Ty2)	DEL	5.9 kb	S288C	✅ confirmed	202	<200 MB	2 min
chrVII multi-gene	DEL	55.7 kb	S288C	⚠️ partial	581	<300 MB	4 min
SX2 chrII	INV	430 kb	SX2	✅ confirmed	6,448	<500 MB	13 min
BJ4 chrXII	INV	205 kb	BJ4	✅ confirmed	3,552	<400 MB	10 min
How it works: LAR extracts reads mapping to the SV region ± 3 kb, assembles them with Flye, aligns the contig back to the reference with minimap2, and parses the CIGAR string for deletion gaps or opposite-strand alignment (inversion). An assembled contig spanning a breakpoint is the computational equivalent of a PCR product.

chrVII note: The 55.7 kb region contains multiple deletions; LAR found the largest single gap (23.9 kb). The 35 affected genes include TIF4631, GTR2, ERG1, BUB1, UBR1, and CRH1.

🏗️ Pipeline Architecture
PacBio HiFi reads (*.fastq.gz)
│
▼
[minimap2] ← map-hifi preset
│
▼
┌─────────────────────────────────────────┐
│ PHASE 1: ICB CONSENSUS │
│ │
│ Sniffles2 │ cuteSV │ SVIM │
│ └────────┴─────────┘ │
│ │ │
│ Consensus: ≥2 callers, 0.5 overlap │
└─────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────┐
│ PHASE 2: TRIANGULATION VALIDATION │
│ │
│ Layer 1: Local Assembly (Flye) 0.20 │
│ Layer 2: Depth (DHFFC+DHBFC) 0.35 │
│ Layer 3: k-mer Spectrum (Jellyfish) 0.15 │
│ Layer 4: Breakpoint Junction 0.30 │
│ Layer 5: Ploidy Confirmation FILTER │
│ │
│ T = Σ(score × weight) / Σ(weights) │
└─────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────┐
│ PHASE 3: CONFIDENCE TIERS (v0.9.3) │
│ │
│ TRIPLE_TRIANGULATED: T ≥ 0.80 │
│ DOUBLE_CONFIRMED: T ≥ 0.60 │
│ SINGLE_LINE: T ≥ 0.40 │
│ WEAK: T ≥ 0.20 │
│ CONTRADICTED: T < 0.20 │
└─────────────────────────────────────────┘

Calibrated Weights (Haploid Fungi)
Layer	Original	Calibrated	Rationale
Depth (DHFFC + DHBFC)	0.25	0.35	100% depth loss unambiguous in haploids
Breakpoint Junction	0.20	0.30	Only signal for inversions; critical for all
Local Assembly	0.30	0.20	Computationally expensive; standalone tool available
k-mer Spectrum	0.25	0.15	Redundant with depth in haploids
Haploid DHFFC Thresholds
SV Type	Human Diploid	Haploid Fungi	Source
Deletion	< 0.7	< 0.3	Pedersen & Quinlan (2019); this study
Duplication	> 1.3	> 2.0	Pedersen & Quinlan (2019); this study
Size-Stratified Scoring
SV Size	Factor	Basis
≥ 5,000 bp	1.00	AUC ~1.0 (Pedersen & Quinlan 2019)
1,000–4,999 bp	0.95	AUC ~0.97
500–999 bp	0.85	Moderate penalty
100–499 bp	0.75	Significant penalty
< 100 bp	0.60	Heavy penalty
🚀 Quick Start
Prerequisites
Linux (Ubuntu 20.04+) | ≥16 GB RAM | Conda/Mamba

PacBio HiFi reads (≥20× coverage)

Installation
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv

# Create environments
conda create -n sv_align -c bioconda -c conda-forge minimap2 samtools -y
conda create -n sv_call -c bioconda -c conda-forge sniffles=2.2 cutesv svim bcftools -y
conda create -n sv_valid -c conda-forge -c bioconda python=3.11 numpy scipy pandas pysam pyyaml matplotlib -y
conda create -n sv_lar -c bioconda flye minimap2 samtools -y
Usage
# 1. Align reads to reference
conda activate sv_align
minimap2 -t 8 -ax map-hifi -R '@RG	ID:sample	SM:sample' ref.fasta reads.fastq.gz     | samtools sort -@ 4 -o sample.sorted.bam -
samtools index sample.sorted.bam

# 2. Detect SVs (ICB consensus — 3 callers)
conda activate sv_call
python fungus_sv/core/icb.py --bam sample.sorted.bam --reference ref.fasta     --output results/ --threads 4

# 3. Validate (orthogonal evidence triangulation)
conda activate sv_valid
PYTHONPATH=. python -m valid_sv.run_validation     --consensus-vcf results/consensus_svs.vcf     --bam sample.sorted.bam --reference ref.fasta --fastq reads.fastq.gz     --output results/validation/ --threads 4

# 4. LAR: Definitive proof for key SVs (standalone)
conda activate sv_lar
python valid_sv/evidence/layer_lar.py
📚 Parameter Sources
Parameter	Value	Source
ICB min_overlap	0.5	Liu et al. (2024) Genome Biology
ICB flank	200 bp	Kronenberg et al. (2025) Nature Methods
ICB min_callers	2	Liu et al. (2024) Genome Biology
DHFFC method	—	Pedersen & Quinlan (2019) GigaScience
DHBFC method	—	Pedersen & Quinlan (2019) GigaScience
Haploid DEL threshold	< 0.3	This study (calibrated from CICC-1445)
Haploid DUP threshold	> 2.0	This study (calibrated from CICC-1445)
Size stratification	0.60–1.00	Pedersen & Quinlan (2019); this study
distance_support	0.2×len + 2000/len	Zheng & Shang (2024) PLOS ONE
MAPQ filter	≥ 20	Zheng & Shang (2024); Liu et al. (2024)
k-mer size	31	PAV (Ebert 2021), SV-JIM (Todd 2025)
Depth min size	100 bp	Liu et al. (2024) Nature Communications
Haploid max het	7%	Xing et al. (2025) BMC Genomics
Confidence tiers	T≥0.80/0.60/0.40/0.20	SMaHT (Zhang et al. 2025) bioRxiv
Weights (calibrated)	0.20/0.35/0.15/0.30	Liu et al. (2024); this study
Inversion scoring	Breakpoint-only	This study (11/11 INVs confirmed)
LAR flank size	3,000 bp	This study (tested on 4 SVs)
📊 Visualization & Reports
Interactive HTML report: fungus-sv.netlify.app — genome tracks, gene annotations, validation dashboard

Samplot images: Manual curation of candidate SVs (Belyeu et al. 2021)

Comparative Samplot: Same SV region across all 5 references (YBL005W-B, chrVII)

DHFFC/DHBFC plots: Size-stratified depth validation

Split-read bar charts: Inversion breakpoint evidence

⚠️ Known Limitations
Limitation	Detail
Inversions	Now reported as INV_SPLIT_READ_CONFIRMED when split-reads support (11/11 confirmed)
Duplications	Only 2/18 confirmed by depth; most score CONTRADICTED
Small SVs (<100 bp)	Size factor 0.60; depth signal unreliable
Repetitive regions	FLO genes, rDNA, Ty elements produce complex signals
Cross-species	k-mer layer fails (control k-mers absent in distant species)
No spike-in calibration	64 LAR-validated DELs provide empirical calibration; spike-in pending
| DHFFC cross-reference	DHFFC thresholds are reference-dependent; 0.000 DHFFC can be FALSE in some strains |
📄 Citation
Guimarães, K.H.A. et al. (2026). FUNGUS-SV: A triangulation-based structural variant discovery and validation pipeline for haploid genomes. In preparation.

📧 GitHub: @keltonjenkovguimaraes-alt

🔬 Previous Validation: Acinetobacter baumannii (Bacterial)
The pipeline was originally validated on 12 clinical A. baumannii strains against ATCC 19606 reference.

Metric	Bacteria (A. baumannii)	Fungi (S. cerevisiae)
Strains tested	12 clinical	5 reference
Mean SVs/strain	72	259
HIGH confidence rate	85%	12–18% (v2 calibrated)
Dominant SV type	DEL	DEL
Genome size	~4 Mb	~12 Mb
Reads	PacBio HiFi	PacBio HiFi
Callers	Sniffles2 + cuteSV	Sniffles2 + cuteSV + SVIM
📁 Repository Structure
fungus-sv/
├── config/                  # Pipeline configuration
├── fungus_sv/core/          # ICB consensus calling
├── valid_sv/                # Triangulation validation
│   ├── evidence/            # Layer implementations (layer_lar.py = standalone LAR)
│   ├── engine/              # Scoring engine
│   └── reporting/           # Report generation
├── data/yeast/              # Validation data (CICC-1445 vs 5 references)
├── docs/                    # Documentation & figures
├── figures/                 # Generated plots, HTML & Samplot images
└── workflow/envs/           # Conda environment specs
📖 References
Pedersen BS, Quinlan AR. Duphold. GigaScience. 2019;8(4):giz040.

Zheng Y, Shang X. SVvalidation. PLOS ONE. 2024;19(1):e0291741.

Belyeu JR, et al. Samplot. Genome Biology. 2021;22(1):161.

David G, et al. Manual curation of SVs. GBE. 2024;16(4):evae049.

Dhakal U, et al. SVs in Fusarium graminearum. G3. 2024;14(6):jkae065.

Li X, et al. Variant calling in Candida auris. Microbial Genomics. 2023;9(4):000979.

Liu Y, et al. SV detection evaluation. Nature Communications. 2024.

Liu Y, et al. Multi-pipeline SV evaluation. Genome Biology. 2024.

Zhang Y, et al. SMaHT benchmark. bioRxiv. 2025.

# FUNGUS-SV: Structural Variant Discovery for Non-Model Haploid Fungi

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-yellow.svg)](https://python.org)

> **A benchmark-free pipeline for structural variant discovery and validation from PacBio HiFi long reads, optimized for haploid fungal genomes.**

---

## Why FUNGUS-SV?

Most SV callers require benchmark truth sets (like GIAB for humans). For non-model organisms — including most pathogenic fungi — no benchmarks exist. FUNGUS-SV solves this through three innovations:

| Innovation | Description |
|-----------|-------------|
| **ICB** | Iterative Consensus Builder — 3 orthogonal SV callers, scored by agreement |
| **LAR** | Local Assembly Refinement — base-pair precision via local de novo assembly |
| **Pan-Atlas** | Cross-species SV conservation analysis |

---

## Pipeline

PacBio HiFi reads → minimap2 → [pbsv + Sniffles2 + cuteSV] → ICB → LAR → Pan-Atlas

---

## Validation: Sporothrix schenckii

| Metric | Value |
|--------|-------|
| Reads | 197,830 PacBio Revio HiFi |
| Coverage | 58× |
| Mapping rate | 99.91% |
| Consensus SVs | 275 |
| High-confidence | 174 |
| Genes affected | 83 |
| Cross-species conservation | 86.7% |

---

## Quick Start

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/fungus-sv.git
cd fungus-sv
conda env create -f workflow/envs/environment.yaml
conda activate snp_svant_pacbio
snakemake --cores 8
Citation
Guimarães, K.H.A;Philippsen, H.k; et al. (2026). FUNGUS-SV: A benchmark-free structural variant discovery pipeline for non-model haploid fungi. In preparation.
License
MIT

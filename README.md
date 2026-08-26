# FUNGUS-SV
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21402263.svg)](https://doi.org/10.5281/zenodo.21402263)
Structural variant detection and validation for haploid fungal genomes using PacBio HiFi reads.

## Version

v1.1

## Overview

FUNGUS-SV calls structural variants (SVs) with a three-caller ensemble (Sniffles2, cuteSV, SVIM), merges them into an ICB consensus, and validates each call through a multi-layer evidence pipeline.

The automated scoring layer is depth-based. Breakpoint presence is treated as a confirmatory flag, not a scored weight. DUP and INV calls are reclassified through overrides that correct the depth layer's known blindness to balanced events in haploid genomes.

LAR (local assembly reconstruction) is available for manual validation of key candidates and uses a two-assembler consensus (Flye + Miniasm).

## Pipeline structure

fungus_sv/core/ ICB consensus building
valid_sv/
run_validation.py main validation entry point
engine/ scoring and tier assignment
evidence/ evidence layers
layer_depth.py
layer_breakpoint.py
layer_lar.py
layer_ploidy.py
layer_genomic_context.py
layer_backtrack.py
benchmarks/ synthetic and spike-in evaluation
config/ pipeline configuration
workflow/ Snakefile and conda envs
docs/ documentation

## Evidence layers

| Layer | Role | Status |
|-------|------|--------|
| ICB consensus | Caller agreement | Informational |
| LAR | Local assembly validation | Manual, key candidates |
| Depth signature | Primary scored layer | Active (weight 0.65) |
| k-mer spectrum | Optional auxiliary | Not run in this study |
| Breakpoint junction | Confirmatory flag | Not scored |
| Ploidy | Haploid consistency hard filter | Active |
| Genomic context | Gene annotation hard filter | Active |
| Backtrack | Read- and reference-based SV validation | Standalone / integrative |

## Backtrack validation

The backtrack layer supports two modes:

- **Read-based**: uses split-read orientation and depth from a BAM file.
- **Reference-based**: simulates the SV in silico on the reference assembly and realigns it with minimap2. Works without reads.

For INV validation, the reference-based mode reverse-complements the candidate region and aligns it back to the reference. A reverse-strand alignment indicates a true inversion.

**Status: heuristic thresholds, not yet calibrated against the full truth set.**

The current INV verdict thresholds were set from a small simulation study (11 simulated inversions) and the CEN5 centromeric inversion. They are reasonable defaults but have not been validated against the complete synthetic truth set or the LAR truth set. Treat backtrack verdicts as advisory, not authoritative.

Current thresholds:

CONFIRMED: reverse strand + identity >= 0.95 + MAPQ >= 50
PARTIAL: reverse strand + identity >= 0.80
AMBIGUOUS: reverse strand + identity < 0.80
CONTRADICTED: forward strand
UNCALLABLE: SV smaller than 500 bp

## Synthetic benchmark evaluation

`valid_sv/benchmarks/evaluate_synthetic_replicates.py` performs one-to-one greedy matching between consensus calls and truth SVs. Each truth SV can be claimed by at most one consensus call, and pooled and per-type metrics are derived from the same match list.

This fixes the double-matching bug in the original `run_calibration.py` evaluator, where two consensus calls near the same truth SV could both be counted as true positives.

## Validation results

### Synthetic benchmark (3 replicates)

| Metric | Value |
|--------|-------|
| Truth SVs | 243 |
| Detected SVs | 165 |
| True positives | 159 |
| False positives | 6 |
| False negatives | 84 |
| Recall | 65.4% |
| Precision | 96.4% |

Per-type recall:

| Type | Recall |
|------|--------|
| DEL | 79.5% |
| DUP | 44.4% |
| INV | 60.3% |

### Self-alignment specificity

| Strain | TRIPLE-tier FDR (v1.1) | 95% CI |
|--------|------------------------|--------|
| CICC-1445 | 0.8% (1/130) | 0.1–4.2% |
| SC5314 | 0.0% (0/200) | 0.0–1.9% |

### LAR truth set (180 SVs, 12 strains)

| Category | n | Confirmed (strict) | Raw rate | Excl. technical |
|----------|---|-------------------|----------|-----------------|
| Within-species | 151 | 78 | 51.7% | 69.6% |
| Self-alignment | 20 | 10 | 50.0% | 76.9% |
| Oggenfuss external | 7 | 5 | 71.4% | — |

### Spike-in benchmark

| Metric | Value |
|--------|-------|
| Recall | 76.2% (16/21) |

## Known limitations

- **DUP recall is low (44.4%)**. DUPs are the hardest SV type for haploid depth-based scoring. Treat DUP calls with more caution than DEL calls.
- **INV calls under 2 kb are unreliable**. The backtrack layer flags INVs below 500 bp as uncallable. INV detection in general remains challenging (60.3% recall).
- **Backtrack thresholds are heuristic**, not empirically calibrated. See the Backtrack section above.
- **The 5× coverage flag was heuristic**. It has since been replaced by the Li (2014) max-depth filter. See the v1.1 audit for details.

## Installation

The pipeline is distributed as a Docker image containing five conda environments:

sv_align minimap2, samtools
sv_call sniffles2, cutesv, svim, bcftools
sv_valid python, numpy, scipy, pandas, pysam, matplotlib
sv_lar flye, minimap2, samtools, miniasm, racon
sv_kmers jellyfish

## Running the validation pipeline

```bash
python -m valid_sv.run_validation \
  --consensus-vcf results/consensus_svs.vcf \
  --bam reads_vs_ref.bam \
  --reference reference.fasta \
  --output validation_output/
License
MIT. See LICENSE.

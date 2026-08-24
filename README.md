# FUNGUS-SV

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

For INV validation, the reference-based mode reverse-complements the candidate region and aligns it back to the reference. A reverse-strand alignment with high identity and MAPQ indicates a true inversion.

CONFIRMED: reverse strand + identity >= 0.95 + MAPQ >= 50
PARTIAL: reverse strand + identity >= 0.80
AMBIGUOUS: reverse strand + identity < 0.80
CONTRADICTED: forward strand
UNCALLABLE: SV smaller than 500 bp

## Synthetic benchmark evaluation

`valid_sv/benchmarks/evaluate_synthetic_replicates.py` performs one-to-one greedy matching between consensus calls and truth SVs. Each truth SV can be claimed by at most one consensus call, and pooled and per-type metrics are derived from the same match list.

This fixes the double-matching bug in the original `run_calibration.py` evaluator, where two consensus calls near the same truth SV could both be counted as true positives.

## Key results

Current synthetic benchmark results across three replicates:

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

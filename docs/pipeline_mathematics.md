# FUNGUS-SV Mathematical Framework

## Overview

The pipeline uses a weighted triangulation scoring system with five evidence layers,
each producing a score from 0.0 (contradicts) to 1.0 (strongly supports).

---

## Layer 1: Depth Signature (weight = 0.35)

### Formula: DHFFC (Duphold Flank Fold-Change)

depth_ratio = median(depth_SV_region) / median(depth_flank_region)

### Scoring (Haploid-Calibrated)

**Deletions:**
If depth_ratio < 0.3: score = min(1.0, (0.25 - depth_ratio) × 4 + 0.8) → [0.8, 1.0]
If 0.3 ≤ depth_ratio < 0.7: score = 0.6 + (0.5 - depth_ratio) × 1.6 → [0.28, 0.92]
If depth_ratio > 0.80: score = 0.0
Else: score = 0.3

**Duplications:**
If depth_ratio > 2.0: score = min(1.0, (depth_ratio - 1.0) × 0.8) → [0.8, 1.0]
If 1.3 < depth_ratio ≤ 2.0: score = 0.5 + (depth_ratio - 1.2) × 1.67 → [0.67, 1.0]
If depth_ratio < 0.80: score = 0.0
Else: score = 0.3

**Source:** Pedersen & Quinlan (2019), *GigaScience*; thresholds adapted for haploid fungi based on empirical data (CICC-1445 vs S288C)

---

## Layer 2: Breakpoint Junction (weight = 0.30)

### Split-Read Detection

split_reads = count of reads with CIGAR containing 'S' at SV breakpoints (±100 bp window)
total_reads = all reads mapped to breakpoint window
support_ratio = split_reads / total_reads

### Size Tolerance (from Zheng & Shang 2024)

distance_support = 0.2 × SV_length + 2000 / SV_length

A read supports the SV if:
abs(observed_SV_length - called_SV_length) ≤ distance_support

### Scoring

score = f(support_ratio, split_reads)
- Requires ≥3 split reads at each breakpoint for minimum score
- Score increases with proportion of supporting reads
- ≥80% support → score ≥ 0.8

**Source:** Zheng & Shang (2024), *PLOS ONE*; Belyeu et al. (2021), *Genome Biology*

---

## Layer 3: k-mer Spectrum (weight = 0.15)

### Method

k-mer size = 31
Jellyfish count on reads → compare k-mer frequencies in SV region vs flanking

### Scoring

If k-mer pattern matches expected for SV type → score = 0.7–1.0
If k-mer pattern ambiguous → score = 0.3–0.6
If k-mer pattern contradicts → score = 0.0–0.2

**Source:** PAV (Ebert 2021); SV-JIM (Todd 2025); Liu et al. (2024)

---

## Layer 4: Local Assembly (weight = 0.20)

### Method

Flye assembly on reads overlapping the SV region
→ compare assembled contig to reference

### Scoring

If assembly confirms SV → score = 0.8–1.0
If assembly partial → score = 0.4–0.7
If assembly contradicts → score = 0.0–0.3
If assembly not run → layer unavailable

**Source:** Liu et al. (2024), *Nature Communications*; DeBreak (Chen et al. 2023)

---

## Layer 5: Ploidy Confirmation (Hard Filter)

### Method

heterozygosity_rate = heterozygous_SNPs / total_SNPs

If heterozygosity_rate > 0.07 → FLAG (possible diploid/contamination)

**Source:** Xing et al. (2025), *BMC Genomics*

---

## Triangulation Engine

### T-Score Formula

T = Σ(layer_score_i × weight_i) / Σ(weight_i)

Where `i` iterates over all available layers.

### Completeness Penalty

If completeness < 0.5:
T = T × completeness

completeness = available_layers / total_possible_layers

### Score Uncertainty

σ = std_dev(all_layer_scores)

### Inversion-Specific Handling

Inversions skip the depth and k-mer layers (both structurally silent for balanced inversions):

If SV_type == "INV":
available_layers = [breakpoint_junction, local_assembly] # only these apply
No completeness penalty applied

### Confidence Tiers

TRIPLE_TRIANGULATED: T ≥ 0.80
DOUBLE_CONFIRMED: T ≥ 0.60
SINGLE_LINE: T ≥ 0.40
WEAK: T ≥ 0.20
CONTRADICTED: T < 0.20

**Source:** Zhang et al. (2025), *bioRxiv* (SMaHT convention)

---

## FDR Estimation

estimated_fdr = f(T_score)
- T ≥ 0.80 → FDR < 5%
- T ≥ 0.60 → FDR 5–20%
- T ≥ 0.40 → FDR 20–50%
- T ≥ 0.20 → FDR > 50%
- T < 0.20 → FDR > 90%

**Note:** FDR estimates are uncalibrated. Actual FDR from orthogonal validation:
- DEL: ~45% overall (T ≥ 0.60 subset: lower)
- INV: ~0% (all confirmed by split reads regardless of T-score)
- DUP: ~80% overall

---

## Layer Weights

| Layer | Original Weight | Calibrated Weight | Rationale |
|-------|----------------|-------------------|-----------|
| Depth | 0.25 | **0.35** | Haploid: 100% depth loss is unambiguous |
| Breakpoint | 0.20 | **0.30** | Only signal for inversions; critical for all types |
| Assembly | 0.30 | **0.20** | Optional; computationally expensive |
| k-mer | 0.25 | **0.15** | Redundant with depth in haploids |

**Source:** Liu et al. (2024) original framework; recalibrated from CICC-1445 vs S288C empirical validation

---

## References

1. Pedersen BS, Quinlan AR. Duphold. *GigaScience*. 2019;8(4):giz040.
2. Zheng Y, Shang X. SVvalidation. *PLOS ONE*. 2024;19(1):e0291741.
3. Belyeu JR, et al. Samplot. *Genome Biology*. 2021;22(1):161.
4. Liu Y, et al. Multi-pipeline evaluation. *Genome Biology*. 2024.
5. Liu Y, et al. ICB thresholds. *Nature Communications*. 2024.
6. Zhang Y, et al. SMaHT benchmark. *bioRxiv*. 2025.
7. Ebert P, et al. PAV. *Science*. 2021.
8. Todd J, et al. SV-JIM. *Methods*. 2025.
9. Chen Y, et al. DeBreak. *Nature Communications*. 2023.
10. Xing Y, et al. k-mer ploidy. *BMC Genomics*. 2025.

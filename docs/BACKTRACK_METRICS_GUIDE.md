# Backtrack Layer — Metrics Guide

## Overview

The backtrack layer reports raw depth metrics across a candidate SV region and its flanks. It makes **no verdicts** — it only measures and reports numbers. The triangulation scorer interprets these numbers.

## How It Works

Left Flank (1kb) | SV Region (variable) | Right Flank (1kb)
─────────────────────────────────────────────────────────────────────
Depth: ~350x | Depth: drops or changes | Depth: ~350x



Three `samtools depth` calls fetch per-base depth for:
1. Left flank (default 1000 bp before SV start)
2. SV region (from start to end)
3. Right flank (default 1000 bp after SV end)

Breakpoint transitions are computed by slicing the already-fetched arrays — no extra calls.

---

## Metrics Dictionary

### Flank Metrics
| Metric | What it is | Example | Meaning |
|--------|-----------|---------|---------|
| `left_flank_mean` | Average depth in left flank | 372.2 | Normal coverage left of the SV |
| `left_flank_median` | Median depth in left flank | 373.0 | Robust to outliers |
| `right_flank_mean` | Average depth in right flank | 361.9 | Normal coverage right of the SV |
| `right_flank_median` | Median depth in right flank | 362.0 | Robust to outliers |

**What to look for:**
- Symmetric flanks (~equal) = clean region
- Asymmetric flanks (very different) = possible repeat, CNV, or contig boundary
- One flank < 10x = SV at contig edge or assembly gap

---

### SV Region Metrics
| Metric | What it is | Example (real DEL) | Example (false DEL) |
|--------|-----------|-------------------|---------------------|
| `sv_region_mean` | Average depth inside SV | 97.1 | 327.5 |
| `sv_region_median` | Median depth inside SV | 98.0 | 326.0 |
| `sv_region_p10` | 10th percentile depth | 89.0 | 255.0 |
| `sv_region_p90` | 90th percentile depth | 103.0 | 387.0 |

**What to look for:**
- **DEL**: median much lower than flanks; p10 near zero
- **DUP**: median ~2x flanks
- **INV**: median ~same as flanks (balanced event)
- p10 vs p90 spread: wide spread = uneven depth (repeats, noise); narrow spread = clean signal

---

### Depth Ratios
| Metric | Formula | Example (real DEL) | Meaning |
|--------|---------|-------------------|---------|
| `mean_ratio` | sv_mean / flank_mean | 0.265 | SV region has 26.5% of flank depth |
| `median_ratio` | sv_median / flank_median | 0.267 | More robust: SV region has 26.7% of flank depth |

**What to look for:**
- **DEL**: ratio < 0.5 (depth dropped by >50%)
- **Strong DEL**: ratio < 0.2 (depth dropped by >80%)
- **DUP**: ratio > 1.5 (depth increased by >50%)
- **No SV / false call**: ratio ~1.0 (no change)
- Use `median_ratio` as primary — it's robust to outlier spikes

---

### Sparsity Metrics
| Metric | What it is | Example (real DEL) | Example (false) |
|--------|-----------|-------------------|-----------------|
| `zero_fraction` | % of bases with depth < 5 | 0.444 (44.4%) | 0.000 (0%) |
| `low_fraction` | % of bases with depth < 10 | 0.511 (51.1%) | 0.000 (0%) |

**What to look for:**
- **DEL**: high zero_fraction (>30%) = strong signal. Reads cannot align because sequence is absent.
- **DEL with long reads**: zero_fraction may be low even for real DELs. PacBio HiFi reads (15-20kb) can span a deletion, maintaining partial coverage via flanking homology. The depth drops but not to zero.
- **False DEL**: zero_fraction near 0% with high median_ratio = no deletion
- zero_fraction is most informative for deletions larger than read length

---

### Breakpoint Transition Metrics
| Metric | What it is | Example | Meaning |
|--------|-----------|---------|---------|
| `left_drop_ratio` | depth_after / depth_before at left breakpoint | 0.014 | After start, depth is 1.4% of before = **98.6% drop** |
| `right_drop_ratio` | depth_after / depth_before at right breakpoint | 1.075 | After end, depth is 107.5% of before = slight increase |
| `left_drop_sharpness` | Bases until depth crosses 50% threshold | 2 | Depth halves within 2 bp = **razor-sharp breakpoint** |
| `right_drop_sharpness` | Same for right breakpoint | 0 | Immediate transition |

**What to look for:**
- **drop_ratio < 0.3**: >70% depth drop = strong DEL signal at that breakpoint
- **drop_ratio < 0.1**: >90% drop = very strong
- **drop_ratio ~1.0**: no transition = no SV at that breakpoint
- **drop_ratio > 2.0**: depth spike = reads piling up at deletion boundary (common in long-read data)
- **sharpness < 50**: clean, biological breakpoint
- **sharpness > 200**: gradual transition = likely mapping artifact or repeat

**Why both breakpoints?** A true deletion should show a drop at the left breakpoint (entering the deletion) and recovery at the right breakpoint (exiting). If only one side drops, it may be a contig boundary or complex rearrangement.

---

### INV-Specific Metrics
| Metric | What it is | Example | Meaning |
|--------|-----------|---------|---------|
| `split_reads` | Total reads with SA tag in region | 48 | Number of split reads spanning the inversion |
| `split_reads_forward` | Split reads on forward strand | 8 | |
| `split_reads_reverse` | Split reads on reverse strand | 40 | |
| `strand_bias` | max(fwd_ratio, rev_ratio) | 0.833 | 83.3% of split reads on one strand = strong INV signal |

**What to look for:**
- **strand_bias > 0.65**: >65% of split reads favor one strand = supports inversion
- **strand_bias ~0.50**: 50/50 split = no strand preference = no INV signal
- **split_reads < 5**: insufficient evidence regardless of bias
- Depth should be ~normal (median_ratio ~1.0) since inversions are balanced events

---

## Interpreting a Real Case

### IMX2600 DEL 6,229 bp (LAR CONFIRMED by both assemblers)
Left flank: 372x
Right flank: 362x
SV median: 98x
p10: 89x
p90: 103x
zero_frac: 0.02%
median_ratio: 0.267
Left drop: 0.219 (78% drop, sharpness=0)
Right drop: 2.858 (spike, sharpness=0)

**Interpretation:**
- Flanks are symmetric (~370x) → clean region, no boundary issues
- SV median is 98x vs 370x flanks → **73% depth reduction** → real deletion
- zero_frac is 0.02% → reads still map inside the deletion. Why? PacBio reads are 15-20kb and span the 6.2kb deletion. They align using flanking homology, maintaining ~98x coverage.
- Left drop 0.219 → at the breakpoint, depth drops from ~370x to ~81x instantly (sharpness=0)
- Right drop 2.858 → depth spikes after the deletion (reads pile up at the boundary)
- p90 is only 103x → even the noisiest bases are well below flank depth

**Conclusion:** This is a real deletion. The depth drops 73-78% with sharp breakpoints. The residual coverage (98x) comes from long reads spanning the deletion — not from the deleted sequence being present. The previous threshold-based version incorrectly called this CONTRADICTS because median_ratio (0.267) was above the 0.05 cutoff. The numbers clearly tell the real story.

---

### Makgeolli CP025104.1 (LAR CONFIRMED but anomalous)
Left flank: 377x
Right flank: 405x
SV median: 95x
p90: 387x
zero_frac: 0.0%
median_ratio: 0.244
Left drop: ~1.0 (no drop)

**Interpretation:**
- zero_frac is 0.0% → every base has reads
- p90 is 387x → some areas have HIGHER depth than flanks
- No breakpoint drop → depth is continuous
- Despite median_ratio suggesting a drop, the p90 and zero_frac tell a different story

**Conclusion:** This is likely NOT a deletion. The 3 callers and LAR may have been misled by a complex rearrangement. Backtrack's numbers provide the counter-evidence.

---

### SX2 DEL (LAR CONTRADICTED — both assemblers found only 196-307bp)
Left flank: 127x
Right flank: 3x ← WARNING: near-zero right flank
SV median: 9x
p10: 0x
zero_frac: 44.4%
Left drop: 0.014 (98.6% drop, sharpness=2)
Right drop: 1.075

**Interpretation:**
- Right flank at 3x → **contig boundary or assembly gap** on the right side
- Left drop 0.014 → 98.6% depth drop at left breakpoint (very strong)
- zero_frac 44.4% → nearly half the region at zero depth
- SV median 9x → very low coverage inside the deletion

**Conclusion:** This looks like a real deletion at a contig boundary. LAR may have contradicted because the right side couldn't be assembled (no flanking sequence). The left-side signal is very strong. The scorer should weigh this as "likely real but at assembly boundary — treat with caution."

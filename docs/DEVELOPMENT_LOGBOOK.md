# FUNGUS-SV Development Logbook

## Project Overview
- **Pipeline:** FUNGUS-SV — Structural variant discovery and triangulation-based validation for haploid fungal genomes using PacBio HiFi reads.
- **Query strain:** *Saccharomyces cerevisiae* CICC-1445 (Chinese industrial)
- **Primary reference:** S288C (laboratory strain)
- **Additional references:** BJ4, IMX2600, Makgeolli, SX2
- **Data:** PacBio HiFi (SRR18210299), 274,915 reads, ~20 kb N50, ~100× coverage
- **GitHub:** https://github.com/keltonjenkovguimaraes-alt/fungus-sv
- **Lead:** Guimarães, K.H.A. et al. (2026, in preparation)

---

## Timeline

---

### 2026-05-27 — v0.8.0 Baseline

**Status at start:**
- ICB consensus: Sniffles2 + cuteSV + SVIM (≥2 callers, 0.5 overlap, 200 bp flank)
- 5-layer triangulation with original weights (Liu et al. 2024): depth=0.25, breakpoint=0.20, assembly=0.30, k-mer=0.25
- DHFFC thresholds from human diploid (Pedersen & Quinlan 2019): DEL < 0.7, DUP > 1.3
- 277 SVs called in CICC-1445 vs S288C
- LAR layer was a placeholder — never actually ran
- All 11 INVs scored CONTRADICTED (T=0.167) despite 100% split-read confirmation
- All 18 DUPs scored CONTRADICTED (T=0.167)

**Known issues identified:**
1. Weights not calibrated for haploid genomes
2. DHFFC thresholds too permissive (diploid values)
3. No DHBFC integration
4. No size-stratified scoring
5. LAR non-functional
6. INV scoring broken
7. DUP scoring broken
8. No gene context validation

---

### 2026-05-28 — Calibration Tasks 1-7

**Task 1: Weights recalibrated for haploid fungi**
- Depth: 0.25 → 0.35 (100% depth loss unambiguous in haploids)
- Breakpoint: 0.20 → 0.30 (only signal for inversions, critical for all)
- Assembly: 0.30 → 0.20 (computationally expensive, standalone tool)
- k-mer: 0.25 → 0.15 (redundant with depth in haploids)
- Source: Empirical data from CICC-1445 vs S288C

**Task 2: Inversion-specific scoring**
- INVs skip depth/k-mer layers
- Scored on breakpoint_junction + local_assembly only
- `INV_SPLIT_READ_CONFIRMED` tier for INVs with breakpoint ≥ 0.6
- All 11 S288C INVs confirmed by split reads (32-1,091 junction reads per breakpoint)

**Task 3: Distance-support formula implemented**
- `distance_support = 0.2 × sv_size + 2000 / sv_size`
- Source: Zheng & Shang (2024) PLOS ONE

**Task 4: MAPQ ≥ 20 filter**
- Applied via `samtools view -q 20`
- Source: Zheng & Shang (2024), Liu et al. (2024)

**Task 5: Coverage cap and haploid DHFFC thresholds**
- DEL strong: DHFFC < 0.3 (was < 0.7)
- DEL weak: 0.3 ≤ DHFFC < 0.7
- DUP strong: DHFFC > 2.0 (was > 1.3)
- DUP weak: 1.3 < DHFFC ≤ 2.0
- Source: rDNA DEL (DHFFC=0.03) and FLO1 (DHFFC=1.24) empirical data

**Task 6: BED files**
- S288C BED file exists
- Other strains (BJ4, IMX2600, Makgeolli, SX2) deferred — different chromosome naming (LR/CP vs NC_)

**Task 7: Genome-wide DHFFC validation**
- 266 DEL/DUP calls analyzed
- 115 confirmed (43.2%), 67 weak (25.2%), 84 contradicted (31.6%)
- FLO1, FLO9, HXT7 "deletions" in NOT_DEL/DUP category
- Ty2 deletions (YBL005W-B, YBR012W-B) confirmed with DHFFC 0.23-0.30

---

### 2026-05-29 — Multi-Strain LAR Calibration Begins

**LAR (Local Assembly Refinement) development:**
- Standalone tool at `valid_sv/evidence/layer_lar.py`
- Environment: `conda activate sv_lar` (Flye 2.9.6 + minimap2 + samtools)
- Regional assembly: extracts reads from SV region ± 3 kb, assembles with Flye, aligns contig to reference, parses CIGAR
- RAM: <500 MB per SV
- Time: 2-15 min per SV

**Initial LAR tests (4 SVs, 30 May 2026):**
| Locus | Type | Size | Strain | Result | Reads | RAM | Time |
|-------|------|------|--------|--------|-------|-----|------|
| YBL005W-B (Ty2) | DEL | 5.9 kb | S288C | ✅ CONFIRMED | 202 | <200 MB | 2 min |
| chrVII multi-gene | DEL | 55.7 kb | S288C | ⚠️ PARTIAL | 581 | <300 MB | 4 min |
| SX2 chrII | INV | 430 kb | SX2 | ✅ CONFIRMED | 6,448 | <500 MB | 13 min |
| BJ4 chrXII | INV | 205 kb | BJ4 | ✅ CONFIRMED | 3,552 | <400 MB | 10 min |

**Key finding:** 4/4 candidate SVs confirmed. LAR is viable as a standalone validation tool.

---

### 2026-06-01 — LAR Validation of Chromosome IV CONTRADICTED Calls

**Method:** Regional Flye assembly on 3 CONTRADICTED calls on chrIV.

**Results:**
| Test | Type | Size | Pipeline | LAR | Conclusion |
|------|------|------|----------|-----|------------|
| 1 | DEL | 6,269 bp | CONTRADICTED | ✅ REAL | Pipeline too conservative |
| 2 | DUP | 1,701 bp | CONTRADICTED | ⚠️ COMPLEX | Pipeline correct — not simple DUP |
| 3 | DEL | 120 bp | CONTRADICTED | ❌ FALSE | Pipeline correct |

**Key finding:** 2/3 CONTRADICTED calls correctly classified. One real 6.2 kb deletion missed.

**Multi-strain LAR calibration (64 CONTRADICTED DELs, 5.0-6.5 kb):**
| Strain | Tested | Real | False | % Real |
|--------|--------|------|-------|--------|
| S288C | 15 | 8 | 7 | 53% |
| BJ4 | 1 | 0 | 1 | 0% |
| IMX2600 | 28 | 5 | 23 | 18% |
| Makgeolli | 16 | 3 | 13 | 19% |
| SX2 | 4 | 1 | 3 | 25% |
| **Total** | **64** | **17** | **47** | **27%** |

**Key finding:** 27% of CONTRADICTED large DELs are real. DHFFC alone cannot predict reality across different reference genomes. S288C best reference quality (53% real); IMX2600 worst (engineered strain, 18% real).

---

### 2026-06-01 — Peter et al. (2018) Cross-Reference Analysis Plan

**Objective:** Validate FUNGUS-SV SVs against 1,011 yeast genomes population-scale data.

**5-Step Plan:**
| Step | Task | Status |
|------|------|--------|
| 1 | Build gene name mapping (4-letter codes → systematic ORF names) | ✅ Completed 2026-06-03 |
| 2 | Cross-reference affected genes with CN data | ✅ Completed 2026-06-03 |
| 3 | Validate Ty2/S288C-specific features | ⬜ Pending |
| 4 | Population context for CICC-1445 | ⬜ Pending |
| 5 | Sensitivity benchmark vs known CNVs | ⬜ Pending |

**Data downloaded:**
- `genesMatrix_CopyNumber.tab.gz` (980 KB, 7,796 ORFs × 1,011 isolates)
- `allORFs_pangenome.fasta.gz`
- `genesMatrix_PresenceAbsence.tab.gz`
- `1011DistanceMatrixBasedOnSNPs.tab.gz`
- `1011DistanceMatrixBasedOnORFs.tab.gz`

---

### 2026-06-02 — v0.9.3 Release

**New features:**
1. DHBFC integration (GC-corrected depth from Pedersen & Quinlan 2019)
2. Inversion-specific scoring with `INV_SPLIT_READ_CONFIRMED` tier
3. Duplication split-read rescue (+0.15 boost)
4. Repeat region flag `[REPEAT_REGION]` for FLO/rDNA/Ty
5. Translocation flag `[NEAR-ZERO_DEPTH]` for possible rearrangements
6. Small SV split-read minimum (≥6 reads for <100 bp)

**S288C v3 validation results:**
| Tier | Count | % |
|------|-------|---|
| TRIPLE_TRIANGULATED | 5 | 1.8% |
| DOUBLE_CONFIRMED | 17 | 6.1% |
| SINGLE_LINE | 1 | 0.4% |
| WEAK | 92 | 33.2% |
| CONTRADICTED | 151 | 54.5% |
| INV_SPLIT_READ_CONFIRMED | 11 | 4.0% |

**Known issues remaining:**
1. DHFFC reference-dependent (IMX2600 vs S288C)
2. 27% false negative rate in CONTRADICTED tier for large DELs
3. No BND/translocation caller
4. k-mer layer often unavailable (2.4 GB Jellyfish DB)
5. Spike-in calibration pending
6. Config weights not read by scorer code (old defaults hardcoded)

---

### 2026-06-03 — Code Audit: Weights and DHBFC Not Actually Implemented

**Discovery:** The v0.9.3 calibrated parameters exist in `config/config.yaml` but the Python code still uses original uncalibrated values.

**Bugs found:**
1. `scorer.py` DEFAULT_WEIGHTS: old values (0.30/0.25/0.25/0.20) — config ignored
2. `run_validation.py`: all LayerResult calls hardcoded to weight=0.25
3. `layer_depth.py`: `combined_ratio`, `dhbfc`, `size_factor` used but never computed
4. `layer_depth.py` `DepthEvidence` dataclass: missing `dhbfc` and `combined_ratio` fields
5. Ploidy weight: 0.0 in config, 0.15 in code, missing from scorer DEFAULT_WEIGHTS

**Fixes applied:**
- ✅ `scorer.py`: Weights updated to calibrated values (depth=0.35, breakpoint=0.30, LAR=0.20, k-mer=0.15, ploidy=0.0)
- ✅ `run_validation.py`: All LayerResult weights updated to match config
- ✅ `run_validation.py`: `import yaml` added, config weights loaded
- ✅ `layer_depth.py`: DHBFC computation added (10 kb local context)
- ✅ `layer_depth.py`: Size-stratified scoring factor added
- ✅ `layer_depth.py`: Repeat region and translocation warnings added
- ✅ `layer_depth.py`: `DepthEvidence` dataclass updated with `dhbfc` and `combined_ratio` fields
- ✅ `layer_depth.py`: All return statements updated with new fields

**Post-fix test (29 kb DUP, LAR-confirmed REAL):**
| Metric | Before Fix | After Fix |
|--------|-----------|-----------|
| DHFFC | 0.8489 | 0.8489 |
| DHBFC | N/A | 0.9082 |
| Combined ratio | N/A | 0.8786 |
| T-score | 0.1667 | 0.5538 |
| Tier | CONTRADICTED | SINGLE_LINE |

**DUPs still under-scored** — depth layer returns 0.3 (ambiguous) even with DHBFC. The depth signal doesn't reach >2.0 for these real DUPs. Rescue boost helps but doesn't reach DOUBLE_CONFIRMED.

---

### 2026-06-03 — LAR on 5 Largest S288C DUPs

**Results:**
| # | SV ID | Size | Reads | Contigs | Verdict |
|---|-------|------|-------|---------|---------|
| 1 | ICB_NC_001146.8_562192_221 | 40,203 bp | 246 | 3 | ✅ CONFIRMED |
| 2 | ICB_NC_001224.1_4153_256 | 38,812 bp | 1,013 | 1 | ✅ CONFIRMED |
| 3 | ICB_NC_001134.8_197384_269 | 29,549 bp | 193 | 3 | ✅ CONFIRMED |
| 4 | ICB_NC_001146.8_547075_220 | 21,163 bp | 168 | 3 | ✅ CONFIRMED |
| 5 | ICB_NC_001144.5_451418_275 | 17,511 bp | 0 | 0 | ❌ TIMEOUT |

**Key finding:** 4/5 largest DUPs are real despite all scoring T=0.167 CONTRADICTED (estimated FDR 1.0) by the pipeline. The depth layer systematically fails DUPs.

**Samplot images generated** for all 4 confirmed DUPs in `figures/samplot_DUPs/`.

---

### 2026-06-03 — Layer 6: Genomic Context Filters

**New file:** `valid_sv/evidence/layer_genomic_context.py`

**Design:**
- SV-type-specific hard filters (PASS/FLAG/FAIL)
- Does NOT contribute to T-score
- Positioned after triangulation, before ploidy
- Uses pre-computed SV annotation TSVs

**Filter rules:**
| SV Type | Filters | Source |
|---------|---------|--------|
| DUP | Gene continuity, ORF integrity, S288C-specific check, mitochondrial exception | FocalSV-adapted + this study |
| DEL | Size-stratified reliability (<500 bp = FLAG), essential gene check, repeat region flag | Pedersen & Quinlan 2019 + this study |
| INV | Gene disruption at breakpoints | This study |
| TRA/BND | Gene fusion detection | This study |

**Essential genes catalog:** BDP1, GPI15, COG6, SSN8, SAM50, HHT2, HHF2, RRN6, FMT1, SCT1, HIR1, SLA1, PDR3, COX1, ATP8, ATP6, COB, IDH1, NCE103, BOP3

**S288C-specific features:** YBL005W-B (Ty2), YBR012W-B (Ty2), YBL005W-A, YBR012W-A

**DEL test results:**
| Test | SV | Size | LAR Truth | Layer 6 Verdict |
|------|-----|------|-----------|-----------------|
| 1 | chrIV DEL | 6,269 bp | ✅ REAL | PASS |
| 2 | chrIV DEL | 120 bp | ❌ FALSE | FLAG (small DEL) |
| 3 | Ty2 DEL | 5,921 bp | ✅ REAL | PASS + REPEAT FLAG |
| 4 | FLO9 DEL | 516 bp | Complex | PASS + REPEAT FLAG |

**DUP test results:**
| DUP | Size | LAR Truth | Layer 6 Verdict |
|-----|------|-----------|-----------------|
| 1 | 40 kb chrVIII | ✅ REAL | FLAG (essential genes: HHT2, HHF2, SAM50, SSN8) |
| 2 | 38 kb chrM | ✅ REAL | FLAG (mitochondrial — multicopy exception) |
| 3 | 29 kb chrII | ⚠️ COMPLEX | **FAIL** (S288C-specific Ty2 insertion) |
| 4 | 21 kb chrVIII | ✅ REAL | FLAG (essential genes: BDP1, GPI15, COG6, IDH1) |

**Key finding:** Layer 6 correctly catches the S288C-specific Ty2 DUP as FAIL and flags small DELs. Essential gene warnings are informational (FLAG), not blocking.

---

### 2026-06-03 — Peter et al. (2018) Step 1 & 2 Completed

**Step 1: Gene name mapping built**
- BLAST: 7,348 pangenome ORFs → S288C chromosomes
- GFF intersection: 7,161 Peter ORFs mapped to S288C gene names (97.5% success)
- 692 out of 1,011 Peter 4-letter codes mapped to S288C systematic names
- Mapping saved to `data/peter2018_mapping/peter_code_to_s288c.json`

**Step 2: Cross-reference with CN data**
- 66 FUNGUS-SV affected genes cross-referenced with 1,011 strains
- **Critical finding:** The `genesMatrix_CopyNumber.tab.gz` contains pangenome presence/absence, not quantitative copy number
- All genes show ~22% CN=0 (mean 22.1%, min 21.8%, max 22.6%) — impossibly narrow distribution
- CN=0 means "gene family not detected in that strain's assembly" rather than "gene biologically deleted"
- Cannot be used for quantitative population CNV validation

**Decision:** Peter et al. data limited to presence/absence validation. Alternative needed for true CNV benchmarking (SGD, spike-in, or new dataset).

---

### 2026-06-03 — FocalSV (Luo et al. 2025) Paper Review

**Key relevance to FUNGUS-SV:**
- FocalSV performs target region assembly-based SV detection — the publication-validated version of LAR
- Diploid-specific features: haplotype phasing, genotype refinement, dual assembly — NOT applicable
- Genome-agnostic features applicable to FUNGUS-SV:
  - DUP via insertion realignment (take INS allele, realign — if maps near breakpoint → DUP)
  - Split-read DUP signatures (FocalSV Eq. 7)
  - Split-read INV signatures (FocalSV Eq. 15)
  - Split-read TRA detection (FocalSV Eq. 11)
  - Intra + inter-alignment DEL merging
  - Signature clustering (100-500 bp thresholds)

**LAR DUP logic critique:** Current `total_ins >= sv_size * 0.5 or alignment_count >= 2` is too simplistic. FocalSV's insertion realignment method is more accurate.

**FocalSV somatic DUP F1 scores:** 63-70% (PacBio), 55-69% (ONT) — DUPs are hard even for dedicated tools.

---

### 2026-06-03 — Samplot (Belyeu et al. 2021) Paper Review

**Key relevance:**
- Samplot displays 3 SV evidence categories: split reads, discordant pairs, coverage depth
- No computational decision-making — leaves interpretation to human
- Samplot-ML CNN reduces false positives by 51.4% (short-read), 27.8% (long-read)
- **Limitation:** Samplot-ML trained only on deletions — not applicable to DUPs, INVs
- True-negative strategy: sample from "exclude regions" (problematic genomic regions)
- Visual curation: 91% false positive rate before review, 7% after — human eye still beats automated filters

---

### 2026-06-03 — Li (2014) Artifacts Paper Review

**Key findings applicable to FUNGUS-SV:**

1. **Haploid = built-in FDR estimator:** Li used CHM1 haploid human cell line — heterozygous calls = errors. Same logic as FUNGUS-SV ploidy layer (het < 7%).

2. **Low-complexity regions (LCRs) = 80-90% of INDEL errors:** 2% of genome harbors majority of false calls. FUNGUS-SV's `[REPEAT_REGION]` flag should be expanded to include homopolymer runs and DUST/mdust-identified LCRs.

3. **Max depth filter most effective:** `depth > d + 3√d` was the single best filter. For ~100× coverage: threshold = 130×. FUNGUS-SV currently uses `max_coverage_multiple: 5` — less precise.

4. **Incomplete reference causes false variants:** When reference misses sequences present in sample, reads mismap. Explains IMX2600's low LAR confirmation rate (18%) — engineered strain with assembly differences.

5. **Local assembly beats realignment:** Li's fermi assembler outperformed all variant callers. Supports LAR as gold standard.

6. **PCR duplicates = INDEL artifacts:** PacBio HiFi is PCR-free — less concern for FUNGUS-SV.

---

### 2026-06-03 — Multi-Assembler LAR Concept

**Rationale:** Assembly consensus (Flye + miniasm + hifiasm) is truly orthogonal — different algorithms, different error modes. Caller consensus shares BAM input (correlated errors).

**Tiered assembly strategy:**
| Tier | Assembler | RAM | Time | Use Case |
|------|-----------|-----|------|----------|
| 1 | Flye | <500 MB | 2-15 min | All SVs (default) |
| 2 | Miniasm | <200 MB | <1 min | Contradicted or partial Flye results |
| 3 | hifiasm | ~5 GB | 5-20 min | High-priority unresolved calls (NOT viable on current 5.7 GB system) |

**Hardware constraint:** System has 5.7 GB RAM. hifiasm OOM-killed at 5.5 GB even on 1,561 reads. Multi-assembler consensus limited to Flye + miniasm on current hardware.

---

### 2026-06-03 — Layer 6 Architecture Decision

**Position in pipeline:** After LAR, before ploidy.

**Rationale:** LAR is the strongest single piece of evidence. If LAR says CONFIRMED, Layer 6 flags become informational, not gatekeeping. If LAR not run, Layer 6 identifies suspicious calls for LAR prioritization.

**Flow:**
Triangulation (T-score) → LAR (manual) → Layer 6 (PASS/FLAG/FAIL) → Ploidy → Tier

---

### 2026-06-04 — SX2 chrV 80.7 kb INV: Flye vs Miniasm

**SV:** ICB_LR813589.2_362102_92, INV, 80,695 bp, 3 callers (cutesv + sniffles2 + svim)

**Flye result:**
- Reads: 781, contigs: 1
- Verdict: CONTRADICTED — "no strand change detected"

**Miniasm result (NEW — second assembler):**
- Reads: 1,561, contigs: 1 (126,972 bp)
- RAM: 0.092 GB (92 MB)
- Time: 0.4 sec
- Verdict: PARTIAL — 4 alignments on LR813589.2, all reverse strand, positions spanning 337,096-459,395

**Two-assembler consensus:**
| Assembler | Verdict |
|-----------|---------|
| Flye | ❌ CONTRADICTED |
| Miniasm | ⚠️ PARTIAL (complex rearrangement) |

**Interpretation:** NOT a clean simple inversion. Complex rearrangement on chrV — possibly inversion with duplications, or translocated chrV sequence. The 3-caller consensus detected real structural complexity, but it's more complex than a single INV.

**Multi-assembler LAR proven viable:** Miniasm runs at 92 MB RAM / 0.4 sec — feasible as Tier 2 assembler on 5.7 GB system.

---

### 2026-06-04 — Flye Parameter Audit

**Current command (layer_lar.py line 104-107):**
flye --pacbio-hifi <fastq> --genome-size <window_size> --read-error 0.005 --meta --threads <threads> --out-dir <dir>

**Status:** Not yet applied. Pending re-run of SX2 chrV INV with improved parameters.

---

## Environment Structure

| Environment | Tools | Purpose |
|-------------|-------|---------|
| `sv_align` | minimap2, samtools | Read alignment |
| `sv_call` | sniffles2, cutesv, svim, bcftools | SV detection |
| `sv_valid` | python, numpy, scipy, pandas, pysam, matplotlib | Triangulation |
| `sv_lar` | flye, minimap2, samtools, miniasm | LAR assembly |
| `sv_samplot` | samplot | Visualization |
| `assembly` | hifiasm | High-accuracy assembly (OOM on 5.7 GB) |

---

## Files Modified (v0.9.3 → v0.9.4-dev)

| File | Change | Date |
|------|--------|------|
| `valid_sv/engine/scorer.py` | Calibrated DEFAULT_WEIGHTS (depth=0.35, breakpoint=0.30, LAR=0.20, k-mer=0.15, ploidy=0.0) | 2026-06-03 |
| `valid_sv/run_validation.py` | Config-loaded weights, yaml import, all LayerResult weights updated | 2026-06-03 |
| `valid_sv/evidence/layer_depth.py` | DHBFC, size_stratification, repeat/translocation flags, updated dataclass | 2026-06-03 |
| `valid_sv/evidence/layer_genomic_context.py` | **NEW** — SV-type-specific gene annotation filters | 2026-06-03 |
| `config/config.yaml` | Calibrated weights (already correct, now code reads it) | 2026-05-28 |
| `docs/DEVELOPMENT_LOGBOOK.md` | **NEW** — This file | 2026-06-04 |

---

## Next Steps (Priority Order)

1. ⬜ Apply improved Flye parameters (`--read-error 0.005 --meta`) and re-run SX2 chrV INV
2. ⬜ Run LAR on remaining 13 S288C DUPs for complete DUP calibration
3. ⬜ Wire `layer_genomic_context.py` into `run_validation.py` as post-scoring hard filter
4. ⬜ LAR-test ~10 random HIGH/DOUBLE_CONFIRMED SVs (measure false positive rate at top)
5. ⬜ LAR-test ~10 random WEAK SVs (characterize middle tier)
6. ⬜ Implement max depth filter: `depth > mean_depth + 3*sqrt(mean_depth)` (Li 2014)
7. ⬜ Add mdust/DSUST low-complexity region detection
8. ⬜ Expand `[REPEAT_REGION]` to include homopolymer runs
9. ⬜ Peter et al. Step 3: Validate Ty2/S288C-specific features against pangenome
10. ⬜ Build fungal HiFi read simulator for spike-in benchmarking
11. ⬜ Add BND/translocation detection (FocalSV Eq. 11)
12. ⬜ Integrate FocalSV DUP insertion realignment into LAR CIGAR logic
13. ⬜ Add miniasm as Tier 2 assembler in `layer_lar_tiered.py`
14. ⬜ Re-run all 5 strains with v3 fixes (only S288C done)
15. ⬜ Publish to Zenodo

---

## Papers Integrated

1. Pedersen & Quinlan (2019) — Duphold: DHFFC/DHBFC metrics
2. Zheng & Shang (2024) — SVvalidation: distance_support, MAPQ filter
3. Belyeu et al. (2021) — Samplot: visualization, ML curation
4. David et al. (2024) — Manual curation FDR estimates
5. Dhakal et al. (2024) — Fungal SV landscape, repeat regions
6. Li et al. (2023) — Haploid fungal variant calling
7. Nkouamedjo et al. (2025) — SV-MeCa: XGBoost meta-caller
8. Sedlazeck et al. (2017) — SURVIVOR: SV annotation
9. Liu et al. (2024) — ICB consensus, triangulation weights
10. Zhang et al. (2025) — SMaHT benchmark, confidence tiers
11. Luo et al. (2025) — FocalSV: target region assembly-based SV detection
12. Li (2014) — Artifacts in variant calling from high-coverage samples
13. Peter et al. (2018) — 1011 yeast genomes
14. Xing et al. (2025) — Haploid fungi heterozygosity thresholds

---

*Last updated: 2026-06-04*


---

### 2026-06-04 (Session 2) — All Remaining Implementations Complete

**Queue completed:**

| # | Implementation | File | Status |
|---|---------------|------|--------|
| 1 | Flye `--read-error 0.005 --meta` | layer_lar.py | ✅ |
| 2 | Assembly ploidy check | layer_lar.py | ✅ |
| 3 | Miniasm+Racon pipeline (`run_lar_miniasm`) | layer_lar.py | ✅ |
| 4 | INV repeat escalation (Cheng & Sedlazeck 2025) | layer_genomic_context.py | ✅ |
| 5 | Size-stratified INV breakpoint precision | layer_genomic_context.py | ✅ |
| 6 | Max depth filter (Li 2014): `depth > mean + 3√mean` | layer_depth.py | ✅ |
| 7 | Low-complexity region flag (Li 2014) | layer_depth.py | ✅ |
| 8 | Layer 6 wired into run_validation.py | run_validation.py | ✅ |
| 9 | `samtools fastq -o` bug fixed (stdout redirect) | layer_lar.py | ✅ |

**SX2 chrV 80.7 kb INV — Final Three-Assembler Consensus:**

| Assembler | Verdict | Evidence |
|-----------|---------|----------|
| Flye (--meta) | CONTRADICTED | No strand change |
| Miniasm (raw) | PARTIAL | Multiple reverse-strand alignments |
| **Miniasm+Racon** | **CONFIRMED** | Opposite strand + ploidy=1.49 |

**Conclusion:** Real complex rearrangement (inversion + partial duplication).
3-caller consensus correct. LAR with multiple assemblers resolves ambiguity.

**All papers integrated:**
1. Pedersen & Quinlan (2019) — DHFFC/DHBFC
2. Zheng & Shang (2024) — distance_support, MAPQ
3. Belyeu et al. (2021) — Samplot visualization
4. David et al. (2024) — Manual curation FDR
5. Dhakal et al. (2024) — Fungal SV landscape
6. Li et al. (2023) — Haploid fungal benchmark
7. Liu et al. (2024) — ICB consensus, weights
8. Zhang et al. (2025) — SMaHT tiers
9. Luo et al. (2025) — FocalSV
10. Li (2014) — Artifacts, max depth, LCR
11. Peter et al. (2018) — 1011 yeast genomes
12. Mochizuki et al. (2023) — Assembly guideline, ploidy, polishing
13. Cheng & Sedlazeck (2025) — Inversion benchmark

**Ready for next session:**
- Run LAR on remaining 13 S288C DUPs
- LAR-test HIGH/DOUBLE_CONFIRMED SVs
- Test full pipeline with Layer 6 integrated
- Peter et al. Step 3 (Ty2 validation)


---

### 2026-06-05 — Complete DUP LAR Calibration Dataset

**S288C DUPs (18/18 complete):**
- 16/18 confirmed (88.9%)
- 1 timeout (chrXI 17.5 kb)
- 1 insufficient reads (chrI 913 bp, only 12 reads)
- All 18 scored CONTRADICTED (T=0.167) by pipeline

**SX2 DUPs (18/18 complete):**
- 15/18 confirmed (83.3%)
- 1 mtDNA assembly failed (chrM 42 kb)
- 1 timeout (chrXIV 3.2 kb)
- 1 timeout (chrI 61 bp)
- All 18 scored CONTRADICTED by pipeline

**Cross-strain DUP truth set:**
| Strain | Tested | Confirmed | % Real |
|--------|--------|-----------|--------|
| S288C | 18 | 16 | 88.9% |
| SX2 | 18 | 15 | 83.3% |
| **Total** | **36** | **31** | **86.1%** |

**Remaining DUPs:** BJ4 (13), IMX2600 (20), Makgeolli (17) = 50 remaining

**Key finding:** Pipeline systematically fails DUPs across strains.
86% of CONTRADICTED DUPs are real. Depth layer (DHFFC >2.0) 
cannot confirm duplications in haploids. LAR is essential for DUP validation.

**LAR performance metrics:**
- Average reads extracted: 150-300 per DUP
- Flye RAM: <500 MB, Time: 2-10 min per DUP
- Success rate: 86% (31/36 confirmed)
- Failure modes: timeout (3), insufficient reads (1), mtDNA assembly (1)

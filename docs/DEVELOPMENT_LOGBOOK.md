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

---

## 2026-06-16 to 2026-06-18 — v1.0.0 Publication Readiness Sprint

### Docker Container
- Created Dockerfile with all 5 conda environments (sv_align, sv_call, sv_valid, sv_lar, sv_kmers)
- Added .dockerignore to exclude 55 GB of data from build context
- Fixed lar.yaml: added miniasm=0.3 and racon=1.5.0 (present locally but missing from YAML)
- Fixed validation.yaml: added samtools and bcftools (present locally via system PATH, missing in container)
- Fixed kmers.yaml: unpinned jellyfish version (2.3.0 not available in Docker architecture)
- Added int() casts in run_validation.py for pos/end parameters (Python 3.11 type safety)
- Added HEALTHCHECK to Dockerfile
- Built and tested Docker image — produces identical T-scores to local installation
- Docker image size: 4.3 GB (all 5 environments with tools)

### Benchmarking
- Ran MUMmer4 on all 5 Saccharomyces and 5 Candida strains
- Ran SVIM-asm on all 10 strains
- Compared individual callers (Sniffles2, cuteSV, SVIM) vs ICB consensus
- Key findings:
  - ICB consensus removes 80–95% of single-caller noise
  - SVIM-asm missed 6/7 LAR-confirmed SVs on S288C
  - SVIM-asm found 0 DUPs across all 10 strains
  - FUNGUS-SV found 9–70 DUPs per strain
  - Assembly fragmentation inflates SV counts (UAB012: 893, ATCC64124: 111)
- Created scripts/run_benchmarking.sh for reproducibility

### LAR Truth Set Expansion
- Ran 10 additional LAR SVs (5 Saccharomyces + 5 Candida) — batch incomplete
- 2 completed: S288C INV 168 kb (PARTIAL), BJ4 INV 65 kb (PARTIAL)
- 3 timeouts resolved by Samplot: all confirmed real
  - S288C INV 168 kb: Real INV (Samplot zoom confirmed)
  - SX2 DUP 59 kb: Real, complex (repetitive/multicopy context)
  - 101 INV 508 kb: Real, complex (INV with internal deletion)
- Total LAR truth set: 140 SVs across 12 strains

### Samplot Visual Validation
- Completed SC5314 self-alignment Samplot review (10/10 SVs)
  - 8 confirmed real, 1 false positive, 1 real but mis-sized
  - LAR + Samplot accuracy: 10/10 (100%)
- Key findings:
  - --zoom flag essential for SVs >50 kb
  - Two-assembler design prevented false negatives (Flye missed 2 INVs, Miniasm caught them)
  - Flye alone produced 1 false positive (Miniasm contradicted, Samplot confirmed)
  - ~32% of SVs require visual tie-breaking

### Snakefile Fixes
- Fixed callers: "pbsv" → reads from config.yaml (sniffles2, cutesv, svim)
- Fixed conda paths: absolute (/home/kelto/...) → relative YAML (workflow/envs/*.yaml)
- Commented out dead lar_refine rule (local_assembly.py was deleted)
- 6 active rules remain

### Test Dataset
- Created minimal test dataset: S288C mitochondrial chromosome (NC_001224.1, 85 kb, 8,789 reads)
- Full pipeline runs in <5 minutes
- Verified end-to-end in both local and Docker

### Documentation
- BENCHMARKING_COMPLETE.md — full comparison across all tools
- SAMPLOT_VALIDATION_COMPLETE.md — 19-SV cross-genus visual validation
- PRESENTATION_FUNGUS_SV.md — presentation guide with 10 schemes
- SCHEMES_PIPELINE.md — detailed pipeline diagrams
- PUBLICATION_SUMMARY.md — final summary of everything built
- README.md — comprehensive v1.0.0 README

### GitHub
- Pushed all core code changes (Snakefile, env YAMLs, Dockerfile, .dockerignore)
- Clean commit history
- Local and GitHub in sync

### Storage
- Docker data disk (docker_data.vhdx) compacted from 88 GB → 2.22 GB
- Reclaimed ~86 GB on C: drive

### Status
- Pipeline: PUBLICATION READY ✅
- Docker: TESTED AND VERIFIED ✅
- Benchmarking: COMPLETE ✅
- Documentation: COMPLETE ✅
- Manuscript: PENDING


---

## 2026-06-18 to 2026-06-20 — v1.0.0 Final Sprint: Figures, External Validation, Manuscript

### Docker Pipeline Completion
- Built and tested Docker container with all 5 conda environments
- Discovered samtools missing from validation.yaml — added samtools and bcftools
- Discovered miniasm and racon missing from lar.yaml — added both
- Added .dockerignore to prevent 55 GB of data from being copied into image
- Docker image size: ~4.3 GB
- End-to-end Docker test produces identical T-scores to local installation
- Added HEALTHCHECK to Dockerfile
- Fixed jellyfish version pinning (unpinned for Docker compatibility)
- Fixed int() casts in run_validation.py for Python 3.11 compatibility

### Storage Crisis
- Docker build cache consumed 88 GB on C: drive
- Compressed docker_data.vhdx from 88 GB → 2.22 GB using diskpart
- Freed ~86 GB on Windows C: drive
- Moved Docker storage recommendation to D: drive

### Publication Figures Created (8 main + supplementary)

**Main Figures:**
1. `figure1_pipeline_workflow.png` — End-to-end pipeline schematic
2. `figure_sv_counts.png` — SV counts per strain (both genera, stacked bars)
3. `figure_icb_noise.png` — ICB consensus noise reduction across all strains
4. `figure_lar_truth.png` — LAR truth set results (145 SVs, 13 strains)
5. `figure_self_alignment.png` — Self-alignment empirical FDR baseline
6. `figure_benchmarking.png` — Tool comparison (FUNGUS-SV vs 6 other tools)
7. `figure_assembly_quality.png` — Assembly quality impact on SV detection
8. `figure_false_positive_categories.png` — Why LAR rejects SVs (38 rejected)
9. `figure_ogg_comparison.png` — FUNGUS-SV vs Oggenfuss et al. (2025)

**Scripts:**
- All figures generated via matplotlib scripts in `figures/plot_*.py`
- Samplot annotation scripts created for publication-quality images

### Math Audit — Complete Verification
- All percentages recalculated from raw data
- ICB noise reduction: S288C=92.7%, UAB012=85.4%, L26=90.1%, P75063=90.3%
- LAR within-species: Saccharomyces=72% (36/50), Candida=74% (37/50), Combined=73% (73/100)
- DUP calibration: 39/48=81.2% confirmed
- Self-alignment TRIPLE: CICC-1445=3.1% (4/130), SC5314=2.5% (5/200)
- SVIM-asm LAR recovery: 1/7=14.3%
- Assembly quality: UAB012 TRIPLE=46.5% (416/893), ATCC64124 TRIPLE=9.0% (10/111)
- Total SVs: Saccharomyces=1,527, Candida=2,083, Grand=3,610

### External Tool Comparison Expanded
- Added DeBreak to benchmarking: 521 SVs on S288C (232 DEL, 16 DUP, 22 INV, 221 INS, 30 TRA)
- DeBreak detects insertions and translocations that FUNGUS-SV doesn't call
- Ran Samplot on 5 largest DeBreak SVs: 2 confirmed false, 3 pending
- Ran LAR on DeBreak calls: DEL 514kb contradicted (44kb actual), DUP 1Mb timed out
- DeBreak's 1Mb DUP shown by Samplot to be a deletion with progressive coverage drop
- Total comparison tools: 7 (Sniffles2, cuteSV, SVIM, DeBreak, FUNGUS-SV ICB, SVIM-asm, MUMmer4)
- Attempted NanoSV and SVision — NanoSV failed (BED requirement), SVision not attempted

### LAR Truth Set Expansion
- Batch 3: 10 largest untested SVs across both genera
  - 3 CONFIRMED (BJ4 INV 44kb, WO1 INV 35kb, UAB012 DUP 69kb)
  - 3 PARTIAL (S288C INV 139kb, FDAARGOS656 INV 70kb, Makgeolli DEL 13kb)
  - 4 CONTRADICTED (101 INV 85kb, S288C DEL 55kb, IMX2600 DUP 49kb, SX2 DEL 19kb)
- Total LAR-validated: 145 SVs (across all batches)
- LAR rejection categories: Size Mismatch (12), Timeout (9), Other (7), Assembly Failed (6), Contradicted (4)

### External Validation — Oggenfuss et al. (2025)
- Downloaded 3 C. albicans genome assemblies from NCBI (PRJNA967712)
- Strains: SC5314, L26, P75063 — originally analyzed with ONT + DELLY
- Ran FUNGUS-SV (PacBio HiFi) on SC5314 vs L26 and SC5314 vs P75063
- Results:
  - L26: 230 SVs (FUNGUS-SV) vs 679 (DELLY) — 90.1% ICB reduction
  - P75063: 333 SVs (FUNGUS-SV) vs 864 (DELLY) — 90.3% ICB reduction
- Identified centromeric inversions matching Oggenfuss findings:
  - CEN4 (JBIBQQ010000014.1, 3.1 kb): LAR CONFIRMED (both assemblers)
  - CEN5 (JBIBQQ010000005.1, 20.8 kb): LAR PARTIAL (Flye confirmed, Miniasm contradicted)
- First external validation of FUNGUS-SV against published data

### Critical Self-Review
- Identified overconfident claims needing revision:
  1. "72% confirmation rate" → Add technically-successful rate (88%)
  2. "FDR ~0% for DUP/INV" → "No biologically contradictory cases observed"
  3. "100% Samplot accuracy" → "All 19 inspected consistent with consensus"
  4. Weight selection → Acknowledge heuristic, not formal optimization
  5. Layer independence → Acknowledge partial dependence
  6. TRIPLE FDR ~3% → Note wide confidence intervals with 4-5 events
- All changes will be reflected in final manuscript

### Manuscript Status
- Complete draft with all sections
- Figure legends written for all 9 main figures
- Supplementary figure plan defined (7 figures)
- [FILL] sections identified for author input
- External validation section strengthened with Oggenfuss comparison
- Honest limitations section drafted

### GitHub
- All core code committed and pushed
- Dockerfile, .dockerignore, updated env YAMLs pushed
- Clean commit history

### Next Steps
- Final manuscript revisions
- Zenodo upload
- Submit to BMC Bioinformatics or Genome Biology


2026-06-26 — Samplot Review & Backtrack Layer Development

Precision-Recall Analysis
- Ran precision-recall on all 133 LAR-matched SVs (including batch4: 20 new SVs from June 25)
- 75 CONFIRMED, 25 BIOLOGICAL_FP, 33 TECHNICAL_FAILURE
- Excluding technical failures: 100 biologically resolved SVs
- DEL tier system works: TRIPLE (T≥0.80) = 87.5% precision (14/16), DOUBLE = 88.2% (15/17)
- DUP/INV tiers inverted by raw T-score: all INVs T<0.20, all DUPs initially CONTRADICTED
- CONTRADICTED-but-CONFIRMED: 36 SVs (21 DEL, 15 INV) — overrides fix this
- Generated figure_precision_recall.png (4 panels: PR curve, T-score distribution, 3-category bar, precision by SV type)

TRIPLE-Tier False Positive Investigation
- Identified 19 unique TRIPLE (T≥0.80) non-CONFIRMED SVs across all LAR batches
- Classified into 3 categories:
  - Category A (Technical Failures): 4 SVs — timeout/insufficient reads → excluded from analysis
  - Category B (PARTIAL): 6 SVs — one assembler confirmed, one failed/contradicted → Samplot needed
  - Category C (Both contradicted): 4 SVs — both assemblers agree SV is false → likely true FPs
  - Category D (Size mismatch): 2 SVs — real SV but pipeline overestimated size → reclassify as CONFIRMED
- Generated 12 Samplots for Category B, C, D SVs
- Fixed contig name mapping: CP127* = IMX2600, CP025* = Makgeolli, LR813523.2 = BJ4, CM000309.1 = WO1
- Samplots saved to figures/TRIPLE_FP_samplot/

Batch 4 LAR (20 SVs, June 25)
- Extracted and classified all 20 SVs from data/LAR_batch4/LAR_batch4.log
- 13 CONFIRMED (both assemblers), 1 PARTIAL, 2 CONTRADICTED, 4 technical failures
- Key finding: DEL sizes systematically underestimated by pipeline (called 6.5-14.7 kb, LAR finds 12.5-60 kb)
- DUP confirmation: 4/4 technically successful DUPs confirmed (100%)
- INV confirmation: 3/3 confirmed (100%)
- PLOIDY_WARNING on most Miniasm runs (assembly_ploidy 1.41-12.40)
- Integrated into precision-recall dataset

Backtrack Layer Development (layer_backtrack.py)
- New evidence layer: aligns reads to original vs SV-modified reference sequences
- Fast, orthogonal to depth/breakpoint/LAR, works for all SV types
- Initial bug: capture_output=True failed on binary BAM output from samtools
- Fix 1: Changed to stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
- Second bug: minimap2 not installed in sv_valid environment
- Fix 2: Hardcoded path to sv_align minimap2 (/home/kelto/miniforge3/envs/sv_align/bin/minimap2)
- Third issue: Read count metric couldn't distinguish large DELs (609 vs 609 aligned reads)
- Fix 3: Changed metric to split-read count (SA tag) for DELs — more sensitive to breakpoint-spanning reads
- Fix 4: Inverted ratio logic for DELs (ratio = orig_splits / mod_splits; high ratio = supports deletion)
- Fix 5: Adjusted verdict thresholds for DELs: ratio >5.0 = STRONG_SUPPORT, >2.0 = WEAK, <1.1 = CONTRADICTS

Backtrack Validation Results
- Makgeolli DEL 7kb (CP025104.1_471195_112, LAR CONFIRMED both):
  - 107 split-reads in original, 6 in modified → ratio 17.83 → STRONG_SUPPORT (score 1.000)
  - Quadruple-confirmed: ICB + Triangulation + LAR + Backtrack
- Makgeolli DEL 19.3kb (CP025112.1_916639_231, LAR PARTIAL):
  - 480 split-reads in original, 33 in modified → ratio 14.55 → STRONG_SUPPORT
- IMX2600 DEL 12.5kb (CP127210.1_446292_298, LAR CONFIRMED both):
  - 43 split-reads in original, 24 in modified → ratio 1.79 → AMBIGUOUS
  - Tested multiple flank sizes (5kb, 10kb, 15kb) — same result

Makgeolli CP025104.1 Anomaly
- Samplot review shows DUP-like signal (high insert size), not DEL
- Depth shows 327.5x in SV region, zero zero-depth bases — inconsistent with haploid DEL
- All 3 callers + LAR + backtrack say DEL → contradiction with raw coverage/insert size
- Possible complex SV (inverted duplication?) — requires further investigation

Files Modified
- valid_sv/evidence/layer_backtrack.py — major revisions (minimap2 path, split-read metric, DEL-specific logic)
- docs/DEVELOPMENT_LOGBOOK.md — this update
- figures/figure_precision_recall.png — generated
- figures/TRIPLE_FP_samplot/ — 12 new Samplots

Pending
- Samplot review of 12 TRIPLE-tier SVs
- Backtrack calibration against full LAR truth set
- Integration of backtrack layer into run_validation.py scoring
- Resolution of Makgeolli CP025104.1 anomaly
- Final precision-recall figure after Samplot reclassification


2026-07-01 — Backtrack Layer Rewrite: Pure Depth Reporting

Design Decision
- Removed all verdict logic and hard thresholds from backtrack layer
- Backtrack now reports raw depth metrics only — no STRONG_SUPPORT/WEAK_SUPPORT/CONTRADICTS
- Scoring and tier assignment belong in the triangulation scorer, not in individual evidence layers
- This makes backtrack truly orthogonal: it measures depth facts, not interpretations

Architecture Change
- Before: 7 samtools depth calls (3 regions + 4 breakpoint windows), then verdict logic
- After: 3 samtools depth calls (left flank, SV region, right flank), breakpoints use sliced pre-fetched arrays
- Runtime reduced by ~60% (7→3 external calls)

Metrics Reported (BacktrackReport dataclass)
- Flanks: left_flank_mean, left_flank_median, right_flank_mean, right_flank_median
- SV region: sv_region_mean, sv_region_median, sv_region_p10, sv_region_p90
- Ratios: mean_ratio (sv_mean/flank_mean), median_ratio (sv_median/flank_median)
- Sparsity: zero_fraction (bases with depth<5), low_fraction (bases with depth<10)
- Breakpoints: left_drop_ratio (depth_after/depth_before at start), right_drop_ratio (at end),
  left_drop_sharpness (bases until 50% threshold crossed), right_drop_sharpness
- INV-specific: split_reads count, strand_bias (max of fwd_ratio, rev_ratio)

How to Interpret the Numbers

  DEL confirmation signals (stronger → weaker):
  1. left_drop_ratio < 0.3: >70% depth drop at left breakpoint (strongest single signal)
  2. zero_fraction > 0.30: at least 30% of bases at near-zero depth
  3. median_ratio < 0.30: median depth in SV region less than 30% of flank median
  4. sharpness < 50: depth transition occurs within 50bp (clean breakpoint)
  5. p90 < flank_mean * 0.5: even the noisiest 10% of bases are below half flank depth

  DEL contradiction signals:
  1. median_ratio > 0.80: depth unchanged (not a deletion)
  2. zero_fraction < 0.01 and p90 > flank_mean * 0.5: consistent depth throughout
  3. Both left_drop_ratio and right_drop_ratio near 1.0: no breakpoint transitions

  DUP confirmation signals:
  1. mean_ratio > 1.5: depth increased >50% in SV region
  2. median_ratio > 1.3: consistent duplication signal

  INV confirmation signals:
  1. strand_bias > 0.65: >65% of split reads on one strand
  2. split_reads > 10: sufficient evidence
  3. median_ratio ~1.0: depth unchanged (balanced event, as expected)

  Boundary/edge artifacts to watch for:
  - One flank has very low depth (<10x): SV near contig boundary or assembly gap
  - right_drop_ratio > 2.0: depth spike at right breakpoint (reads piling up at deletion edge)
  - asymmetric flanks: left_flank ≠ right_flank (possible repeat or copy number variation)

Validated Cases

  Case 1: IMX2600 DEL 6.2kb (LAR CONFIRMED both)
  - Left flank 372x, Right flank 362x, SV median 98x
  - Left drop 0.219 (78% drop), sharpness=0 (immediate)
  - median_ratio 0.267, zero_fraction 0.02%
  - Interpretation: Real deletion. Depth drops 78% but not to zero because
    PacBio HiFi reads (15-20kb) span the 6.2kb deletion and maintain partial coverage.
    The sharp left breakpoint and 73% median reduction confirm the DEL.
  - Previous version called CONTRADICTS (threshold too strict at <0.05 median_ratio).

  Case 2: Makgeolli CP025104.1 (LAR CONFIRMED both, but anomalous)
  - Left flank 377x, Right flank 405x, SV median 95x
  - Left drop 0.0 (no drop), zero_fraction 0.0%, p90=387x
  - Interpretation: NOT a deletion. No depth change at all. Despite 3 callers
    and LAR both confirming, the raw depth contradicts. Possible complex SV
    misclassified as DEL. Backtrack caught this immediately.

  Case 3: SX2 DEL (LAR CONTRADICTED, both assemblers found only 196-307bp)
  - Left flank 127x, Right flank 3x, SV median 9x
  - Left drop 0.014 (98.6% drop), zero_fraction 44.4%, right flank at 3x (boundary)
  - Interpretation: Real deletion near contig boundary. The right flank at 3x
    suggests the deletion extends to a contig end or assembly gap. LAR may have
    failed because of the boundary. Backtrack numbers suggest this is real.

Files Modified
- valid_sv/evidence/layer_backtrack.py — complete rewrite (375→275 lines)
- valid_sv/evidence/layer_backtrack_verdicts.py — backup of version with verdicts
- valid_sv/evidence/layer_backtrack_complex.py — original minimap2 version backup

Next Steps
- Integrate backtrack report into run_validation.py scoring
- Define scoring weights for each backtrack metric per SV type
- Run backtrack on full LAR truth set to calibrate signal thresholds
- Test on DUPs and INVs

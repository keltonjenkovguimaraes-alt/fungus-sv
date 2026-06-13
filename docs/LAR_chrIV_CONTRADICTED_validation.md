# LAR Validation of CONTRADICTED Calls — Chromosome IV

**Date:** 1 June 2026  
**Method:** Regional Flye assembly of reads from SV loci  
**Objective:** Determine whether CONTRADICTED calls (T < 0.20) are false positives or real SVs with insufficient triangulation evidence

---

## Background

The FUNGUS-SV v3 pipeline classified 151 of 277 SVs as CONTRADICTED (T < 0.20). This tier indicates insufficient or contradictory evidence from the triangulation layers. However, it does not distinguish between:

- **False positives** — the SV does not exist
- **Real SVs with weak signal** — the SV exists but depth/k-mer layers cannot confirm it
- **Complex SVs** — the called type is incorrect but a structural change exists

LAR (Local Assembly Refinement) can resolve this ambiguity by assembling the actual DNA sequence at each locus.

---

## Methods

Three CONTRADICTED calls on chromosome IV were selected representing different SV types and sizes:

| Test | SV ID | Type | Size | T-Score |
|------|-------|------|------|---------|
| 1 | ICB_NC_001136.10_645234_46 | DEL | 6,269 bp | 0.148 |
| 2 | ICB_NC_001136.10_758078_69 | DUP | 1,701 bp | 0.167 |
| 3 | ICB_NC_001136.10_384877_37 | DEL | 120 bp | 0.148 |

For each SV, reads mapping to the SV region ± 3 kb were extracted, assembled with Flye, and aligned back to the S288C reference with minimap2 (asm5 preset). The CIGAR string was inspected for evidence of the called structural variant.

---

## Results

### Test 1: Large Deletion (6,269 bp) — CONFIRMED

| Metric | Value |
|--------|-------|
| Reads extracted | 189 |
| Contigs assembled | 2 |
| CIGAR evidence | **6269D** in both contigs |

The assembled CICC-1445 contigs align to S288C with a precise 6,269 bp deletion gap, matching the called SV size exactly.

**Verdict:** ✅ **REAL** — The deletion exists. The pipeline incorrectly classified it as CONTRADICTED.

**Why the pipeline missed it:** The T-score of 0.148 results from the size-stratified depth penalty combined with k-mer layer noise. The deletion is real but the triangulation evidence was insufficient to reach WEAK (T ≥ 0.20).

---

### Test 2: Duplication (1,701 bp) — COMPLEX

| Metric | Value |
|--------|-------|
| Reads extracted | 192 |
| Contigs assembled | 2 |
| CIGAR evidence | Multiple large insertions and overlapping alignments |

The assembled contigs show structural complexity in this region (gaps of 6–9 kb, multiple overlapping alignments) but no clean 1.7 kb tandem duplication.

**Verdict:** ⚠️ **COMPLEX** — A structural change exists but it is not a simple tandem duplication. The pipeline correctly flagged this as uncertain.

**Why the pipeline flagged it:** DUP calls require DHFFC > 2.0 for confirmation in haploid genomes. This region shows only modest depth change, consistent with a dispersed duplication or complex rearrangement rather than a tandem duplication.

---

### Test 3: Small Deletion (120 bp) — FALSE POSITIVE

| Metric | Value |
|--------|-------|
| Reads extracted | 201 |
| Contigs assembled | 3 |
| CIGAR evidence | Multiple small gaps (30–400 bp), none matching 120 bp |

The assembled contigs do not contain a 120 bp deletion. The region has many small indels but none correspond to the called SV.

**Verdict:** ❌ **FALSE POSITIVE** — The called 120 bp deletion does not exist. The pipeline correctly scored it CONTRADICTED.

**Why the pipeline flagged it:** Small SVs (<500 bp) have unreliable depth signals (Pedersen & Quinlan 2019: AUC drops for small events). The size-stratified penalty (0.75 for 100–500 bp) combined with weak depth evidence correctly pushed this call into CONTRADICTED.

---

## Summary

| Test | Type | Size | Pipeline | LAR | Conclusion |
|------|------|------|----------|-----|------------|
| 1 | DEL | 6.3 kb | CONTRADICTED | REAL | Pipeline too conservative — missed a real deletion |
| 2 | DUP | 1.7 kb | CONTRADICTED | COMPLEX | Pipeline correct — not a simple DUP |
| 3 | DEL | 120 bp | CONTRADICTED | FALSE | Pipeline correct — call is spurious |

**2/3 CONTRADICTED calls were correctly classified by the pipeline. One real deletion was missed due to overly strict thresholds.**

---

## Why LAR Must Be Run Separately

LAR is intentionally kept as a **standalone tool** rather than integrated into the automated pipeline for several reasons:

### 1. Computational Cost

| Metric | Per SV |
|--------|--------|
| Time | 2–15 minutes |
| RAM | 200–500 MB |
| Disk (temporary) | 100–500 MB |

Running LAR on all 277 SVs would take **~10–70 hours** and generate **~50 GB of temporary files**. It is impractical as an automated step.

### 2. Scientific Judgment Required

LAR provides assembled contigs and CIGAR strings — not a simple yes/no score. Interpreting the results requires human judgment:

- **Test 1 (DEL):** Clear — the CIGAR contains exactly the called deletion size
- **Test 2 (DUP):** Ambiguous — structural complexity exists but the called type may be wrong
- **Test 3 (DEL):** Clear — no matching deletion found

These nuanced interpretations cannot be automated reliably.

### 3. Targeted Application

LAR is most valuable when applied to:
- **Key candidate SVs** for publication
- **CONTRADICTED calls** that might be real (e.g., large deletions, inversions)
- **SVs in genes of biological interest**
- **SVs selected for experimental validation (PCR primer design)**

Running LAR on the ~10–20 most important SVs provides definitive proof where it matters most, without the computational burden of validating all 277 calls.

### 4. Independent Validation

As a standalone method, LAR serves as **independent validation** — it does not influence the pipeline's T-scores. This separation ensures:

- The triangulation system can be evaluated independently
- LAR results can be reported as orthogonal confirmation
- Users can choose which SVs to validate based on their research goals

---

## Recommendations

1. **Run LAR on all CONTRADICTED calls >5 kb** — These are most likely to be real (Test 1)
2. **Run LAR on all DUP calls** — To distinguish true duplications from complex rearrangements
3. **Use LAR for PCR primer design** — Assembled contigs provide exact breakpoint sequences
4. **Report LAR results separately** — As orthogonal validation, not as part of the automated pipeline

---

## References

- Kolmogorov et al. (2019). Flye. *Nature Biotechnology*, 37:540–546.
- Li H. (2018). Minimap2. *Bioinformatics*, 34:3094–3100.
- Pedersen & Quinlan (2019). Duphold. *GigaScience*, 8(4):giz040.
- David et al. (2024). Manual curation of SVs. *Genome Biology and Evolution*, 16(4):evae049.

# Understanding SV Direction in FUNGUS-SV

## The Reference-Query Relationship

FUNGUS-SV calls structural variants by aligning **query reads** (CICC-1445) to a **reference genome** (e.g., S288C). All SV calls describe differences **in the query relative to the reference**.

| Term | Meaning |
|------|---------|
| **Query** | The sample being studied (CICC-1445) — reads are from this strain |
| **Reference** | The genome used for alignment (S288C, BJ4, etc.) — the baseline |
| **DEL** | Region present in reference but absent/divergent in query |
| **DUP** | Region present in higher copy number in query vs reference |
| **INV** | Region present in both but inverted in query vs reference |

---

## Example: YBL005W-B Ty2 Element (chrII: 221,032–226,953)

### The Finding

The pipeline called a **5.9 kb DEL** in CICC-1445 vs S288C. But the same deletion does NOT appear when CICC-1445 is compared to BJ4, IMX2600, Makgeolli, or SX2.

### What This Means

| Comparison | Gap Visible? | Interpretation |
|-----------|:---:|----------------|
| CICC-1445 vs **S288C** | ✅ Yes | S288C has the Ty2; CICC-1445 does not → called DEL |
| CICC-1445 vs **BJ4** | ❌ No | BJ4 lacks the Ty2 (same as CICC-1445) → no variant |
| CICC-1445 vs **IMX2600** | ❌ No | IMX2600 lacks the Ty2 |
| CICC-1445 vs **Makgeolli** | ❌ No | Makgeolli lacks the Ty2 |
| CICC-1445 vs **SX2** | ❌ No | SX2 lacks the Ty2 |

### The Biological Reality

**S288C is the outlier** — it carries a Ty2 retrotransposon insertion that the other five strains (including CICC-1445) lack. The pipeline correctly identifies the structural difference between query and reference, but the evolutionary interpretation (insertion vs deletion) requires multi-reference comparison.

### Why It's Called a Deletion, Not an Insertion

The pipeline reports SVs from the **reference's perspective only**. Since S288C is the reference and has the Ty2 element, the absence of this element in CICC-1445 reads is reported as a deletion. The pipeline does not infer evolutionary history — it only reports: *"This region exists in the reference but not in the query."*

If the comparison were reversed (CICC-1445 as reference, S288C reads as query), the same event would be called an **insertion**.

---

## The Value of Multi-Reference Comparison

Comparing one query against multiple references reveals:

- **Reference-specific variants** — SVs that appear only against one reference (e.g., S288C Ty2 insertion)
- **Query-specific variants** — SVs that appear against ALL references (true CICC-1445-specific features)
- **Shared alleles** — Regions where query and some references match, others differ

This is why FUNGUS-SV was tested against 5 reference strains — to distinguish query-specific SVs from reference-specific polymorphisms.

---

## Practical Implications

### For PCR Validation

When designing primers to validate a DEL:

1. Check if the "deletion" appears against multiple references
2. If it only appears against one reference, the reference likely has an insertion — design primers to amplify the reference-specific region
3. If it appears against all references, the query truly has a deletion — design primers spanning the deleted region

### For Biological Interpretation

- **DEL called only vs S288C** → likely an S288C-specific insertion (or ancestral loss in other strains)
- **DEL called vs all 5 references** → likely a true CICC-1445 deletion
- **INV called vs all references** → likely a CICC-1445-specific inversion

---

## S288C-Specific Features Identified

| Locus | Type | Size | Present In |
|-------|------|------|------------|
| YBL005W-B (Ty2) | Insertion | 5.9 kb | S288C only |
| YBR012W-B (Ty2) | Insertion | 5.9 kb | S288C only |
| FLO1 internal variation | Complex | ~500 bp | S288C differs from all others |

---

## References

- 1002 Yeast Genomes Project (Peter et al., 2018) — documents natural Ty element variation across *S. cerevisiae* strains
- Yue et al. (2017) — de novo assembly of 12 strains including S288C; documents S288C-specific features
- Jeffares et al. (2017) — transposon insertion catalog across yeast populations

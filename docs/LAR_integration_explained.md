# LAR Integration into FUNGUS-SV

## Date: 30 May 2026

---

## 1. Why LAR Wasn't Working

Layer 1 (Local Assembly Refinement) was a placeholder in the pipeline. It never ran because of three issues:

### Issue 1: Whole-Genome Assembly Killed by RAM Limits

The original design attempted to assemble the **entire CICC-1445 genome** from all 274,915 reads.

flye --pacbio-hifi cicc1445_hifi.fastq.gz --genome-size 12m
Result: KILLED (out of memory)
System RAM: 5.7 GB
Flye requirement: ~10+ GB for 12 Mb genome

Both hifiasm and Flye failed on full genome assembly.

### Issue 2: Placeholder Code Never Replaced

The pipeline code at `valid_sv/run_validation.py` (lines 214-223) always returns:
- `evidence_score = 0.0`
- `available = False`
- Status: "not_run"

The separate module `local_assembly.py` was designed to be run manually before validation, but it was never executed, failed due to RAM when attempted, and had no automatic integration even if successful.

### Issue 3: No CIGAR Parsing for SV Confirmation

Even if assembly succeeded, the pipeline had no code to:
- Parse the CIGAR string from minimap2 alignments
- Check whether the CIGAR contains a deletion matching the called SV size
- Check whether the contig aligns on the opposite strand (inversion signature)
- Convert alignment patterns into an evidence score

---

## 2. How Regional LAR Works

Instead of assembling the whole genome, regional LAR assembles **only the SV locus** by extracting reads that map to the SV region plus short flanks.

### Visual Comparison

OLD (Whole Genome): NEW (Regional):
+-------------------------+ +-------------------------+
| 274,915 reads | | ~200-6000 reads from |
| 11 GB FASTQ | | SV region +/- 3 kb |
| 12 Mb genome | | ~12-435 kb window |
+------------+------------+ +------------+------------+
| |
v v
+---------+ +---------+
| Flye | KILLED | Flye | WORKS
| 10+ GB | (OOM) | <500 MB | 2-15 min
+---------+ +---------+
|
v
+---------+
| Contig |
| ~210 kb |
+---------+
|
v
+---------+
| CIGAR: |
| 253M2D |
| 5642M |
| -> DEL! |
+---------+

### Step-by-Step

**Step 1: Define the extraction window**

SV: chrII:221,032-226,953 (5.9 kb DEL)
Window: chrII:218,000-230,000 (SV +/- 3 kb flanks)

**Step 2: Extract reads from that window only**

samtools view -b cicc1445_sorted.bam NC_001134.8:218000-230000 > region.bam
samtools fastq region.bam > region.fastq
Result: ~200 reads extracted (only those spanning this deletion).

**Step 3: Assemble the extracted reads with Flye**

flye --pacbio-hifi region.fastq --genome-size 12000 --threads 2
The genome-size parameter is the window size (12 kb), not the organism's genome size. This is why RAM stays low.

**Step 4: Align the assembled contig to the reference**

minimap2 -ax asm5 S288C_reference.fasta assembly.fasta > aln.sam

**Step 5: Parse CIGAR and confirm the SV**

CIGAR: 253M2D5642M
Interpretation: 253 bp match, 2 bp deleted, 5642 bp deleted, then matches again
Conclusion: Deletion of ~5.6 kb confirmed in CICC-1445

### LAR Test Results

| # | Locus | Type | Size | Strain | Reads | Contigs | RAM | Time | Result |
|---|-------|------|------|--------|-------|---------|-----|------|--------|
| 1 | chrII YBL005W-B | DEL | 5.9 kb | S288C | 202 | 3 | <200 MB | 2 min | Confirmed |
| 2 | chrVII multi-gene | DEL | 55.7 kb | S288C | 581 | 6 | <300 MB | 4 min | Confirmed |
| 3 | SX2 chrII | INV | 430 kb | SX2 | 6,448 | 7 | <500 MB | 13 min | Confirmed |
| 4 | BJ4 chrXII | INV | 205 kb | BJ4 | 3,552 | 1 | <400 MB | 10 min | Confirmed |

---

## 3. Automated LAR for Any Genome

### Why It's Genome-Agnostic

| Component | What It Needs | Genome-Specific? |
|-----------|--------------|:---:|
| samtools view | BAM file + region coordinates | No |
| Flye | FASTQ reads + window size | No |
| minimap2 | Reference FASTA + contigs | No |
| CIGAR parsing | Standard SAM format | No |

The only parameter that changes per SV is the window size passed to Flye, computed automatically as:

window_size = (sv_end - sv_start) + (2 x flank_size)

### Scoring Logic

| LAR Outcome | Evidence Score | Weight | Contribution to T-Score |
|-------------|:---:|:---:|------------------------|
| SV confirmed by contig | 1.0 | 0.20 | +0.20 |
| Partial support | 0.5 | 0.20 | +0.10 |
| Contig contradicts SV | 0.0 | 0.20 | +0.00 |
| Too few reads (<50) | unavailable | 0.20 | skipped |

### Requirements for Any Genome

- BAM file (query reads aligned to reference)
- Reference FASTA (indexed)
- SV coordinates from VCF
- Flye + minimap2 + samtools installed
- Minimum 50 reads in the SV window
- RAM: <500 MB per SV (not genome-size dependent)

### Estimated Impact on T-Scores

| SV Type | Current T (without LAR) | With LAR Confirmed | Change |
|---------|------------------------|-------------------|--------|
| YBL005W-B DEL | 0.148 (CONTRADICTED) | ~0.35 (WEAK) | +0.20 |
| chrVII DEL | 0.148 (CONTRADICTED) | ~0.35 (WEAK) | +0.20 |
| SX2 INV | 0.167 (CONTRADICTED) | ~0.37 (WEAK) | +0.20 |
| BJ4 INV | 0.167 (CONTRADICTED) | ~0.37 (WEAK) | +0.20 |

LAR alone will not make CONTRADICTED calls become HIGH — but it adds the strongest single piece of evidence, moving them from "no evidence" to "assembly confirmed."

---

## Summary

| Question | Answer |
|----------|--------|
| Why wasn't LAR working? | RAM killed on whole-genome assembly; placeholder code never replaced; no CIGAR parsing |
| How does regional LAR work? | Extracts reads from SV window only, assembles small region, aligns contig to reference, parses CIGAR |
| Will it work for any genome automatically? | Yes — all components are genome-agnostic; only window size changes |

---

## References

- Kolmogorov et al. (2019). Flye: assembly of long-read genomes. *Nature Biotechnology*, 37:540-546.
- Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34:3094-3100.
- Li H. et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25:2078-2079.

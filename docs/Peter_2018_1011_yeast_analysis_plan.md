# Analysis Plan: Cross-Referencing FUNGUS-SV with the 1011 Yeast Genomes (Peter et al. 2018)

## Objective

Validate FUNGUS-SV SVs by comparing them against the largest available population-scale yeast genomic dataset. This provides external, independent evidence for whether our called SVs represent genuine polymorphic loci in *S. cerevisiae*.

## Data Already Downloaded

- **File:** `genesMatrix_CopyNumber.tab.gz` (980 KB)
- **Location:** `/home/kelto/fungus-sv/data/genesMatrix_CopyNumber.tab.gz`
- **Content:** Copy number estimates for 7,796 pangenomic ORFs across 1,011 isolates

## Key Challenge

The gene names in the Peter et al. matrix use **4-letter codes** (e.g., BFC, AGL, ABE) that do not match standard *S. cerevisiae* gene names (e.g., FLO1, HXT7) or systematic ORF names (e.g., YAL063C). A mapping table is required.

## Files Available from Peter et al. That Can Help

| File | Description | Potential Use |
|------|-------------|---------------|
| `allORFs_pangenome.fasta.gz` | Sequences of 7,796 pangenomic ORFs | BLAST against S288C to map 4-letter codes → standard names |
| `genesMatrix_PresenceAbsence.tab.gz` | Presence/absence of ORFs per isolate | Cross-reference our DELs with known absent genes |
| `1011DistanceMatrixBasedOnSNPs.tab.gz` | SNP-based distances | Place CICC-1445 in population context |
| `genesMatrix_Frameshift.tab.gz` | Frameshift mutations per gene | Functional impact of SVs affecting genes |
| `1011GWASMatrix.tar.gz` | GWAS matrix with CNVs | Independent CNV calls to compare |

## Analysis Plan

### Step 1: Build Gene Name Mapping

**Goal:** Map Peter et al. 4-letter codes to S288C systematic names (YAL063C → FLO9, etc.)

**Method:**
1. Download `allORFs_pangenome.fasta.gz` from the Peter et al. supplement
2. BLAST the pangenomic ORF sequences against the S288C reference genome (or proteome)
3. For each 4-letter code, find the best S288C hit
4. Build a mapping dictionary: `{BFC: YAL063C, AGL: YAR050W, ...}`

**Commands:**
```bash
# Download the pangenome ORFs
wget -O data/peter2018_allORFs_pangenome.fasta.gz <URL from Nature supplement>

# Make BLAST database of S288C genes
makeblastdb -in data/yeast/S288C_reference.fasta -dbtype nucl

# BLAST pangenome ORFs against S288C
blastn -query data/peter2018_allORFs_pangenome.fasta.gz -db data/yeast/S288C_reference.fasta -outfmt 6 -out data/peter_to_s288c.blast
Step 2: Cross-Reference Our Affected Genes
Goal: For each gene affected by a FUNGUS-SV DEL/DUP, check its copy number variation across the 1,011 strains.

Method:

Load our list of genes overlapping SVs (from figures/genes.bed and manual curation)

Map each gene to its Peter et al. 4-letter code

Extract copy number values from genesMatrix_CopyNumber.tab.gz

Compute population statistics: mean CN, variance, % of strains with CN=0, % with CN≠1

Expected Output Table:
Gene     S288C_Name    Mean_CN    %_Strains_CN=0    %_Strains_CN≠1    Our_SV_Type    Interpretation
BFC      FLO9          0.95       2%                15%               DEL (516bp)    Common variant
AGL      HXT7          0.88       8%                22%               DEL (5.4kb)    Frequent DEL
Step 3: Validate Specific SV Calls
Goal: Determine if our called SVs are known population variants or likely false positives.

For DEL calls:

If a gene has CN=0 in >5% of the 1,011 strains → our DEL is a common population variant → VALIDATED

If a gene has CN=1 in 100% of strains → our DEL is likely false or reference-specific → FLAG FOR REVIEW

If a gene has CN>1 in some strains → our region has known amplification → CHECK IF DUP/DEL MISCLASSIFICATION

For DUP calls:

If a gene has CN>1 in >5% of strains → our DUP is a known amplification → VALIDATED

If a gene has CN=1 in 100% of strains → our DUP is likely false → FLAG

For S288C-specific features (Ty2 insertions):

We found YBL005W-B and YBR012W-B are present in S288C but absent in CICC-1445

Check if these Ty2 elements are present in the Peter et al. pangenome

If most strains lack them → confirms S288C is the outlier → VALIDATES our interpretation

Step 4: Population Context for CICC-1445
Goal: Determine where CICC-1445 sits in the yeast population.

Method:

If possible, find CICC-1445 or a close relative in the 1,011 strains

Download 1011DistanceMatrixBasedOnSNPs.tab.gz

Check if our finding (CICC-1445 closest to BJ4/Chinese industrial strains) matches SNP-based distances

Step 5: Compare CNV Detection Sensitivity
Goal: Benchmark FUNGUS-SV against the Peter et al. CNV calls.

Method:

Extract all genes with CN≠1 in >10% of the 1,011 strains (these are common CNVs)

Check if FUNGUS-SV called SVs overlapping these genes

Compute: What fraction of known population CNVs did FUNGUS-SV detect?

This gives an independent sensitivity estimate (no spike-in needed!)

Priority Order
Priority	Task	Effort	Impact
1	Download pangenome ORFs + build name mapping	1-2 hours	Unlocks all other analyses
2	Cross-reference our affected genes with CN data	1 hour	Direct validation of our calls
3	Validate Ty2/S288C-specific features	30 min	Confirms our interpretation
4	Population context for CICC-1445	1 hour	Biological context
5	Sensitivity benchmark vs known CNVs	2 hours	Independent performance metric
Expected Outcomes
External validation — SVs matching known population variants are likely real

False positive detection — SVs in genes that never vary are likely false

Biological interpretation — Are our SVs common or rare in the species?

Pipeline calibration — Use population frequencies to adjust confidence scores

Paper figure — Population frequency histogram of our validated SVs

References
Peter, J. et al. (2018). Genome evolution across 1,011 Saccharomyces cerevisiae isolates. Nature, 556(7701), 339-344.

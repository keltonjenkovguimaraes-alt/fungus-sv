# Internal Benchmark Methodology for FUNGUS-SV

## Assembly-Based Truth Set Construction for Non-Model Organisms

### Authors: Kelton Jenkov Guimarães et al.
### Pipeline: FUNGUS-SV v0.5.1
### Date: May 2026

---

## 1. Background

Structural variant (SV) detection in non-model organisms lacks the gold-standard benchmark sets available for human genomes (GIAB, HGSVC, SMaHT). Without a truth set, it is impossible to measure precision, recall, or false discovery rates for any SV caller. This document describes a computational methodology for constructing an internal truth set using de novo genome assembly, enabling rigorous validation of SV detection pipelines in any organism.

The methodology is based on approaches used by:
- **Kronenberg et al. (2025)** *Nature Methods* — Platinum Pedigree truth set construction
- **Zhang et al. (2025)** *bioRxiv* — SMaHT mosaic SV benchmark
- **Ebert et al. (2021)** *Science* — HGSVC assembly-based SV discovery
- **Chen et al. (2023)** *Nature Communications* — DeBreak assembly validation

---

## 2. Principle

The core principle is: **SVs detected by comparing a de novo assembly to a reference genome are independent of read-alignment artifacts that affect read-based SV callers.** When both an assembly-based approach and a read-based approach (FUNGUS-SV) agree on an SV, the probability that it is a true variant is substantially higher than either method alone.

### Why Assembly-Based SVs Are More Reliable

| Evidence Type | Source | Artifacts |
|--------------|--------|-----------|
| Read-alignment SVs | Mapped reads (BAM) | Mapping errors, chimeric alignments, repeat collapse |
| Assembly-based SVs | Contigs from de novo assembly | Assembly errors (less frequent with HiFi) |
| **Concordant SVs** | **Both agree** | **Minimal — orthogonal methods confirm each other** |

The assembly-based approach has been validated in human genomics, where assembly-based SV calls from hifiasm + PAV achieve >90% concordance with GIAB truth sets (Ebert et al. 2021).

---

## 3. Methodology

### 3.1 Required Data

| Item | Format | Purpose |
|------|--------|---------|
| PacBio HiFi reads | FASTQ (.fastq.gz) | From the isolate to be analyzed |
| Reference genome | FASTA (.fasta) | Different strain/species for comparison |

### 3.2 Step 1: De Novo Genome Assembly

Assemble the HiFi reads into a haplotype-resolved genome assembly using hifiasm:

```bash
# Install hifiasm
conda create -n assembly -c bioconda hifiasm -y
conda activate assembly

# Run assembly (haploid mode for fungi)
hifiasm -o isolate_asm -t 16 --primary isolate_reads.fastq.gz

# Convert to FASTA
awk '/^S/{print ">"$2;print $3}' isolate_asm.p_ctg.gfa > isolate_asm.p_ctg.fasta
Expected output: A haploid assembly with contig N50 > 1 Mb for fungal genomes at >50× HiFi coverage.

Quality control:


# Check assembly statistics
quast isolate_asm.p_ctg.fasta -r reference.fasta -o quast_results/
3.3 Step 2: Assembly-Based SV Calling
Compare the assembled genome to the reference using SVIM-asm in haploid mode:
# Align assembly to reference
minimap2 -x asm5 -t 8 --cs reference.fasta isolate_asm.p_ctg.fasta \
    | samtools sort -@ 4 -o asm_to_ref.sorted.bam -
samtools index asm_to_ref.sorted.bam

# Call SVs from assembly comparison (haploid mode)
svim-asm haploid --min_sv_size 50 asm_to_ref.sorted.bam reference.fasta \
    results/asm_svs/
Why SVIM-asm:

Designed for both haploid and diploid assemblies (Heller & Vingron 2020)

Validated on human HGSVC data with high precision

Reports DEL, INS, INV, DUP, and complex SVs

3.4 Step 3: Read-Based SV Detection (FUNGUS-SV)
Run FUNGUS-SV on the same HiFi reads aligned to the same reference:


# Align reads to reference
minimap2 -t 8 -ax map-hifi -R '@RG\tID:isolate\tSM:sample' \
    reference.fasta isolate_reads.fastq.gz \
    | samtools sort -@ 4 -o reads_to_ref.sorted.bam -
samtools index reads_to_ref.sorted.bam

# Run FUNGUS-SV ICB consensus
python fungus_sv/core/icb.py \
    --bam reads_to_ref.sorted.bam \
    --reference reference.fasta \
    --output results/icb/ \
    --callers sniffles2 cutesv svim \
    --min-callers 2

# Run FUNGUS-SV validation
python -m valid_sv.run_validation \
    --consensus-vcf results/icb/consensus_svs.vcf \
    --bam reads_to_ref.sorted.bam \
    --reference reference.fasta \
    --fastq isolate_reads.fastq.gz \
    --output results/validation/
3.5 Step 4: Benchmark Comparison
Compare FUNGUS-SV results to the assembly-based truth set using Truvari:
# Compare ICB consensus to assembly SVs
truvari bench \
    -b results/asm_svs/variants.vcf \
    -c results/icb/consensus_svs.vcf \
    --reference reference.fasta \
    -o benchmark_results/ \
    --pctseq 0 --pctsize 0.7 --refdist 500

# Generate summary metrics
truvari summarize --input benchmark_results/ --output summary.txt
3.6 Step 5: Stratified Performance Analysis
Compute precision, recall, and F1 score stratified by:

Stratification	Why
SV type (DEL, INS, INV, DUP)	Different types have different detection difficulty
SV size (50-100 bp, 100-500 bp, 500-5000 bp, >5 kb)	Size-dependent accuracy
T-score tier (HIGH ≥0.6, MEDIUM 0.4-0.6, WEAK <0.4)	Validates T-score calibration
Repeat context (in repeat vs. unique)	Repeat regions are challenging
Number of supporting ICB callers (2 vs. 3)	Validates ICB consensus
4. Expected Outcomes
4.1 Performance Targets (Based on Human Benchmarks)
Metric	Target	Source
ICB Precision	>80%	SV-JIM (Todd et al. 2025): 3-caller consensus
ICB Recall (DEL)	>60%	Liu et al. (2024): best callers on PacBio
ICB Recall (INS)	>40%	Liu et al. (2024): insertions are harder
T≥0.6 FDR	<15%	SV-MeCa (Nkouamedjo et al. 2025)
Ablation: multi-layer > single-layer	≥15% improvement	Expected from triangulation theory
4.2 T-Score Calibration
The benchmark enables mapping T-score thresholds to empirical false discovery rates:

T-Score	Expected FDR	Action
≥0.8	<5%	Suitable for functional follow-up
≥0.6	<15%	Prioritized candidate
≥0.4	<35%	Exploratory only
<0.4	>50%	Likely false positive
These thresholds are estimates pending empirical calibration.

5. Limitations
Assembly errors can create false SVs. hifiasm assemblies of HiFi data have high accuracy (QV > 50), but errors in repetitive regions can introduce false SVs that will be misclassified as true positives.

Assembly-based callers miss complex SVs. SVIM-asm does not detect translocations or complex rearrangements. The truth set will be incomplete.

This is a silver standard, not gold. Unlike GIAB (which uses multiple technologies, pedigree information, and manual curation), this approach uses a single technology and automated methods.

Validated only for haploid genomes. Diploid or polyploid genomes require phasing before assembly, which introduces additional complexity.

6. Validation Against Published Methods
To ensure the benchmark is reliable, cross-validate against orthogonal approaches:

Method	Tool	Comparison
Different assembly-based caller	PAV	Compare to SVIM-asm results
Different assembler	Flye	Compare to hifiasm assembly
Read-only validation	SVvalidation (Zheng 2024)	Independent breakpoint confirmation
High concordance (>80%) between methods increases confidence in the truth set.

7. Software Requirements
Tool	Version	Installation
hifiasm	≥0.19	conda install -c bioconda hifiasm
minimap2	≥2.30	conda install -c bioconda minimap2
samtools	≥1.17	conda install -c bioconda samtools
svim-asm	≥1.0	conda install -c bioconda svim-asm
Truvari	≥5.0	conda install -c bioconda truvari
quast	≥5.0	conda install -c bioconda quast
8. References
Kronenberg, Z. et al. (2025). The Platinum Pedigree: A long-read benchmark for genetic variants. Nature Methods.

Zhang, Y. et al. (2025). Comprehensive benchmarking of somatic structural variant detection at ultra-low allele fractions. bioRxiv.

Ebert, P. et al. (2021). Haplotype-resolved diverse human genomes and integrated analysis of structural variation. Science, 372.

Chen, Y. et al. (2023). Deciphering the exact breakpoints of structural variations using long sequencing reads with DeBreak. Nature Communications, 14.

Heller, D. & Vingron, M. (2020). SVIM-asm: structural variant detection from haploid and diploid genome assemblies. Bioinformatics, 36.

Liu, Y.H. et al. (2024). Tradeoffs in alignment and assembly-based methods for structural variant detection with long-read sequencing data. Nature Communications, 15.

Todd, C. et al. (2025). SV-JIM: detailed pairwise structural variant calling using long-reads and genome assemblies. Methods, 234.

Nkouamedjo Fankep, R.C. et al. (2025). SV-MeCa: an XGBoost-based meta-caller approach for structural variant calling. BMC Bioinformatics, 26.

Zheng, Y. & Shang, X. (2024). SVvalidation: A long-read-based validation method for genomic structural variation. PLOS ONE, 19.

This methodology is part of the FUNGUS-SV pipeline documentation. The pipeline is available at https://github.com/keltonjenkovguimaraes-alt/fungus-sv

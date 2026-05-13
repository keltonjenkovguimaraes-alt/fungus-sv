
FUNGUS-SV Methods
Sequencing
Platform: PacBio Revio HiFi

Read length: 15-21 kb

Coverage: 58× (197,830 reads)

Alignment
minimap2 v2.30, map-hifi preset

Mapping rate: 99.91%

ICB (Iterative Consensus Builder)
Callers: pbsv v2.11.0, Sniffles2 v2.8.0, cuteSV v1.0.8

Consensus: reciprocal overlap >= 0.5, >= 2 callers

Iterations: 3

LAR (Local Assembly Refinement)
Read extraction: samtools (SV +/- 5kb)

Assembly: Flye v2.9.6

Validation: 99.85% identity

Annotation
Gene overlap: custom Python parser

Impact: HIGH/MODERATE/LOW

Software Versions
Tool	Version
minimap2	2.30
pbsv	2.11.0
Sniffles2	2.8.0
cuteSV	1.0.8
Flye	2.9.6
Snakemake	8.0+
EOF	
Create .gitignore
cat > .gitignore << 'EOF'
data/
results/
.fastq
*.fasta
*.fna
*.bam
*.bai
.vcf
pycache/
*.pyc
.env/
.vscode/
.idea/

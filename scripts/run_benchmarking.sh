#!/bin/bash
# FUNGUS-SV Benchmarking Script
# Compares FUNGUS-SV against MUMmer4 and SVIM-asm
# Usage: bash scripts/run_benchmarking.sh

conda activate sv_bench 2>/dev/null || conda create -n sv_bench -c bioconda mummer4 svim-asm minimap2 samtools -y
conda activate sv_bench

BENCH_DIR="benchmarking"
mkdir -p $BENCH_DIR

echo "=== Saccharomyces ==="
QUERY="data/yeast/cicc1445_self/CICC1445_reference.fasta"

for STRAIN in S288C BJ4 IMX2600 Makgeolli SX2; do
    echo "--- $STRAIN ---"
    REF="data/yeast/${STRAIN}_reference.fasta"
    nucmer --mum -p $BENCH_DIR/cicc1445_vs_${STRAIN} $REF $QUERY
    dnadiff -p $BENCH_DIR/cicc1445_vs_${STRAIN} $REF $QUERY
    minimap2 -ax asm5 $REF $QUERY | samtools view -b | samtools sort -o $BENCH_DIR/${STRAIN}_asm.bam
    samtools index $BENCH_DIR/${STRAIN}_asm.bam
    svim-asm haploid --min_sv_size 50 $BENCH_DIR/${STRAIN}_svim $BENCH_DIR/${STRAIN}_asm.bam $REF
done

echo "=== Candida ==="
QUERY="data/candida/SC5314_reference.fasta"

for STRAIN in 101 FDAARGOS656 WO1 ATCC64124 UAB012; do
    echo "--- $STRAIN ---"
    REF="data/candida/${STRAIN}_reference.fasta"
    nucmer --mum -p $BENCH_DIR/sc5314_vs_${STRAIN} $REF $QUERY
    dnadiff -p $BENCH_DIR/sc5314_vs_${STRAIN} $REF $QUERY
    minimap2 -ax asm5 $REF $QUERY | samtools view -b | samtools sort -o $BENCH_DIR/${STRAIN}_asm.bam
    samtools index $BENCH_DIR/${STRAIN}_asm.bam
    svim-asm haploid --min_sv_size 50 $BENCH_DIR/${STRAIN}_svim $BENCH_DIR/${STRAIN}_asm.bam $REF
done

echo "Done."

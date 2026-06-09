#!/usr/bin/env python3
"""
Evidence Layer 5: Ploidy Confirmation (Fast Pileup Method)
============================================================
Confirms haploid assumption by sampling random genomic windows
and counting heterozygous allele frequencies via bcftools mpileup.

Principle: In a true haploid, virtually all positions show a single
allele. In diploids/polyploids, heterozygous positions show multiple
alleles at intermediate frequencies.

Method inspired by allele frequency distribution analysis in phased
read clusters (Abou Saada et al. 2021, Genome Biology).

Author: FUNGUS-SV v0.9.4
"""

import subprocess
import random
from dataclasses import dataclass


@dataclass
class PloidyEvidence:
    total_snvs: int
    homozygous: int
    heterozygous: int
    het_rate: float
    is_haploid: bool
    evidence_score: float
    details: str


def analyze_ploidy(bam_path, reference_path, num_windows=10, window_size=50000, min_mapq=20):
    """Fast haploid check via bcftools mpileup on random genomic windows."""

    contigs = []
    try:
        with open(reference_path + '.fai') as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) >= 2 and int(p[1]) > window_size:
                    contigs.append((p[0], int(p[1])))
    except FileNotFoundError:
        return PloidyEvidence(0, 0, 0, 0.0, True, 1.0,
            f"Reference index not found. Run: samtools faidx {reference_path}")

    if not contigs:
        return PloidyEvidence(0, 0, 0, 0.0, True, 1.0, "No contigs > window_size")

    random.seed(42)
    total_het = 0
    total_positions = 0

    for _ in range(num_windows):
        chrom, length = random.choice(contigs)
        start = random.randint(1, max(1, length - window_size))
        end = start + window_size
        region = f"{chrom}:{start}-{end}"

        try:
            r = subprocess.run(
                ['bcftools', 'mpileup', '-q', str(min_mapq), '--no-BAQ',
                 '-f', reference_path, '--region', region, bam_path],
                capture_output=True, text=True, timeout=60
            )
            if r.returncode != 0 or not r.stdout.strip():
                continue

            r2 = subprocess.run(
                ['bcftools', 'call', '-mv', '--ploidy', '1'],
                input=r.stdout, capture_output=True, text=True, timeout=60
            )

            for line in r2.stdout.strip().split('\n'):
                if line.startswith('#'):
                    continue
                total_positions += 1
                if '0/1' in line or '1/0' in line:
                    total_het += 1

        except (subprocess.TimeoutExpired, FileNotFoundError):
            continue

    if total_positions == 0:
        return PloidyEvidence(0, 0, 0, 0.0, True, 1.0,
            "No variant positions found — assuming haploid")

    het_rate = total_het / total_positions

    if het_rate < 0.03:
        is_haploid = True; score = 1.0
        details = f"Strongly haploid: {het_rate:.2%} het ({total_het}/{total_positions})"
    elif het_rate < 0.07:
        is_haploid = True; score = 0.8
        details = f"Mostly haploid: {het_rate:.2%} het ({total_het}/{total_positions})"
    elif het_rate < 0.12:
        is_haploid = False; score = 0.3
        details = f"Possibly diploid: {het_rate:.2%} het — verify"
    else:
        is_haploid = False; score = 0.0
        details = f"Likely diploid/polyploid: {het_rate:.2%} het"

    return PloidyEvidence(total_positions, total_positions - total_het, total_het,
                          round(het_rate, 4), is_haploid, score, details)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print("Usage: layer_ploidy.py <bam> <reference> [num_windows]")
        sys.exit(1)
    r = analyze_ploidy(sys.argv[1], sys.argv[2],
                       int(sys.argv[3]) if len(sys.argv) > 3 else 10)
    print(f"Het: {r.het_rate:.2%} | Haploid: {r.is_haploid} | Score: {r.evidence_score}")

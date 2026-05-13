#!/usr/bin/env python3
"""
Evidence Layer: Ploidy Confirmation via SNV Analysis
======================================================
Confirms the haploid assumption by analyzing heterozygous
SNV rates across the genome. A truly haploid organism
should show <2% heterozygous calls (consistent with
paralogous regions or assembly errors, not true diploidy).

Uses Longshot for SNV calling from PacBio HiFi data.

Author: VALID-SV / FUNGUS-SV
"""

import subprocess
import re
from dataclasses import dataclass
from typing import Dict, Tuple


@dataclass
class PloidyEvidence:
    """Ploidy analysis results."""
    total_snvs: int
    homozygous: int
    heterozygous: int
    het_rate: float
    is_haploid: bool
    evidence_score: float
    details: str


def run_longshot(bam_path: str, reference_path: str, 
                 output_vcf: str, threads: int = 4) -> str:
    """Run Longshot SNV caller for PacBio HiFi data."""
    result = subprocess.run([
        'longshot', '--bam', bam_path, '--ref', reference_path,
        '--out', output_vcf, '--threads', str(threads)
    ], capture_output=True, text=True, timeout=600)
    
    if result.returncode != 0:
        raise RuntimeError(f"Longshot failed: {result.stderr}")
    
    return output_vcf


def analyze_ploidy(vcf_path: str) -> PloidyEvidence:
    """
    Analyze heterozygous SNV rate to confirm haploidy.
    
    In a haploid organism, >98% of SNV calls should be homozygous.
    Heterozygous calls >2% suggest diploid regions or sample mixture.
    """
    total = 0
    homozygous = 0
    heterozygous = 0
    
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            
            # Check genotype field
            fmt = parts[8].split(':')
            gt_idx = None
            for i, field in enumerate(fmt):
                if field == 'GT':
                    gt_idx = i
                    break
            
            if gt_idx is None:
                continue
            
            gt = parts[9].split(':')[gt_idx]
            
            total += 1
            if '0/1' in gt or '1/0' in gt or '0|1' in gt or '1|0' in gt:
                heterozygous += 1
            elif '1/1' in gt or '0/0' in gt or '1|1' in gt or '0|0' in gt:
                homozygous += 1
    
    if total == 0:
        return PloidyEvidence(
            total_snvs=0, homozygous=0, heterozygous=0,
            het_rate=0.0, is_haploid=False, evidence_score=0.0,
            details="No SNVs detected — check data quality"
        )
    
    het_rate = heterozygous / total
    
    # Haploid expectation: <2% heterozygous
    if het_rate < 0.02:
        is_haploid = True
        score = 1.0
        details = (f"Strongly haploid: {het_rate:.2%} heterozygous "
                  f"({heterozygous}/{total} SNVs)")
    elif het_rate < 0.05:
        is_haploid = True
        score = 0.7
        details = (f"Mostly haploid: {het_rate:.2%} heterozygous "
                  f"({heterozygous}/{total} SNVs) — some paralogous regions")
    elif het_rate < 0.10:
        is_haploid = False
        score = 0.3
        details = (f"Possibly diploid regions: {het_rate:.2%} heterozygous "
                  f"({heterozygous}/{total} SNVs)")
    else:
        is_haploid = False
        score = 0.0
        details = (f"Likely diploid or mixed sample: {het_rate:.2%} heterozygous "
                  f"({heterozygous}/{total} SNVs)")
    
    return PloidyEvidence(
        total_snvs=total,
        homozygous=homozygous,
        heterozygous=heterozygous,
        het_rate=round(het_rate, 4),
        is_haploid=is_haploid,
        evidence_score=round(score, 4),
        details=details
    )


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print("Usage: layer_ploidy.py <bam> <reference> [threads]")
        sys.exit(1)
    
    vcf = 'results/variants/snvs/longshot_snvs.vcf'
    run_longshot(sys.argv[1], sys.argv[2], vcf,
                int(sys.argv[3]) if len(sys.argv) > 3 else 4)
    
    result = analyze_ploidy(vcf)
    print(f"\n{'='*60}")
    print(f"  VALID-SV: Ploidy Analysis")
    print(f"{'='*60}")
    print(f"  Total SNVs: {result.total_snvs:,}")
    print(f"  Homozygous: {result.homozygous:,} ({result.homozygous/result.total_snvs*100:.1f}%)")
    print(f"  Heterozygous: {result.heterozygous:,} ({result.het_rate:.2%})")
    print(f"  Haploid: {result.is_haploid}")
    print(f"  Score: {result.evidence_score:.3f}")
    print(f"  {result.details}")
    print(f"{'='*60}\n")

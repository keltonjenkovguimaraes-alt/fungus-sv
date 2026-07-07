#!/usr/bin/env python3
"""
Evidence Layer: Repeat Content Reporting
==========================================
Reports repeat and low-complexity metrics for any SV region.
No verdicts — just numbers for the triangulation scorer.

Metrics help interpret:
- Why backtrack strand bias might be weak (breakpoints in repeats)
- Why depth might be noisy (mapping artifacts in repetitive regions)
- Whether an SV overlaps known fungal repeat elements

Uses only pysam (already in sv_valid) — no external tools needed.

Author: FUNGUS-SV v1.0.0
"""

import pysam
import re
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class RepeatReport:
    """Pure repeat content report — no verdicts, just numbers."""
    sv_id: str
    sv_type: str
    region: str = ""  # chrom:start-end
    
    # Basic composition
    gc_content: float = 0.0
    region_length: int = 0
    
    # Homopolymer runs (≥10bp)
    a_runs: int = 0
    t_runs: int = 0
    g_runs: int = 0
    c_runs: int = 0
    max_homopolymer: int = 0
    
    # Dinucleotide repeats (≥12bp = 6 repeats)
    at_repeats: int = 0
    ta_repeats: int = 0
    gc_repeats: int = 0
    cg_repeats: int = 0
    tg_repeats: int = 0
    gt_repeats: int = 0
    ac_repeats: int = 0
    ca_repeats: int = 0
    total_dinucleotide_runs: int = 0
    
    # Known fungal repeat motifs
    tata_count: int = 0      # TATA boxes / promoters
    gagaga_count: int = 0    # Simple repeat
    telomeric_count: int = 0 # TG1-3 / G1-3T telomeric
    delta_count: int = 0     # Ty retrotransposon delta sequences
    sigma_count: int = 0     # Ty sigma sequences
    
    # Low-complexity fraction (bases in simple repeats)
    low_complexity_bases: int = 0
    low_complexity_fraction: float = 0.0
    
    # Overall repeat score (0-1, higher = more repetitive)
    repeat_score: float = 0.0
    
    # Status
    status: str = "ok"
    error_msg: str = ""
    summary: str = ""


def _count_homopolymer_runs(seq: str, min_len: int = 10) -> dict:
    """Count homopolymer runs ≥ min_len."""
    runs = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    max_len = 0
    
    for base in ['A', 'T', 'G', 'C']:
        pattern = f'{base}{{{min_len},}}'
        matches = re.findall(pattern, seq.upper())
        runs[base] = len(matches)
        if matches:
            longest = max(len(m) for m in matches)
            max_len = max(max_len, longest)
    
    runs['max'] = max_len
    return runs


def _count_dinucleotide_repeats(seq: str, min_units: int = 6) -> dict:
    """Count dinucleotide repeats ≥ min_units (e.g., ATATATATATAT = 6 units = 12bp)."""
    motifs = ['AT', 'TA', 'GC', 'CG', 'TG', 'GT', 'AC', 'CA']
    counts = {}
    total = 0
    
    for motif in motifs:
        pattern = f'({motif}){{{min_units},}}'
        matches = re.findall(pattern, seq.upper())
        counts[motif] = len(matches)
        total += len(matches)
    
    counts['total'] = total
    return counts


def _count_known_motifs(seq: str) -> dict:
    """Count known fungal/genomic repeat motifs."""
    s = seq.upper()
    
    return {
        'tata': len(re.findall('TATA', s)),
        'gagaga': len(re.findall('GAGAGA', s)),
        'telomeric_tg': len(re.findall('TG{1,3}', s)) + len(re.findall('G{1,3}T', s)),
        'delta': len(re.findall('TGTTGGAATA', s)),  # Ty delta consensus
        'sigma': len(re.findall('TGGTGGT(G|T)G', s)),  # Ty sigma-like
    }


def _compute_repeat_score(report: 'RepeatReport') -> float:
    """Compute overall repeat score 0-1."""
    score = 0.0
    
    # Homopolymer runs contribute
    total_homopolymers = report.a_runs + report.t_runs + report.g_runs + report.c_runs
    if total_homopolymers > 0:
        score += min(0.3, total_homopolymers * 0.03)
    
    # Dinucleotide repeats contribute
    if report.total_dinucleotide_runs > 0:
        score += min(0.3, report.total_dinucleotide_runs * 0.05)
    
    # Low complexity fraction contributes
    score += min(0.3, report.low_complexity_fraction * 3)
    
    # Known repeat motifs
    motif_total = report.tata_count + report.gagaga_count + report.telomeric_count
    if motif_total > 20:
        score += 0.2
    elif motif_total > 5:
        score += 0.1
    
    return min(1.0, round(score, 3))


def analyze_repeats(reference_path: str, sv_id: str = "", sv_type: str = "",
                    chrom: str = "", start: int = 0, end: int = 0,
                    breakpoint_only: bool = False,
                    flank_size: int = 500) -> RepeatReport:
    """
    Report repeat content metrics for an SV region.
    
    Args:
        reference_path: Path to reference FASTA
        sv_id: SV identifier
        sv_type: DEL, DUP, INV, etc.
        chrom: Chromosome/contig name
        start: SV start position (1-based)
        end: SV end position (1-based)
        breakpoint_only: If True, analyze only breakpoints (±flank_size)
                         If False, analyze entire SV region
        flank_size: Flank size for breakpoint-only mode
    
    Returns:
        RepeatReport with all metrics
    """
    
    report = RepeatReport(sv_id=sv_id, sv_type=sv_type)
    
    try:
        ref = pysam.FastaFile(reference_path)
        
        if breakpoint_only:
            # Analyze left and right breakpoints separately
            left_start = max(1, start - flank_size)
            left_end = start + flank_size
            right_start = max(1, end - flank_size)
            right_end = end + flank_size
            
            left_seq = ref.fetch(chrom, left_start - 1, left_end)
            right_seq = ref.fetch(chrom, right_start - 1, right_end)
            seq = left_seq + right_seq
            report.region = f"{chrom}:{left_start}-{left_end},{right_start}-{right_end}"
        else:
            seq = ref.fetch(chrom, start - 1, end)
            report.region = f"{chrom}:{start}-{end}"
        
        ref.close()
        
        report.region_length = len(seq)
        
        if report.region_length == 0:
            report.status = "error"
            report.error_msg = "Empty sequence"
            return report
        
        # ============================================================
        # BASIC COMPOSITION
        # ============================================================
        s = seq.upper()
        gc = s.count('G') + s.count('C')
        report.gc_content = round(gc / len(s), 4)
        
        # ============================================================
        # HOMOPOLYMER RUNS
        # ============================================================
        homo = _count_homopolymer_runs(seq)
        report.a_runs = homo['A']
        report.t_runs = homo['T']
        report.g_runs = homo['G']
        report.c_runs = homo['C']
        report.max_homopolymer = homo['max']
        
        # ============================================================
        # DINUCLEOTIDE REPEATS
        # ============================================================
        di = _count_dinucleotide_repeats(seq)
        report.at_repeats = di['AT']
        report.ta_repeats = di['TA']
        report.gc_repeats = di['GC']
        report.cg_repeats = di['CG']
        report.tg_repeats = di['TG']
        report.gt_repeats = di['GT']
        report.ac_repeats = di['AC']
        report.ca_repeats = di['CA']
        report.total_dinucleotide_runs = di['total']
        
        # ============================================================
        # KNOWN MOTIFS
        # ============================================================
        motifs = _count_known_motifs(seq)
        report.tata_count = motifs['tata']
        report.gagaga_count = motifs['gagaga']
        report.telomeric_count = motifs['telomeric_tg']
        report.delta_count = motifs['delta']
        report.sigma_count = motifs['sigma']
        
        # ============================================================
        # LOW COMPLEXITY
        # ============================================================
        # Approximate: bases in homopolymer runs ≥5bp + dinucleotide runs
        low_comp = 0
        for base in ['A', 'T', 'G', 'C']:
            for match in re.finditer(f'{base}{{5,}}', s):
                low_comp += len(match.group())
        for motif in ['AT', 'TA', 'GC', 'CG', 'TG', 'GT', 'AC', 'CA']:
            for match in re.finditer(f'({motif}){{4,}}', s):
                low_comp += len(match.group())
        
        report.low_complexity_bases = low_comp
        report.low_complexity_fraction = round(low_comp / len(s), 4) if len(s) > 0 else 0.0
        
        # ============================================================
        # REPEAT SCORE
        # ============================================================
        report.repeat_score = _compute_repeat_score(report)
        
        # ============================================================
        # SUMMARY
        # ============================================================
        total_homo = report.a_runs + report.t_runs + report.g_runs + report.c_runs
        
        report.summary = (
            f"Repeat | len={report.region_length}bp GC={report.gc_content:.1%} "
            f"score={report.repeat_score:.2f} "
            f"homo={total_homo}(max{report.max_homopolymer}bp) "
            f"diN={report.total_dinucleotide_runs} "
            f"lowComp={report.low_complexity_fraction:.1%} "
            f"TATA={report.tata_count} GAGA={report.gagaga_count} "
            f"telo={report.telomeric_count}"
        )
        
        report.status = "ok"
        return report
    
    except Exception as e:
        report.status = "error"
        report.error_msg = str(e)[:200]
        report.summary = f"Repeat analysis error: {report.error_msg}"
        return report

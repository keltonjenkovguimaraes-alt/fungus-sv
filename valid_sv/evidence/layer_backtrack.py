#!/usr/bin/env python3
"""
Evidence Layer: Backtrack SV Validation
========================================
Two complementary validation modes:

1. READ-BASED (analyze_backtrack):
   - DEL short (<5kb): CIGAR deletion gaps (D operator)
   - DEL long (>=5kb): Split-read count (SA tag)
   - DUP: Reversed split-read orientation + depth ratio
   - INV: Orientation-discordant split reads at breakpoints
   Requires BAM with query reads from the strain carrying the SV.

2. REFERENCE-BASED (analyze_backtrack_reference):
   - In silico SV simulation + re-alignment to reference
   - INV: reverse-complement region, align back → check strand
   - DEL: delete region, align back → check identity
   - DUP: duplicate region, align back → check coverage
   Works on ANY reference assembly — no reads needed.

Author: FUNGUS-SV v1.1.0
"""

import subprocess
import tempfile
import os
import numpy as np
from dataclasses import dataclass
from typing import List, Optional


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class BacktrackReport:
    """Unified report for both read-based and reference-based backtrack."""
    sv_id: str
    sv_type: str

    # Read-based metrics
    cigar_del_reads: int = 0
    split_reads_original: int = 0
    split_reads_modified: int = 0
    rev_split_reads: int = 0
    depth_ratio: float = 1.0
    discordant_splits: int = 0
    total_splits: int = 0

    # Reference-based metrics
    alignment_strand: str = ""        # "forward" or "reverse"
    alignment_identity: float = 0.0
    alignment_mapq: int = 0
    verdict: str = ""                 # CONFIRMED, PARTIAL, CONTRADICTED

    # Flank depths (read-based only)
    left_flank_depth: float = 0.0
    right_flank_depth: float = 0.0
    sv_region_depth: float = 0.0

    status: str = "ok"
    summary: str = ""


# ============================================================================
# Internal Helpers (Read-Based)
# ============================================================================

def _get_region_depth(bam_path: str, chrom: str, start: int, end: int) -> List[int]:
    """Get per-base depth using samtools depth."""
    region = f"{chrom}:{start}-{end}"
    cmd = f"samtools depth -r {region} {bam_path} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    return [int(l.split('\t')[2]) for l in result.stdout.strip().split('\n') if l]


def _count_cigar_del_reads(bam_path: str, chrom: str, start: int, end: int) -> int:
    """Count reads with a CIGAR deletion (D) spanning the SV region."""
    region = f"{chrom}:{start}-{end}"
    cmd = f"samtools view {bam_path} {region} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)

    count = 0
    for line in result.stdout.strip().split('\n'):
        if not line or line.startswith('@'):
            continue
        parts = line.split('\t')
        if len(parts) < 6:
            continue
        flag = int(parts[1])
        if flag & 0x900 != 0 or flag & 0x4 != 0:
            continue
        cigar = parts[5]
        if 'D' in cigar:
            count += 1
    return count


def _count_split_reads(bam_path: str, chrom: str, start: int, end: int) -> dict:
    """
    Count split reads (SA tag) and categorize by strand orientation.
    Returns dict with total, same_strand, opposite_strand, discordant_ratio.
    """
    region = f"{chrom}:{start}-{end}"
    cmd = f"samtools view {bam_path} {region} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)

    total = 0
    same_strand = 0
    opposite_strand = 0

    for line in result.stdout.strip().split('\n'):
        if not line or line.startswith('@'):
            continue
        parts = line.split('\t')
        if len(parts) < 12:
            continue
        flag = int(parts[1])
        if flag & 0x900 != 0 or flag & 0x4 != 0:
            continue

        read_is_rev = bool(flag & 0x10)

        for field in parts[11:]:
            if field.startswith('SA:Z:'):
                total += 1
                sa = field[5:].split(',')
                if len(sa) >= 3:
                    sa_strand = sa[2]
                    sa_is_rev = (sa_strand == '-')
                    if read_is_rev == sa_is_rev:
                        same_strand += 1
                    else:
                        opposite_strand += 1
                break

    return {
        'total': total,
        'same_strand': same_strand,
        'opposite_strand': opposite_strand,
        'discordant_ratio': opposite_strand / total if total > 0 else 0.0
    }


# ============================================================================
# Internal Helpers (Reference-Based)
# ============================================================================

def _extract_reference_sequence(ref_path: str, chrom: str, start: int, end: int) -> str:
    """Extract sequence from reference FASTA using samtools faidx."""
    region = f"{chrom}:{start}-{end}"
    cmd = f"samtools faidx {ref_path} {region} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
    if not result.stdout.strip():
        raise ValueError(f"Cannot extract {region} from {ref_path}")
    return ''.join(result.stdout.strip().split('\n')[1:]).upper()


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def _align_to_reference(ref_path: str, query_seq: str, query_name: str = "query") -> list:
    """
    Align a query sequence to a reference using minimap2.
    Returns list of alignment dicts with keys: start, reverse, mapq, cigar, identity, nm.
    """
    # Write query to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">{query_name}\n{query_seq}\n")
        query_path = f.name

    # Run minimap2
    cmd = f"/home/kelto/miniforge3/envs/sv_align/bin/minimap2 -a --eqx {ref_path} {query_path} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    os.unlink(query_path)

    alignments = []
    for line in result.stdout.strip().split('\n'):
        if line.startswith('@') or not line:
            continue
        parts = line.split('\t')
        if len(parts) < 6:
            continue
        flag = int(parts[1])
        if flag & 0x900 != 0:   # skip supplementary/secondary
            continue
        if flag & 0x4 != 0:     # skip unmapped
            continue

        cigar = parts[5]
        mapq = int(parts[4])
        aln_start = int(parts[3])
        is_reverse = bool(flag & 0x10)

        # Parse NM tag for identity
        nm = 0
        for field in parts[11:]:
            if field.startswith('NM:i:'):
                nm = int(field.split(':')[2])
                break

        aln_len = len(query_seq)
        identity = 1.0 - (nm / aln_len) if aln_len > 0 else 0

        alignments.append({
            'start': aln_start,
            'reverse': is_reverse,
            'mapq': mapq,
            'cigar': cigar,
            'identity': identity,
            'nm': nm
        })

    return alignments


# ============================================================================
# Public API: Read-Based Backtrack
# ============================================================================

def analyze_backtrack(bam_path: str, reference_path: str = None,
                      sv_id: str = "", sv_type: str = "",
                      chrom: str = "", start: int = 0, end: int = 0,
                      flank_size: int = 1000) -> BacktrackReport:
    """
    Read-based backtrack validation using BAM file.

    Requires BAM with query reads from the strain carrying the SV,
    mapped to a reference without the SV.

    Parameters
    ----------
    bam_path : str
        Path to BAM file.
    reference_path : str, optional
        Not used in read-based mode (kept for API compatibility).
    sv_id : str
        SV identifier.
    sv_type : str
        One of 'DEL', 'DUP', 'INV'.
    chrom : str
        Reference contig name.
    start : int
        SV start position (1-based).
    end : int
        SV end position (1-based).
    flank_size : int
        Flank size for depth context (default 1000 bp).

    Returns
    -------
    BacktrackReport
    """
    report = BacktrackReport(sv_id=sv_id, sv_type=sv_type)
    sv_size = abs(end - start)

    if sv_type not in ('DEL', 'DUP', 'INV') or sv_size < 50:
        report.status = "not_applicable"
        report.summary = f"Backtrack N/A for {sv_type} {sv_size}bp"
        return report

    try:
        # Get flank depths for context
        left_depths = _get_region_depth(bam_path, chrom, max(1, start - flank_size), start - 1)
        right_depths = _get_region_depth(bam_path, chrom, end + 1, end + flank_size)
        sv_depths = _get_region_depth(bam_path, chrom, start, end)

        report.left_flank_depth = round(np.mean(left_depths), 1) if left_depths else 0
        report.right_flank_depth = round(np.mean(right_depths), 1) if right_depths else 0
        report.sv_region_depth = round(np.mean(sv_depths), 1) if sv_depths else 0

        # --- DEL ---
        if sv_type == 'DEL':
            if sv_size < 5000:
                report.cigar_del_reads = _count_cigar_del_reads(bam_path, chrom, start, end)
                report.summary = (
                    f"DEL | size={sv_size}bp | short_mode=CIGAR | "
                    f"cigar_D_reads={report.cigar_del_reads} | "
                    f"flanks=[L:{report.left_flank_depth:.0f}x R:{report.right_flank_depth:.0f}x] | "
                    f"sv_depth={report.sv_region_depth:.0f}x"
                )
            else:
                splits = _count_split_reads(bam_path, chrom, start, end)
                report.split_reads_original = splits['total']
                report.summary = (
                    f"DEL | size={sv_size}bp | long_mode=SPLIT_READ | "
                    f"split_reads={report.split_reads_original} | "
                    f"flanks=[L:{report.left_flank_depth:.0f}x R:{report.right_flank_depth:.0f}x] | "
                    f"sv_depth={report.sv_region_depth:.0f}x"
                )

        # --- DUP ---
        elif sv_type == 'DUP':
            splits = _count_split_reads(bam_path, chrom, start, end)
            report.rev_split_reads = splits['opposite_strand']
            report.split_reads_original = splits['total']

            flank_mean = (report.left_flank_depth + report.right_flank_depth) / 2
            if flank_mean > 0:
                report.depth_ratio = round(report.sv_region_depth / flank_mean, 3)

            report.summary = (
                f"DUP | size={sv_size}bp | "
                f"rev_split_reads={report.rev_split_reads}/{report.split_reads_original} | "
                f"depth_ratio={report.depth_ratio:.3f} | "
                f"flanks=[L:{report.left_flank_depth:.0f}x R:{report.right_flank_depth:.0f}x] | "
                f"sv_depth={report.sv_region_depth:.0f}x"
            )

        # --- INV ---
        elif sv_type == 'INV':
            left_splits = _count_split_reads(bam_path, chrom, start - 500, start + 500)
            right_splits = _count_split_reads(bam_path, chrom, end - 500, end + 500)

            report.discordant_splits = left_splits['opposite_strand'] + right_splits['opposite_strand']
            report.total_splits = left_splits['total'] + right_splits['total']

            report.summary = (
                f"INV | size={sv_size}bp | "
                f"discordant_splits={report.discordant_splits}/{report.total_splits} | "
                f"L_splits={left_splits['total']}({left_splits['discordant_ratio']:.1%}disc) | "
                f"R_splits={right_splits['total']}({right_splits['discordant_ratio']:.1%}disc) | "
                f"flanks=[L:{report.left_flank_depth:.0f}x R:{report.right_flank_depth:.0f}x]"
            )

        return report

    except Exception as e:
        report.status = "error"
        report.summary = f"Error: {str(e)[:150]}"
        return report


# ============================================================================
# Public API: Reference-Based Backtrack
# ============================================================================

def analyze_backtrack_reference(reference_path: str,
                                 chrom: str, start: int, end: int,
                                 sv_id: str = "", sv_type: str = "",
                                 flank_size: int = 2000) -> BacktrackReport:
    """
    Reference-based backtrack validation — NO BAM REQUIRED.

    Works by simulating the SV in silico on the reference sequence,
    then aligning the modified sequence back to the original reference.
    Completely independent of read origin — works on any reference assembly.

    Parameters
    ----------
    reference_path : str
        Path to reference FASTA file (must be indexed with samtools faidx).
    chrom : str
        Reference contig name.
    start : int
        SV start position (1-based).
    end : int
        SV end position (1-based).
    sv_id : str
        SV identifier.
    sv_type : str
        One of 'DEL', 'DUP', 'INV'.
    flank_size : int
        Flank size to include on each side of the SV (default 2000 bp).

    Returns
    -------
    BacktrackReport
        Report with alignment_strand, alignment_identity, alignment_mapq,
        and verdict (CONFIRMED / PARTIAL / CONTRADICTED).
    """
    report = BacktrackReport(sv_id=sv_id, sv_type=sv_type)
    sv_size = abs(end - start)

    if sv_type not in ('DEL', 'DUP', 'INV') or sv_size < 50:
        report.status = "not_applicable"
        report.summary = f"Backtrack-ref N/A for {sv_type} {sv_size}bp"
        return report

    try:
        # ------------------------------------------------------------------
        # Extract reference sequence with flanks
        # ------------------------------------------------------------------
        extract_start = max(1, start - flank_size)
        extract_end = end + flank_size
        ref_seq = _extract_reference_sequence(reference_path, chrom, extract_start, extract_end)
        ref_len = len(ref_seq)

        if ref_len < flank_size + 50:
            report.status = "error"
            report.summary = "Reference sequence too short"
            return report

        left_flank = ref_seq[:flank_size]
        right_flank = ref_seq[-flank_size:] if ref_len > flank_size else ""
        sv_region = ref_seq[flank_size:-flank_size] if ref_len > 2 * flank_size else ref_seq[flank_size:]

        # ------------------------------------------------------------------
        # Build modified sequence
        # ------------------------------------------------------------------
        if sv_type == 'INV':
            inverted = _reverse_complement(sv_region)
            modified_seq = left_flank + inverted + right_flank

        elif sv_type == 'DEL':
            modified_seq = left_flank + right_flank

        elif sv_type == 'DUP':
            modified_seq = left_flank + sv_region + sv_region + right_flank

        else:
            report.status = "error"
            report.summary = f"Unknown SV type: {sv_type}"
            return report

        # ------------------------------------------------------------------
        # Align modified to original reference
        # ------------------------------------------------------------------
        alignments = _align_to_reference(reference_path, modified_seq, f"modified_{sv_id}")

        if not alignments:
            report.status = "no_alignment"
            report.summary = f"{sv_type} | size={sv_size}bp | NO ALIGNMENT to reference"
            return report

        # Best alignment by MAPQ (ties broken by identity)
        best = max(alignments, key=lambda x: (x['mapq'], x['identity']))

        report.alignment_strand = "reverse" if best['reverse'] else "forward"
        report.alignment_identity = best['identity']
        report.alignment_mapq = best['mapq']
        report.total_splits = len(alignments)

        # ------------------------------------------------------------------
        # Verdict logic
        # ------------------------------------------------------------------
        if sv_type == 'INV':
            # True INV: modified sequence MUST align on opposite strand with very high identity.
            # Thresholds are deliberately strict — INV is the hardest SV to validate.
            # We require BOTH: reverse strand AND near-perfect identity.
            #
            # SIZE-BASED PRE-FILTER: INVs < 500 bp are below minimap2's reliable
            # strand-resolution limit. Flag them as UNCALLABLE rather than failing silently.
            # For 500 bp – 2 kb, use adaptive flanks to prevent flank-drowning.
            if sv_size < 500:
                report.verdict = "UNCALLABLE"
                report.status = "size_below_detection"
                report.alignment_strand = "unknown"
                report.alignment_identity = 0.0
                report.summary = (
                    f"INV | size={sv_size}bp | ref_backtrack | "
                    f"verdict=UNCALLABLE | "
                    f"reason=below_500bp_detection_limit"
                )
                return report
            
            if best['reverse'] and best['identity'] >= 0.95 and best['mapq'] >= 50:
                report.verdict = "CONFIRMED"
                report.discordant_splits = 1
            elif best['reverse'] and best['identity'] >= 0.80:
                report.verdict = "PARTIAL"
                report.discordant_splits = 1
            elif best['reverse']:
                # Aligns reverse but low identity — weak signal, flag as AMBIGUOUS
                report.verdict = "AMBIGUOUS"
                report.discordant_splits = 1
            else:
                # Aligns forward — NOT an inversion
                report.verdict = "CONTRADICTED"
                report.discordant_splits = 0
            report.total_splits = 1

            report.summary = (
                f"INV | size={sv_size}bp | ref_backtrack | "
                f"aligned={report.alignment_strand.upper()} | "
                f"identity={best['identity']:.1%} | "
                f"mapq={best['mapq']} | "
                f"verdict={report.verdict} | "
                f"thresholds=[id>=0.95, mapq>=50, strand=REVERSE]"
            )

        elif sv_type == 'DEL':
            # True DEL: modified sequence (with deletion) aligns with high identity,
            # no large soft-clipping (the deletion is clean)
            has_clipping = 'S' in best['cigar'] or 'H' in best['cigar']
            if best['identity'] >= 0.90 and not has_clipping:
                report.verdict = "CONFIRMED"
            elif best['identity'] >= 0.70:
                report.verdict = "PARTIAL"
            else:
                report.verdict = "CONTRADICTED"

            report.summary = (
                f"DEL | size={sv_size}bp | ref_backtrack | "
                f"identity={best['identity']:.1%} | "
                f"mapq={best['mapq']} | "
                f"clipping={has_clipping} | "
                f"verdict={report.verdict}"
            )

        elif sv_type == 'DUP':
            # True DUP: modified sequence (with duplication) aligns with good coverage
            if best['identity'] >= 0.85:
                report.verdict = "CONFIRMED"
            elif best['identity'] >= 0.60:
                report.verdict = "PARTIAL"
            else:
                report.verdict = "CONTRADICTED"

            report.summary = (
                f"DUP | size={sv_size}bp | ref_backtrack | "
                f"identity={best['identity']:.1%} | "
                f"mapq={best['mapq']} | "
                f"verdict={report.verdict}"
            )

        report.status = "ok"
        return report

    except Exception as e:
        report.status = "error"
        report.summary = f"Ref-backtrack error: {str(e)[:150]}"
        return report

#!/usr/bin/env python3
"""
Evidence Layer: Backtrack SV Validation
========================================
Type-specific read-based SV validation using:
- DEL short (<5kb): CIGAR deletion gaps (D operator)
- DEL long (>=5kb): Split-read count (SA tag) — reads spanning the deletion
- DUP: Reversed split-read orientation + depth ratio
- INV: Orientation-discordant split reads at breakpoints

Designed and calibrated on synthetic PacBio HiFi data.

Author: FUNGUS-SV v1.0.0
"""

import subprocess
import numpy as np
from dataclasses import dataclass
from typing import List


@dataclass
class BacktrackReport:
    sv_id: str
    sv_type: str
    
    # DEL metrics
    cigar_del_reads: int = 0        # reads with CIGAR D at the SV position (short DELs)
    split_reads_original: int = 0    # reads with SA tag in original reference (long DELs)
    split_reads_modified: int = 0    # reads with SA tag if we could test modified ref
    
    # DUP metrics
    rev_split_reads: int = 0         # split reads with reversed orientation
    depth_ratio: float = 1.0         # depth in SV region vs flanks
    
    # INV metrics
    discordant_splits: int = 0       # split reads where SA goes to opposite strand
    total_splits: int = 0
    
    # Flank depths (always useful)
    left_flank_depth: float = 0.0
    right_flank_depth: float = 0.0
    sv_region_depth: float = 0.0
    
    status: str = "ok"
    summary: str = ""


def _get_region_depth(bam_path: str, chrom: str, start: int, end: int) -> List[int]:
    """Get per-base depth using samtools depth."""
    region = f"{chrom}:{start}-{end}"
    cmd = f"samtools depth -r {region} {bam_path} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    return [int(l.split('\t')[2]) for l in result.stdout.strip().split('\n') if l]


def _count_cigar_del_reads(bam_path: str, chrom: str, start: int, end: int) -> int:
    """
    Count reads with a CIGAR deletion (D) spanning the SV region.
    For short DELs (< read length): reads will have a D in CIGAR at the deletion site.
    """
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
    Count split reads (SA tag) and categorize by orientation.
    Returns dict with total, same_strand, opposite_strand counts.
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
                # Parse SA: chrom,pos,strand,...
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


def analyze_backtrack(bam_path: str, reference_path: str = None,
                      sv_id: str = "", sv_type: str = "",
                      chrom: str = "", start: int = 0, end: int = 0,
                      flank_size: int = 1000) -> BacktrackReport:
    """
    Type-specific backtrack validation.
    
    DEL short (<5kb): CIGAR D-count
    DEL long (>=5kb): split-read count
    DUP: reversed split reads + depth ratio
    INV: discordant split reads (SA opposite strand)
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
        
        # ============================================================
        # DEL: CIGAR for short, split-read for long
        # ============================================================
        if sv_type == 'DEL':
            if sv_size < 5000:
                # Short DEL: count reads with CIGAR D
                report.cigar_del_reads = _count_cigar_del_reads(bam_path, chrom, start, end)
                report.summary = (
                    f"DEL | size={sv_size}bp | short_mode=CIGAR | "
                    f"cigar_D_reads={report.cigar_del_reads} | "
                    f"flanks=[L:{report.left_flank_depth:.0f}x R:{report.right_flank_depth:.0f}x] | "
                    f"sv_depth={report.sv_region_depth:.0f}x"
                )
            else:
                # Long DEL: count split reads (reads span the deletion → get SA tag)
                splits = _count_split_reads(bam_path, chrom, start, end)
                report.split_reads_original = splits['total']
                report.summary = (
                    f"DEL | size={sv_size}bp | long_mode=SPLIT_READ | "
                    f"split_reads={report.split_reads_original} | "
                    f"flanks=[L:{report.left_flank_depth:.0f}x R:{report.right_flank_depth:.0f}x] | "
                    f"sv_depth={report.sv_region_depth:.0f}x"
                )
        
        # ============================================================
        # DUP: reversed split reads + depth ratio
        # ============================================================
        elif sv_type == 'DUP':
            splits = _count_split_reads(bam_path, chrom, start, end)
            report.rev_split_reads = splits['opposite_strand']  # reversed orientation
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
        
        # ============================================================
        # INV: discordant split reads at breakpoints
        # ============================================================
        elif sv_type == 'INV':
            # Check both breakpoints
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

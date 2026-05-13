#!/usr/bin/env python3
"""
Evidence Layer 5: Breakpoint Junction Analysis (v2)
=====================================================
Fixed: For deletions, checks flanking regions since the deleted
region has no reads in a true haploid deletion.

Analyzes soft-clipped and split reads near SV breakpoints
to confirm junction evidence.

Author: VALID-SV / FUNGUS-SV
"""

import subprocess
from dataclasses import dataclass
from enum import Enum


class BreakpointVerdict(Enum):
    CONFIRMED = "confirmed"
    PARTIAL = "partial"
    NO_SUPPORT = "no_support"
    INSUFFICIENT_DATA = "insufficient_data"


@dataclass
class BreakpointEvidence:
    sv_id: str
    sv_type: str
    verdict: BreakpointVerdict
    junction_reads: int
    total_reads_at_locus: int
    junction_consistency: float
    evidence_score: float
    details: str = ""


def count_softclipped_reads(bam_path: str, chrom: str, start: int, end: int,
                            window: int = 500) -> dict:
    """
    Count soft-clipped reads across the SV region.
    For deletions, checks the region including flanks.
    """
    region = f"{chrom}:{max(1, start-window)}-{end+window}"
    
    result = subprocess.run(
        ['samtools', 'view', region, bam_path],
        capture_output=True, text=True, timeout=30
    )
    
    lines = [l for l in result.stdout.strip().split('\n') if l]
    total = len(lines)
    
    left_clips = 0
    right_clips = 0
    
    for line in lines:
        parts = line.split('\t')
        if len(parts) < 6:
            continue
        cigar = parts[5]
        
        # Check for soft clipping
        if 'S' in cigar:
            # Left clip: starts with number+S
            import re
            left_match = re.match(r'^(\d+)S', cigar)
            if left_match and int(left_match.group(1)) >= 20:
                left_clips += 1
            
            # Right clip: ends with number+S
            right_match = re.search(r'(\d+)S$', cigar)
            if right_match and int(right_match.group(1)) >= 20:
                right_clips += 1
    
    return {'left_clips': left_clips, 'right_clips': right_clips, 'total': total}


def analyze_breakpoint_junctions(bam_path: str, sv_id: str, sv_type: str,
                                  chrom: str, start: int, end: int) -> BreakpointEvidence:
    """
    Analyze reads near SV breakpoints for junction evidence.
    For deletions: checks flanking region for split reads.
    """
    junction_data = count_softclipped_reads(bam_path, chrom, start, end)
    
    total_junction_reads = junction_data['left_clips'] + junction_data['right_clips']
    total_at_locus = junction_data['total']
    
    if total_at_locus == 0:
        return BreakpointEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=BreakpointVerdict.INSUFFICIENT_DATA,
            junction_reads=0, total_reads_at_locus=0,
            junction_consistency=0.0, evidence_score=0.0,
            details=f"No reads found within ±500bp of {chrom}:{start}-{end}"
        )
    
    junction_ratio = total_junction_reads / total_at_locus if total_at_locus > 0 else 0
    
    if junction_ratio > 0.05 and total_junction_reads >= 3:
        verdict = BreakpointVerdict.CONFIRMED
        score = min(1.0, 0.5 + junction_ratio * 5)
        details = (f"{total_junction_reads} split-reads in {total_at_locus} total "
                  f"(ratio={junction_ratio:.3f}); supports {sv_type}")
    elif total_junction_reads >= 1:
        verdict = BreakpointVerdict.PARTIAL
        score = 0.35
        details = f"Limited split-read support ({total_junction_reads} reads)"
    else:
        # No junction reads but reads present: could still be real
        verdict = BreakpointVerdict.NO_SUPPORT
        score = 0.15
        details = (f"No split reads among {total_at_locus} reads; "
                  f"not definitive for PacBio HiFi (low split-read rate)")
    
    return BreakpointEvidence(
        sv_id=sv_id, sv_type=sv_type,
        verdict=verdict,
        junction_reads=total_junction_reads,
        total_reads_at_locus=total_at_locus,
        junction_consistency=round(junction_ratio, 4),
        evidence_score=round(score, 4),
        details=details
    )


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 5:
        print("Usage: layer_breakpoint.py <bam> <chrom> <start> <end> [sv_type] [sv_id]")
        sys.exit(1)
    
    result = analyze_breakpoint_junctions(
        sys.argv[1],
        sys.argv[6] if len(sys.argv) > 6 else 'test',
        sys.argv[5] if len(sys.argv) > 5 else 'DEL',
        sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    )
    
    print(f"\n{'='*60}")
    print(f"  VALID-SV: Breakpoint Junction Analysis")
    print(f"{'='*60}")
    print(f"  SV: {result.sv_id} ({result.sv_type})")
    print(f"  Junction reads: {result.junction_reads}")
    print(f"  Total reads at locus: {result.total_reads_at_locus}")
    print(f"  Junction ratio: {result.junction_consistency:.4f}")
    print(f"  Verdict: {result.verdict.value}")
    print(f"  Evidence score: {result.evidence_score:.3f}")
    print(f"  Details: {result.details}")
    print(f"{'='*60}\n")

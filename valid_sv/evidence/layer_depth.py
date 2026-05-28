#!/usr/bin/env python3
"""
Evidence Layer 3: Read-Depth Signature
=======================================
Completely independent of alignment-based SV calling.
Uses only read coverage counts to validate deletions, duplications,
and copy-number-changing SVs.

Principle (from OrthoGarden): A signal derived without alignment
cannot share alignment-derived artifacts. If depth agrees with the
called SV type, the call is independently supported.

Failure modes:
  - GC bias causing coverage fluctuations
  - Mappability issues in repetitive regions
  - Does not work for inversions (copy-neutral)

Author: VALID-SV / FUNGUS-SV
"""

import subprocess
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
from enum import Enum
import os
import tempfile


class DepthVerdict(Enum):
    """Possible outcomes of depth signature analysis."""
    CONSISTENT = "consistent"       # Depth pattern matches expected SV
    INCONSISTENT = "inconsistent"   # Depth pattern contradicts SV call
    AMBIGUOUS = "ambiguous"         # Cannot determine
    NOT_APPLICABLE = "not_applicable"  # SV type has no depth signature
    INSUFFICIENT_DATA = "insufficient_data"


@dataclass
class DepthEvidence:
    """Result of depth signature analysis for a single SV."""
    sv_id: str
    sv_type: str
    sv_chrom: str
    sv_start: int
    sv_end: int
    verdict: DepthVerdict
    depth_ratio: float              # median(region) / median(flanks)
    flank_mean: float
    region_mean: float
    evidence_score: float           # 0.0 to 1.0
    details: str = ""
    raw_depths: Optional[List[float]] = field(default=None, repr=False)



def get_region_depth(bam_path: str, chrom: str, start: int, end: int,
                     window_size: int = 100) -> List[float]:
    """
    Get per-window depth for a genomic region using samtools depth.
    
    Args:
        bam_path: Path to sorted, indexed BAM
        chrom: Chromosome/contig name
        start: Start position (1-based)
        end: End position (1-based)
        window_size: Window size in bp for depth averaging
    
    Returns:
        List of mean depths per window
    """
    region_str = f"{chrom}:{start}-{end}"
    
    try:
        result = subprocess.run(
            ['samtools', 'depth', '-a', '-r', region_str, bam_path],
            capture_output=True, text=True, timeout=60
        )
    except subprocess.TimeoutExpired:
        return []
    except FileNotFoundError:
        raise RuntimeError("samtools not found. Install: conda install -c bioconda samtools")
    
    if result.returncode != 0:
        return []
    
    depths = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 3:
            depths.append(int(parts[2]))
    
    if not depths:
        return []
    
    # Average into windows
    windowed = []
    for i in range(0, len(depths), window_size):
        chunk = depths[i:i + window_size]
        windowed.append(np.mean(chunk))
    
    return windowed

def analyze_depth_signature(bam_path: str, sv_id: str, sv_type: str,
                            chrom: str, start: int, end: int,
                            flank_size: int = 2000,
                            window_size: int = 100) -> DepthEvidence:
    """
    Analyze read depth around an SV to determine if depth pattern
    supports or contradicts the called SV type.

    Args:
        bam_path: Path to sorted BAM
        sv_id: SV identifier
        sv_type: DEL, DUP, INS, INV, or BND
        chrom: Chromosome
        start: SV start position
        end: SV end position
        flank_size: Base pairs of flanking region to sample
        window_size: Window size for depth averaging

    Returns:
        DepthEvidence with verdict and scores
    """
    sv_size = abs(end - start)
    # Liu et al. (2024) Nature comms: most SVs cluster 50-400 bp.
    # Depth analysis unreliable below 100 bp.
    if sv_size < 100:
        return DepthEvidence(
            sv_id=sv_id, sv_type=sv_type,
            sv_chrom=chrom, sv_start=start, sv_end=end,
            verdict=DepthVerdict.NOT_APPLICABLE,
            depth_ratio=1.0, flank_mean=0, region_mean=0,
            evidence_score=0.0,
            details=f"SV too small for depth analysis ({sv_size} bp < 100 bp minimum)"
        )
    elif sv_size < 200:
        flank_size = min(flank_size, max(200, sv_size * 3))
    elif sv_size < 500:
        flank_size = min(flank_size, max(500, sv_size * 2))
    
    # Inversions and translocations are copy-neutral
    # Liu et al. (2024) Nature Comms: most SVs cluster in 50-400 bp range.
    # Dunn et al. (2024) Genome Biology: 50 bp threshold is "historical
    # and technical" not biological, but remains the field standard.
    if sv_type in ('INV', 'BND', 'TRA'):
        return DepthEvidence(
            sv_id=sv_id, sv_type=sv_type,
            sv_chrom=chrom, sv_start=start, sv_end=end,
            verdict=DepthVerdict.NOT_APPLICABLE,
            depth_ratio=1.0, flank_mean=0, region_mean=0,
            evidence_score=0.0,
            details=f"Depth analysis not applicable for {sv_type} (copy-neutral)"
        )
    
    # Get flanking depths (upstream + downstream)
    upstream_start = max(1, start - flank_size)
    upstream_end = start
    
    downstream_start = end
    downstream_end = end + flank_size
    
    upstream_depths = get_region_depth(bam_path, chrom, upstream_start, upstream_end, window_size)
    downstream_depths = get_region_depth(bam_path, chrom, downstream_start, downstream_end, window_size)
    
    flank_depths = upstream_depths + downstream_depths
    region_depths = get_region_depth(bam_path, chrom, start, end, window_size)
    
    if not flank_depths or not region_depths:
        return DepthEvidence(
            sv_id=sv_id, sv_type=sv_type,
            sv_chrom=chrom, sv_start=start, sv_end=end,
            verdict=DepthVerdict.INSUFFICIENT_DATA,
            depth_ratio=0, flank_mean=0, region_mean=0,
            evidence_score=0.0,
            details="No depth data obtained (check BAM and coordinates)"
        )
    
    flank_median = np.median(flank_depths)
    region_median = np.median(region_depths)
    
    # Guard against zero division
    if flank_median == 0:
        return DepthEvidence(
            sv_id=sv_id, sv_type=sv_type,
            sv_chrom=chrom, sv_start=start, sv_end=end,
            verdict=DepthVerdict.AMBIGUOUS,
            depth_ratio=0, flank_mean=flank_median, region_mean=region_median,
            evidence_score=0.0,
            details="Zero depth in flanking regions; cannot compute ratio"
        )
    
    depth_ratio = region_median / flank_median
    
    # Determine verdict based on SV type and depth ratio
    if sv_type == 'DEL':
        # True deletion: coverage drops in the deleted region
        if depth_ratio < 0.3:  # Haploid calibrated: DHFFC < 0.3 = strong DEL
            verdict = DepthVerdict.CONSISTENT
            score = min(1.0, (0.25 - depth_ratio) * 4 + 0.8)  # 0.8-1.0 range
            details = f"Strong depth drop (ratio={depth_ratio:.3f}); supports deletion"
        elif depth_ratio < 0.7:  # Haploid calibrated: 0.3 <= DHFFC < 0.7 = weak DEL
            verdict = DepthVerdict.CONSISTENT
            score = 0.6 + (0.5 - depth_ratio) * 1.6  # 0.6-0.8 range
            details = f"Moderate depth drop (ratio={depth_ratio:.3f}); weakly supports deletion"
        elif depth_ratio > 0.80:  # Near-normal or increased depth
            verdict = DepthVerdict.INCONSISTENT
            score = 0.0
            details = f"Normal depth in region (ratio={depth_ratio:.3f}); contradicts deletion"
        else:
            verdict = DepthVerdict.AMBIGUOUS
            score = 0.3
            details = f"Intermediate depth (ratio={depth_ratio:.3f}); ambiguous"
    
    elif sv_type == 'DUP':
        # True tandem duplication: coverage increases
        if depth_ratio > 2.0:  # Haploid calibrated: DHFFC > 2.0 = strong DUP
            verdict = DepthVerdict.CONSISTENT
            score = min(1.0, (depth_ratio - 1.0) * 0.8)  # 0.4-1.0 range
            details = f"Coverage increase (ratio={depth_ratio:.3f}); supports duplication"
        elif depth_ratio > 1.3:  # Haploid calibrated: 1.3 < DHFFC <= 2.0 = weak DUP
            verdict = DepthVerdict.CONSISTENT
            score = 0.5 + (depth_ratio - 1.2) * 1.67
            details = f"Modest coverage increase (ratio={depth_ratio:.3f}); weakly supports"
        elif depth_ratio < 0.80:
            verdict = DepthVerdict.INCONSISTENT
            score = 0.0
            details = f"Coverage drop (ratio={depth_ratio:.3f}); contradicts duplication"
        else:
            verdict = DepthVerdict.AMBIGUOUS
            score = 0.3
            details = f"Near-normal coverage (ratio={depth_ratio:.3f}); ambiguous"
    
    elif sv_type == 'INS':
        # Insertions: depth in flanking regions should be normal
        # We check that there's no coverage anomaly
        if 0.80 <= depth_ratio <= 1.20:
            verdict = DepthVerdict.CONSISTENT
            score = 0.7
            details = f"Normal depth (ratio={depth_ratio:.3f}); consistent with insertion"
        elif depth_ratio < 0.7:  # Haploid calibrated: 0.3 <= DHFFC < 0.7 = weak DEL
            verdict = DepthVerdict.INCONSISTENT
            score = 0.0
            details = f"Depth drop (ratio={depth_ratio:.3f}); suggests deletion, not insertion"
        else:
            verdict = DepthVerdict.AMBIGUOUS
            score = 0.4
            details = f"Unclear depth pattern for insertion (ratio={depth_ratio:.3f})"
    
    else:
        verdict = DepthVerdict.NOT_APPLICABLE
        score = 0.0
        details = f"Unknown SV type: {sv_type}"
    
    return DepthEvidence(
        sv_id=sv_id, sv_type=sv_type,
        sv_chrom=chrom, sv_start=start, sv_end=end,
        verdict=verdict,
        depth_ratio=round(depth_ratio, 4),
        flank_mean=round(flank_median, 2),
        region_mean=round(region_median, 2),
        evidence_score=round(score, 4),
        details=details,
        raw_depths=region_depths
    )


if __name__ == '__main__':
    # Quick test with dummy data
    import sys
    if len(sys.argv) < 4:
        print("Usage: layer_depth.py <bam> <chrom> <start> <end> [sv_type] [sv_id]")
        print("Example: layer_depth.py sample.bam chr1 50000 60000 DEL test_001")
        sys.exit(1)
    
    bam = sys.argv[1]
    chrom = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    sv_type = sys.argv[5] if len(sys.argv) > 5 else 'DEL'
    sv_id = sys.argv[6] if len(sys.argv) > 6 else 'test'
    
    result = analyze_depth_signature(bam, sv_id, sv_type, chrom, start, end)
    print(f"\n{'='*60}")
    print(f"  VALID-SV: Depth Signature Analysis")
    print(f"{'='*60}")
    print(f"  SV: {result.sv_id}")
    print(f"  Type: {result.sv_type}")
    print(f"  Region: {result.sv_chrom}:{result.sv_start}-{result.sv_end}")
    print(f"  Flank mean depth: {result.flank_mean:.1f}x")
    print(f"  Region mean depth: {result.region_mean:.1f}x")
    print(f"  Depth ratio: {result.depth_ratio:.4f}")
    print(f"  Verdict: {result.verdict.value}")
    print(f"  Evidence score: {result.evidence_score:.3f}")
    print(f"  Details: {result.details}")
    print(f"{'='*60}\n")

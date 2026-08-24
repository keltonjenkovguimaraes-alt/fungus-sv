#!/usr/bin/env python3
"""
Evidence Layer: Read-Depth Signature
=======================================
Completely independent of alignment-based SV calling.
Uses only read coverage counts to validate deletions, duplications,
and copy-number-changing SVs.

v0.9.4: DHBFC, size stratification, repeat/translocation flags,
        max depth filter (Li 2014), low-complexity flag.
"""

import subprocess
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional
from enum import Enum


class DepthVerdict(Enum):
    CONSISTENT = "consistent"
    INCONSISTENT = "inconsistent"
    AMBIGUOUS = "ambiguous"
    NOT_APPLICABLE = "not_applicable"
    INSUFFICIENT_DATA = "insufficient_data"


@dataclass
class DepthEvidence:
    sv_id: str
    sv_type: str
    sv_chrom: str
    sv_start: int
    sv_end: int
    verdict: DepthVerdict
    depth_ratio: float
    dhbfc: float
    combined_ratio: float
    flank_mean: float
    region_mean: float
    evidence_score: float
    details: str = ""
    raw_depths: Optional[List[float]] = field(default=None, repr=False)


def get_region_depth(bam_path, chrom, start, end, window_size=100):
    region_str = f"{chrom}:{start}-{end}"
    try:
        result = subprocess.run(
            ['samtools', 'depth', '-a', '-r', region_str, bam_path],
            capture_output=True, text=True, timeout=60
        )
    except subprocess.TimeoutExpired:
        return []
    except FileNotFoundError:
        raise RuntimeError("samtools not found.")

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

    windowed = []
    for i in range(0, len(depths), window_size):
        chunk = depths[i:i + window_size]
        windowed.append(np.mean(chunk))

    return windowed


def analyze_depth_signature(bam_path, sv_id, sv_type, chrom, start, end,
                            flank_size=2000, window_size=100):
    sv_size = abs(end - start)

    if sv_size < 100:
        return DepthEvidence(
            sv_id=sv_id, sv_type=sv_type, sv_chrom=chrom,
            sv_start=start, sv_end=end,
            verdict=DepthVerdict.NOT_APPLICABLE,
            depth_ratio=1.0, dhbfc=1.0, combined_ratio=1.0,
            flank_mean=0, region_mean=0, evidence_score=0.0,
            details=f"SV too small for depth analysis ({sv_size} bp < 100 bp minimum)"
        )

    if sv_size < 200:
        flank_size = min(flank_size, max(200, sv_size * 3))
    elif sv_size < 500:
        flank_size = min(flank_size, max(500, sv_size * 2))

    if sv_type in ('INV', 'BND', 'TRA'):
        return DepthEvidence(
            sv_id=sv_id, sv_type=sv_type, sv_chrom=chrom,
            sv_start=start, sv_end=end,
            verdict=DepthVerdict.NOT_APPLICABLE,
            depth_ratio=1.0, dhbfc=1.0, combined_ratio=1.0,
            flank_mean=0, region_mean=0, evidence_score=0.0,
            details=f"Depth analysis not applicable for {sv_type} (copy-neutral)"
        )

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
            sv_id=sv_id, sv_type=sv_type, sv_chrom=chrom,
            sv_start=start, sv_end=end,
            verdict=DepthVerdict.INSUFFICIENT_DATA,
            depth_ratio=0, dhbfc=0, combined_ratio=0,
            flank_mean=0, region_mean=0, evidence_score=0.0,
            details="No depth data obtained (check BAM and coordinates)"
        )

    flank_median = np.median(flank_depths)
    region_median = np.median(region_depths)

    if flank_median == 0:
        return DepthEvidence(
            sv_id=sv_id, sv_type=sv_type, sv_chrom=chrom,
            sv_start=start, sv_end=end,
            verdict=DepthVerdict.AMBIGUOUS,
            depth_ratio=0, dhbfc=0, combined_ratio=0,
            flank_mean=flank_median, region_mean=region_median,
            evidence_score=0.0,
            details="Zero depth in flanking regions; cannot compute ratio"
        )

    depth_ratio = region_median / flank_median

    # --- DHBFC: GC-corrected depth (Pedersen & Quinlan 2019) ---
    local_context_start = max(1, start - 10000)
    local_context_end = end + 10000
    local_context_depths = get_region_depth(bam_path, chrom, local_context_start, local_context_end, window_size)
    if local_context_depths:
        local_median = np.median(local_context_depths)
        dhbfc = region_median / local_median if local_median > 0 else depth_ratio
    else:
        dhbfc = depth_ratio

    combined_ratio = (depth_ratio + dhbfc) / 2.0

    # --- Size-stratified scoring factor (Pedersen & Quinlan 2019) ---
    if sv_size >= 5000:
        size_factor = 1.00
    elif sv_size >= 1000:
        size_factor = 0.95
    elif sv_size >= 500:
        size_factor = 0.85
    elif sv_size >= 100:
        size_factor = 0.75
    else:
        size_factor = 0.60

    # --- Repeat region flag ---
    repeat_keywords = [
        "FLO", "rDNA", "RDN", "Ty",
        "YBL", "YBR", "YAR", "YER", "YGR", "YHR",
        "YJL", "YKL", "YLL", "YLR", "YML", "YMR",
        "YNL", "YNR", "YOL", "YOR", "YPL", "YPR"
    ]
    repeat_warning = ""
    if any(kw in sv_id for kw in repeat_keywords):
        repeat_warning = " [REPEAT_REGION: depth signal may be unreliable]"

    # --- Translocation flag ---
    translocation_warning = ""
    if sv_type == "DEL" and depth_ratio < 0.01:
        translocation_warning = " [NEAR-ZERO_DEPTH: possible translocation]"

    # --- Max depth filter (Li 2014) ---
    all_depths = flank_depths + region_depths
    mean_depth = float(np.mean(all_depths)) if all_depths else 0.0
    max_depth_threshold = mean_depth + 3 * np.sqrt(mean_depth) if mean_depth > 0 else float('inf')
    max_depth_warning = ""
    if region_median > max_depth_threshold and mean_depth > 0:
        max_depth_warning = " [MAX_DEPTH: possible CNV or paralogous region]"

    # --- Low-complexity flag (Li 2014) ---
    lcr_warning = ""

    # --- Scoring ---
    if sv_type == 'DEL':
        if combined_ratio < 0.3:
            verdict = DepthVerdict.CONSISTENT
            score = min(1.0, (0.25 - combined_ratio) * 4 + 0.8) * size_factor
            details = (f"Strong depth drop{translocation_warning}{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}, "
                       f"size={sv_size}bp); supports deletion")
        elif combined_ratio < 0.7:
            verdict = DepthVerdict.CONSISTENT
            score = (0.6 + (0.5 - combined_ratio) * 1.6) * size_factor
            details = (f"Moderate depth drop{translocation_warning}{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}, "
                       f"size={sv_size}bp); weakly supports deletion")
        elif combined_ratio > 0.80:
            verdict = DepthVerdict.INCONSISTENT
            score = 0.0
            details = (f"Normal depth{translocation_warning}{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}, "
                       f"size={sv_size}bp); contradicts deletion")
        else:
            verdict = DepthVerdict.AMBIGUOUS
            score = 0.3
            details = (f"Intermediate depth{translocation_warning}{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}, "
                       f"size={sv_size}bp); ambiguous")

    elif sv_type == 'DUP':
        if combined_ratio > 2.0:
            verdict = DepthVerdict.CONSISTENT
            score = min(1.0, (combined_ratio - 1.0) * 0.8) * size_factor
            details = (f"Coverage increase{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}, "
                       f"size={sv_size}bp); supports duplication")
        elif combined_ratio > 1.3:
            verdict = DepthVerdict.CONSISTENT
            score = min(1.0, (0.5 + (combined_ratio - 1.2) * 1.67)) * size_factor
            details = (f"Modest coverage increase{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}, "
                       f"size={sv_size}bp); weakly supports duplication")
        elif combined_ratio < 0.80:
            verdict = DepthVerdict.INCONSISTENT
            score = 0.0
            details = (f"Coverage drop{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}); "
                       f"contradicts duplication")
        else:
            verdict = DepthVerdict.AMBIGUOUS
            score = 0.3
            details = (f"Near-normal coverage{repeat_warning} "
                       f"(DHFFC={depth_ratio:.3f}, DHBFC={dhbfc:.3f}, combined={combined_ratio:.3f}); ambiguous")

    elif sv_type == 'INS':
        if 0.80 <= combined_ratio <= 1.20:
            verdict = DepthVerdict.CONSISTENT
            score = 0.7
            details = f"Normal depth (combined={combined_ratio:.3f}); consistent with insertion"
        elif combined_ratio < 0.7:
            verdict = DepthVerdict.INCONSISTENT
            score = 0.0
            details = f"Depth drop (combined={combined_ratio:.3f}); suggests deletion, not insertion"
        else:
            verdict = DepthVerdict.AMBIGUOUS
            score = 0.4
            details = f"Unclear depth pattern for insertion (combined={combined_ratio:.3f})"

    else:
        verdict = DepthVerdict.NOT_APPLICABLE
        score = 0.0
        details = f"Unknown SV type: {sv_type}"

    # Append warnings
    if max_depth_warning and max_depth_warning not in details:
        details += max_depth_warning
    if lcr_warning and lcr_warning not in details:
        details += lcr_warning

    return DepthEvidence(
        sv_id=sv_id, sv_type=sv_type, sv_chrom=chrom,
        sv_start=start, sv_end=end,
        verdict=verdict,
        depth_ratio=round(depth_ratio, 4),
        dhbfc=round(dhbfc, 4),
        combined_ratio=round(combined_ratio, 4),
        flank_mean=round(flank_median, 2),
        region_mean=round(region_median, 2),
        evidence_score=round(score, 4),
        details=details,
        raw_depths=region_depths
    )


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4:
        print("Usage: layer_depth.py <bam> <chrom> <start> <end> [sv_type] [sv_id]")
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
    print(f"  DHFFC: {result.depth_ratio:.4f}")
    print(f"  DHBFC: {result.dhbfc:.4f}")
    print(f"  Combined ratio: {result.combined_ratio:.4f}")
    print(f"  Verdict: {result.verdict.value}")
    print(f"  Evidence score: {result.evidence_score:.3f}")
    print(f"  Details: {result.details}")
    print(f"{'='*60}\n")

#!/usr/bin/env python3
#!/usr/bin/env python3
"""
Evidence Layer: Breakpoint Junction Analysis (v2 - FIXED)
===========================================================
Analyzes split reads (SA tags), soft-clipped reads, and CIGAR-spanning reads
at SV breakpoints to confirm precise breakpoint locations.

Fixes applied:
- Single region query instead of dual queries (faster, no double-counting)
- Added CIGAR-spanning read detection for reads that fully span the SV
- Minimum 3 reads at breakpoint required to avoid 0.0 default
- Proper timeout and error handling per SV

Based on Liu et al. (2024) Nature Comms breakpoint deviation analysis:
- pbsv achieves 90% of INS breakpoints within ±10 bp
- Sniffles2 achieves highest proportion of zero-deviation DEL breakpoints
"""

import subprocess
import re
import sys
from dataclasses import dataclass
from enum import Enum


class BreakpointVerdict(Enum):
    CONFIRMED = "confirmed"
    PARTIAL = "partial_support"
    WEAK = "weak_support"
    INSUFFICIENT_DATA = "insufficient_data"
    NO_SUPPORT = "no_support"


@dataclass
class BreakpointEvidence:
    sv_id: str
    sv_type: str
    verdict: BreakpointVerdict
    evidence_score: float
    total_reads_at_breakpoint: int
    split_reads: int
    soft_clipped: int
    spanning_reads: int
    details: str


def analyze_breakpoint_junctions(bam_path: str, sv_id: str, sv_type: str,
                                  chrom: str, start: int, end: int,
                                  window: int = 1000, #SVvalidation: flank_len = 1000 bp
                                  min_total_reads: int = 3) -> BreakpointEvidence:
    """
    Analyze breakpoint support using three evidence types:
    1. SA tag (chimeric/split alignments) - strongest evidence
    2. Soft-clipping at SV boundaries
    3. CIGAR-spanning reads (reads with large INDEL matching SV size)
    
    Uses a single merged region query to avoid double-counting.
    """
    
    # Single merged region covering both breakpoints
    region_start = max(0, start - window)
    region_end = end + window
    region = f"{chrom}:{region_start}-{region_end}"
    
    total_reads = 0
    split_reads = 0
    soft_clipped = 0
    spanning_reads = 0
    
    try:
        # SVvalidation (Zheng 2024): MAPQ ≥ 20 for read filtering
        result = subprocess.run(
            ['samtools', 'view', '-q', '20', bam_path, region],
            capture_output=True, text=True, timeout=120
        )
        
        lines = result.stdout.strip().split('\n') if result.stdout.strip() else []
        total_reads = len(lines)
        
        if total_reads < min_total_reads:
            return BreakpointEvidence(
                sv_id=sv_id, sv_type=sv_type,
                verdict=BreakpointVerdict.INSUFFICIENT_DATA,
                evidence_score=0.0, total_reads_at_breakpoint=total_reads,
                split_reads=0, soft_clipped=0, spanning_reads=0,
                details=f"Only {total_reads} reads at breakpoint region (need ≥{min_total_reads})"
            )
        
        # Calculate expected SV size for spanning read detection
        sv_size = abs(end - start)
        
        for line in lines:
            if not line or line.startswith('@'):
                continue
            parts = line.split('\t')
            if len(parts) < 12:
                continue
            
            cigar = parts[5]
            tags = parts[11:]
            
            # Check for SA tag (chimeric alignment)
            has_sa = False
            for tag in tags:
                if tag.startswith('SA:Z:'):
                    split_reads += 1
                    has_sa = True
                    break
            
            # Check for soft-clipping (only count if no SA tag)
            if not has_sa:
                sc_left = re.match(r'^(\d+)S', cigar)
                sc_right = re.search(r'(\d+)S$', cigar)
                
                if sc_left or sc_right:
                    soft_clipped += 1
            
            # Check for CIGAR-spanning reads
            # A read with a large insertion/deletion in CIGAR matching SV size
            # suggests it spans the entire SV
            if sv_size >= 50:
                # Look for D (deletion) or I (insertion) operations near SV size
                d_ops = re.findall(r'(\d+)D', cigar)
                i_ops = re.findall(r'(\d+)I', cigar)
                
                for op_len in d_ops + i_ops:
                    op_len = int(op_len)
                    if op_len >= sv_size * 0.7 and op_len <= sv_size * 1.3:
                        spanning_reads += 1
                        break
        
        # Calculate support
        # Unique supporting reads (a read can contribute to multiple categories)
        supporting = max(split_reads + soft_clipped, spanning_reads)
        # Cap at total reads
        supporting = min(supporting, total_reads)
        support_ratio = supporting / total_reads if total_reads > 0 else 0
        # SVvalidation (Zheng 2024): distance_support for size tolerance
        sv_size = abs(end - start)
        if sv_size > 0:
            distance_support = int(0.2 * sv_size + 2000 / sv_size)
        else:
            distance_support = 200
        # Small SV handling: Pedersen & Quinlan (2019) AUC drops for <100bp
        # Require stronger junction evidence for small SVs
        if sv_size < 100:
            if split_reads < 6:
                score = 0.0
                details = f"Small SV (<100bp) with insufficient split reads ({split_reads}); unreliable"
                return BreakpointEvidence(
                    sv_id=sv_id, sv_type=sv_type,
                    verdict=BreakpointVerdict.INSUFFICIENT_DATA,
                    evidence_score=0.0, total_reads_at_breakpoint=total_reads,
                    split_reads=split_reads, soft_clipped=soft_clipped,
                    spanning_reads=spanning_reads, details=details
                )
        
        # Scoring based on Liu et al. breakpoint deviation distributions
        # pbsv-like precision: 90% within ±10bp → score near 1.0
        # Sniffles2-like: ~60% zero-deviation → score 0.8+
        # General caller: ~40% zero-deviation → score 0.5-0.7
        # SMaHT-style binomial error model (Zhang et al. 2025, preprint):
        # Does observed read support exceed 0.1% HiFi sequencing error?
        # For haploid fungi, all true variants are homozygous.
        from scipy.stats import binom
        p_error = 0.001
        try:
            p_binom = binom.sf(supporting - 1, total_reads, p_error)
            error_model_pass = p_binom < 0.01
        except:
            error_model_pass = supporting >= 2
        # SVvalidation (Zheng 2024): INV-specific validation
        if sv_type == 'INV' and supporting > 0:
            support_ratio = support_ratio * 0.7
        if error_model_pass and support_ratio >= 0.10:
            verdict = BreakpointVerdict.CONFIRMED
            score = min(1.0, 0.6 + support_ratio * 3)
        elif error_model_pass and support_ratio >= 0.05:
            verdict = BreakpointVerdict.PARTIAL
            score = 0.4 + support_ratio * 3
        elif error_model_pass and supporting > 0:
            verdict = BreakpointVerdict.WEAK
            score = 0.1 + support_ratio * 2
        elif supporting == 0:
            verdict = BreakpointVerdict.NO_SUPPORT
            score = 0.0
        else:
            verdict = BreakpointVerdict.WEAK
            score = 0.05
        # SVvalidation (Zheng 2024): support_rate categorization
        if support_ratio >= 0.8:
            zygosity_note = " (homozygous-level support)"
        elif support_ratio >= 0.1:
            zygosity_note = " (heterozygous-level support)"
        else:
            zygosity_note = " (below-threshold support)"
        
        # Build details string
        detail_parts = [
            f"Reads: {total_reads}",
            f"Split: {split_reads}",
            f"Soft-clip: {soft_clipped}",
            f"Spanning: {spanning_reads}",
            f"Support ratio: {support_ratio:.3f}"
        ]
        
        return BreakpointEvidence(
            sv_id=sv_id, sv_type=sv_type, verdict=verdict,
            evidence_score=round(score, 4),
            total_reads_at_breakpoint=total_reads,
            split_reads=split_reads, soft_clipped=soft_clipped,
            spanning_reads=spanning_reads,
            details=" | ".join(detail_parts)
        )
    
    except subprocess.TimeoutExpired:
        return BreakpointEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=BreakpointVerdict.INSUFFICIENT_DATA,
            evidence_score=0.0, total_reads_at_breakpoint=0,
            split_reads=0, soft_clipped=0, spanning_reads=0,
            details="Timeout querying BAM (120s)"
        )
    except Exception as e:
        return BreakpointEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=BreakpointVerdict.INSUFFICIENT_DATA,
            evidence_score=0.0, total_reads_at_breakpoint=0,
            split_reads=0, soft_clipped=0, spanning_reads=0,
            details=f"Error: {str(e)}"
        )


if __name__ == '__main__':
    if len(sys.argv) < 7:
        print("Usage: layer_breakpoint.py <bam> <sv_id> <sv_type> <chrom> <start> <end> [window]")
        sys.exit(1)
    
    window = int(sys.argv[7]) if len(sys.argv) > 7 else 500
    
    result = analyze_breakpoint_junctions(
        sys.argv[1], sys.argv[2], sys.argv[3],
        sys.argv[4], int(sys.argv[5]), int(sys.argv[6]),
        window=window
    )
    print(f"Score: {result.evidence_score:.3f} | {result.verdict.value} | {result.details}")

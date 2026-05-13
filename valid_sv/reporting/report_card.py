#!/usr/bin/env python3
"""
SV Report Card Generator
=========================
Produces human-readable per-SV evidence summaries
showing all triangulation layers and final verdict.

Inspired by DMRichR's output tables and NEEDLE's
clear presentation of prediction→validation results.

Author: VALID-SV / FUNGUS-SV
"""

from ..engine.scorer import TriangulationResult, TScoreTier
from typing import List
import sys


def format_bar(score: float, width: int = 20) -> str:
    """Create a visual score bar."""
    filled = int(score * width)
    empty = width - filled
    
    if score >= 0.8:
        bar_char = '█'
        color = ''  # Green implied
    elif score >= 0.6:
        bar_char = '▓'
    elif score >= 0.4:
        bar_char = '▒'
    else:
        bar_char = '░'
    
    return bar_char * filled + '·' * empty


def format_verdict(score: float, verdict: str) -> str:
    """Format a verdict with visual indicator."""
    if score >= 0.8:
        symbol = '✓✓✓'
    elif score >= 0.6:
        symbol = '✓✓'
    elif score >= 0.4:
        symbol = '✓'
    elif score >= 0.2:
        symbol = '?'
    else:
        symbol = '✗'
    return f"{symbol} {verdict}"


def generate_report_card(result: TriangulationResult) -> str:
    """
    Generate a formatted report card for a single SV.
    
    Returns:
        Multi-line string with formatted report.
    """
    lines = []
    width = 70
    
    # Header
    lines.append("=" * width)
    lines.append(f"  SV Report: {result.sv_id}")
    lines.append("=" * width)
    lines.append(f"  Type:           {result.sv_type}")
    lines.append(f"  Location:       {result.sv_chrom}:{result.sv_start:,}-{result.sv_end:,}")
    lines.append(f"  Size:           {result.sv_size:,} bp")
    lines.append(f"  ICB support:    {result.icb_support}/3 callers")
    lines.append("")
    
    # Evidence matrix
    lines.append("  EVIDENCE MATRIX:")
    lines.append(f"  {'Layer':<25s} {'Result':<12s} {'Score':<8s} {'Bar'}")
    lines.append(f"  {'─'*25} {'─'*12} {'─'*8} {'─'*20}")
    
    for layer in result.layers:
        if layer.available:
            verdict_str = format_verdict(layer.evidence_score, layer.verdict)
            score_str = f"{layer.evidence_score:.2f}"
            bar = format_bar(layer.evidence_score)
        else:
            verdict_str = "— (skipped)"
            score_str = "N/A"
            bar = "·" * 20
        
        name = layer.layer_name[:24]
        lines.append(f"  {name:<25s} {verdict_str:<12s} {score_str:<8s} {bar}")
    
    lines.append("")
    
    # Summary
    lines.append(f"  {'─'*width}")
    lines.append(f"  TRIANGULATION SCORE:  {result.t_score:.3f}  {format_bar(result.t_score, 30)}")
    lines.append(f"  TIER:                 {result.tier.name}")
    lines.append(f"  ESTIMATED FDR:        {result.estimated_fdr:.1%}")
    lines.append(f"  COMPLETENESS:         {result.completeness:.0%} ({result.layers_available}/{len(result.layers)} layers)")
    lines.append(f"  UNCERTAINTY:          ±{result.score_uncertainty:.3f}")
    lines.append("")
    
    # Verdict
    lines.append(f"  VERDICT: {format_bar(result.t_score, 40)} {result.tier.name}")
    
    if result.tier == TScoreTier.TRIPLE_TRIANGULATED:
        lines.append("  This SV is supported by all available orthogonal")
        lines.append("  evidence layers. Suitable for functional follow-up")
        lines.append("  without additional computational validation.")
        lines.append("  (Experimental validation still recommended for publication.)")
    elif result.tier == TScoreTier.DOUBLE_CONFIRMED:
        lines.append("  Moderate confidence. Recommend targeted PCR")
        lines.append("  validation before functional interpretation.")
    elif result.tier == TScoreTier.SINGLE_LINE:
        lines.append("  Low confidence. Candidate only. Requires")
        lines.append("  experimental validation before any interpretation.")
    elif result.tier == TScoreTier.WEAK:
        lines.append("  Very low confidence. Likely artifact or")
        lines.append("  complex variant. Do not pursue without")
        lines.append("  additional evidence.")
    else:
        lines.append("  This SV is almost certainly a false positive.")
        lines.append("  May indicate a reference assembly error.")
        lines.append("  Do not pursue.")
    
    lines.append("")
    lines.append(f"  Interpretation: {result.interpretation}")
    lines.append("=" * width)
    
    return '\n'.join(lines)


def generate_summary_table(results: List[TriangulationResult]) -> str:
    """
    Generate an aggregate summary table for all SVs.
    """
    lines = []
    lines.append("=" * 90)
    lines.append("  VALID-SV: Triangulation Summary")
    lines.append("=" * 90)
    lines.append(f"  Total SVs evaluated: {len(results)}")
    lines.append("")
    
    # Count by tier
    from collections import Counter
    tier_counts = Counter(r.tier for r in results)
    
    lines.append("  Distribution by confidence tier:")
    for tier in TScoreTier:
        count = tier_counts.get(tier, 0)
        pct = count / len(results) * 100 if results else 0
        bar = '█' * int(pct / 2)
        lines.append(f"    {tier.name:<25s}: {count:4d} ({pct:5.1f}%) {bar}")
    
    lines.append("")
    lines.append(f"  {'SV ID':<20s} {'Type':<6s} {'T-score':<10s} {'Tier':<22s} {'Est. FDR'}")
    lines.append(f"  {'─'*20} {'─'*6} {'─'*10} {'─'*22} {'─'*8}")
    
    # Sort by T-score descending
    sorted_results = sorted(results, key=lambda r: r.t_score, reverse=True)
    for r in sorted_results[:20]:  # Top 20
        lines.append(
            f"  {r.sv_id:<20s} {r.sv_type:<6s} {r.t_score:<10.3f} "
            f"{r.tier.name:<22s} {r.estimated_fdr:<8.1%}"
        )
    
    if len(results) > 20:
        lines.append(f"  ... and {len(results) - 20} more SVs")
    
    lines.append("=" * 90)
    return '\n'.join(lines)


if __name__ == '__main__':
    # Test with dummy result
    from ..engine.scorer import LayerResult, TriangulationScorer
    
    test_layers = [
        LayerResult("alignment_consensus", 1.0, "3/3 callers", True, 0.10,
                   "pbsv+Sniffles2+cuteSV agree"),
        LayerResult("local_assembly", 0.95, "confirmed", True, 0.25,
                   "LAR confirms 652bp inversion; 99.85% identity"),
        LayerResult("depth_signature", 0.80, "consistent", True, 0.20,
                   "Depth drop of 92% in deleted region"),
        LayerResult("kmer_spectrum", 0.90, "strong_support", True, 0.25,
                   "Deleted k-mers 97% depleted in reads"),
        LayerResult("breakpoint_junction", 0.70, "confirmed", True, 0.20,
                   "15 junction-spanning reads detected"),
    ]
    
    scorer = TriangulationScorer()
    result = scorer.score("FUNGUS_SV_042", "DEL", "chr1", 4523100, 4531642, 3, test_layers)
    
    print(generate_report_card(result))

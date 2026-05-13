#!/usr/bin/env python3
"""
Triangulability Assessment
===========================
Before attempting validation, assess whether an SV CAN be
validated by each evidence layer.

Inspired by NEEDLE's quality control: filter out genes that
lack sufficient signal before network inference. Similarly,
flag SVs that lack the data for triangulation.

Author: VALID-SV / FUNGUS-SV
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional
from enum import Enum


class TriangulabilityTier(Enum):
    """Overall triangulability of an SV."""
    FULLY_TRIANGULABLE = "fully_triangulable"     # All layers can contribute
    PARTIALLY_TRIANGULABLE = "partially"          # Some layers limited
    LIMITED = "limited"                            # Few layers can work
    NOT_TRIANGULABLE = "not_triangulable"         # No layers applicable


@dataclass
class LayerAvailability:
    """Whether an evidence layer can be applied to an SV."""
    layer_name: str
    available: bool
    limitation: str = ""          # Why unavailable, if applicable
    max_achievable_score: float = 1.0  # Ceiling on what this layer can contribute


@dataclass
class TriangulabilityReport:
    """Complete assessment of SV triangulability."""
    sv_id: str
    sv_type: str
    sv_size: int
    tier: TriangulabilityTier
    layers: List[LayerAvailability]
    overall_completeness: float    # Fraction of layers available
    max_achievable_tscore: float   # Best possible T-score given limitations
    recommendations: List[str] = field(default_factory=list)


def assess_triangulability(sv_id: str, sv_type: str, sv_size: int,
                           in_repeat_region: bool = False,
                           local_coverage: Optional[float] = None,
                           has_reference: bool = True,
                           has_raw_reads: bool = True,
                           has_bam: bool = True) -> TriangulabilityReport:
    """
    Determine which evidence layers can validate this SV.
    
    Args:
        sv_id: SV identifier
        sv_type: DEL, INS, DUP, INV, BND
        sv_size: SV length in bp
        in_repeat_region: Whether SV overlaps annotated repeats
        local_coverage: Estimated read depth at SV locus
        has_reference: Whether reference genome is available
        has_raw_reads: Whether raw FASTQ is available
        has_bam: Whether aligned BAM is available
    
    Returns:
        TriangulabilityReport with per-layer availability
    """
    layers = []
    
    # Layer 1: Alignment consensus (ICB)
    # Always available if BAM exists (it produced the call)
    if has_bam:
        layers.append(LayerAvailability(
            "alignment_consensus", True, 
            "ICB consensus from 3 callers (already computed)",
            1.0
        ))
    else:
        layers.append(LayerAvailability("alignment_consensus", False, "No BAM available", 0.0))
    
    # Layer 2: Local assembly (LAR)
    if has_bam and local_coverage is not None and local_coverage >= 10:
        layers.append(LayerAvailability(
            "local_assembly", True,
            f"Sufficient coverage ({local_coverage:.0f}x)",
            1.0
        ))
    elif has_bam and local_coverage is not None:
        layers.append(LayerAvailability(
            "local_assembly", False,
            f"Insufficient coverage ({local_coverage:.0f}x, need ≥10x)",
            0.0
        ))
    elif has_bam:
        layers.append(LayerAvailability(
            "local_assembly", True,
            "Coverage unknown; will attempt",
            0.8
        ))
    else:
        layers.append(LayerAvailability("local_assembly", False, "No BAM available", 0.0))
    
    # Layer 3: Depth signature
    if sv_type in ('DEL', 'DUP') and sv_size >= 100 and has_bam:
        layers.append(LayerAvailability(
            "depth_signature", True,
            f"Applicable for {sv_type} ≥100 bp",
            1.0
        ))
    elif sv_type in ('DEL', 'DUP') and sv_size < 100:
        layers.append(LayerAvailability(
            "depth_signature", False,
            f"{sv_type} too small ({sv_size} bp, need ≥100 bp)",
            0.0
        ))
    elif sv_type in ('INS',):
        layers.append(LayerAvailability(
            "depth_signature", True,
            "Limited applicability for insertions",
            0.5
        ))
    elif sv_type in ('INV', 'BND'):
        layers.append(LayerAvailability(
            "depth_signature", False,
            f"Not applicable for {sv_type} (copy-neutral)",
            0.0
        ))
    else:
        layers.append(LayerAvailability("depth_signature", False, f"Unknown type: {sv_type}", 0.0))
    
    # Layer 4: k-mer spectrum
    if sv_type in ('DEL', 'INS') and has_raw_reads and has_reference:
        if in_repeat_region:
            layers.append(LayerAvailability(
                "kmer_spectrum", True,
                "Applicable but k-mers may be non-unique in repeat region",
                0.5  # Reduced ceiling due to repeat confound
            ))
        else:
            layers.append(LayerAvailability(
                "kmer_spectrum", True,
                "Fully applicable",
                1.0
            ))
    elif sv_type in ('DEL', 'INS') and (not has_raw_reads or not has_reference):
        layers.append(LayerAvailability(
            "kmer_spectrum", False,
            f"Missing: {'raw reads' if not has_raw_reads else 'reference'}",
            0.0
        ))
    elif sv_type in ('INV', 'DUP', 'BND'):
        layers.append(LayerAvailability(
            "kmer_spectrum", False,
            f"Not applicable for {sv_type}",
            0.0
        ))
    else:
        layers.append(LayerAvailability("kmer_spectrum", False, f"Unknown type: {sv_type}", 0.0))
    
    # Layer 5: Breakpoint junction
    if has_bam:
        layers.append(LayerAvailability(
            "breakpoint_junction", True,
            "Always applicable with BAM",
            1.0
        ))
    else:
        layers.append(LayerAvailability("breakpoint_junction", False, "No BAM available", 0.0))
    
    # Compute overall metrics
    available_layers = [l for l in layers if l.available]
    completeness = len(available_layers) / len(layers) if layers else 0
    
    # Max achievable T-score with default weights
    weights = {
        'alignment_consensus': 0.10,
        'local_assembly': 0.25,
        'depth_signature': 0.20,
        'kmer_spectrum': 0.25,
        'breakpoint_junction': 0.20,
    }
    
    max_score = sum(weights[l.layer_name] * l.max_achievable_score 
                    for l in layers if l.available)
    
    # Determine tier
    if completeness >= 0.8:
        tier = TriangulabilityTier.FULLY_TRIANGULABLE
    elif completeness >= 0.5:
        tier = TriangulabilityTier.PARTIALLY_TRIANGULABLE
    elif completeness >= 0.2:
        tier = TriangulabilityTier.LIMITED
    else:
        tier = TriangulabilityTier.NOT_TRIANGULABLE
    
    # Generate recommendations
    recommendations = []
    for layer in layers:
        if not layer.available and layer.limitation:
            recommendations.append(f"{layer.layer_name}: {layer.limitation}")
    
    if completeness < 0.5:
        recommendations.append(
            "WARNING: Less than 50% of evidence layers available. "
            "T-score will have wide uncertainty. Consider experimental validation."
        )
    
    return TriangulabilityReport(
        sv_id=sv_id,
        sv_type=sv_type,
        sv_size=sv_size,
        tier=tier,
        layers=layers,
        overall_completeness=round(completeness, 2),
        max_achievable_tscore=round(max_score, 4),
        recommendations=recommendations
    )


if __name__ == '__main__':
    # Test different SV scenarios
    scenarios = [
        ("DEL_large", "DEL", 5000, False, 58.0),
        ("DEL_small", "DEL", 30, False, 58.0),
        ("INV", "INV", 50000, False, 58.0),
        ("DEL_repeat", "DEL", 2000, True, 58.0),
        ("DEL_lowcov", "DEL", 5000, False, 5.0),
    ]
    
    print(f"\n{'='*70}")
    print(f"  VALID-SV: Triangulability Assessment")
    print(f"{'='*70}")
    
    for sv_id, sv_type, size, repeat, cov in scenarios:
        report = assess_triangulability(
            sv_id, sv_type, size,
            in_repeat_region=repeat,
            local_coverage=cov
        )
        print(f"\n  [{report.tier.value.upper()}] {report.sv_id} "
              f"({report.sv_type}, {report.sv_size} bp)")
        print(f"  Completeness: {report.overall_completeness:.0%}")
        print(f"  Max achievable T-score: {report.max_achievable_tscore:.3f}")
        for layer in report.layers:
            status = "✓" if layer.available else "✗"
            print(f"    {status} {layer.layer_name}: {layer.limitation}")
    
    print(f"\n{'='*70}\n")

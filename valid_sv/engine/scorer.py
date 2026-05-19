#!/usr/bin/env python3
"""
Triangulation Scorer
=====================
Combines evidence from multiple orthogonal layers into a single
triangulation score (T-score) per structural variant.

Inspired by OrthoGarden's all-vs-all orthology inference:
When no single truth set exists, combine all available independent
signals. Agreement across uncorrelated evidence layers implies
high probability of a true variant.

Author: VALID-SV / FUNGUS-SV
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from enum import Enum
import json


class TScoreTier(Enum):
    """Interpretation of triangulation score."""
    TRIPLE_TRIANGULATED = (0.80, 1.00, "Supported by ≥3 orthogonal layers; FDR <5%")
    DOUBLE_CONFIRMED = (0.60, 0.79, "Supported by ≥2 layers; FDR 5-20%")
    SINGLE_LINE = (0.40, 0.59, "Single line of evidence; FDR 20-50%")
    WEAK = (0.20, 0.39, "Weak/conflicting evidence; FDR >50%")
    CONTRADICTED = (0.00, 0.19, "Evidence contradicts; likely false positive")

    def __init__(self, low: float, high: float, description: str):
        self.low = low
        self.high = high
        self.description = description

    @classmethod
    def classify(cls, score: float) -> 'TScoreTier':
        for tier in cls:
            if tier.low <= score <= tier.high:
                return tier
        return cls.CONTRADICTED if score < 0 else cls.TRIPLE_TRIANGULATED


@dataclass
class LayerResult:
    """Result from a single evidence layer."""
    layer_name: str
    evidence_score: float       # 0.0 (contradicts) to 1.0 (strongly supports)
    verdict: str                # Text verdict from the layer
    available: bool             # Whether this layer was run
    weight: float               # Assigned weight in triangulation
    details: str = ""


@dataclass
class TriangulationResult:
    """Complete triangulation result for one SV."""
    sv_id: str
    sv_type: str
    sv_chrom: str
    sv_start: int
    sv_end: int
    sv_size: int
    icb_support: int            # Number of callers agreeing (1-3)
    
    # Layer results
    layers: List[LayerResult]
    
    # Computed scores
    t_score: float              # Weighted triangulation score (0-1)
    t_score_unweighted: float   # Simple mean (0-1)
    layers_available: int       # Number of layers that ran
    layers_agreeing: int        # Number of layers supporting
    layers_contradicting: int   # Number of layers contradicting
    
    # Interpretation
    tier: TScoreTier
    estimated_fdr: float        # From mixture model or lookup
    interpretation: str
    
    # Uncertainty
    score_uncertainty: float    # Standard deviation of layer scores
    completeness: float         # Fraction of possible layers run
    
    def to_dict(self) -> dict:
        """Serialize for JSON output."""
        return {
            'sv_id': self.sv_id,
            'sv_type': self.sv_type,
            'sv_chrom': self.sv_chrom,
            'sv_start': self.sv_start,
            'sv_end': self.sv_end,
            'sv_size': self.sv_size,
            'icb_support': self.icb_support,
            't_score': round(self.t_score, 4),
            'tier': self.tier.name,
            'estimated_fdr': round(self.estimated_fdr, 4),
            'layers_available': self.layers_available,
            'layers_agreeing': self.layers_agreeing,
            'layers_contradicting': self.layers_contradicting,
            'score_uncertainty': round(self.score_uncertainty, 4),
            'completeness': round(self.completeness, 2),
            'interpretation': self.interpretation,
            'layer_details': [
                {
                    'name': l.layer_name,
                    'score': l.evidence_score,
                    'verdict': l.verdict,
                    'available': l.available,
                    'weight': l.weight,
                    'details': l.details
                }
                for l in self.layers
            ]
        }


# Updated weights based on Liu et al. (2024) Nature Comms:
# - LAR/assembly provides breakpoint precision (near-zero shift) → highest weight
# - Depth is independent of alignment → high weight
# - k-mer is alignment-free → high weight, but lower for small SVs
# - Breakpoint junction relies on alignment → moderate weight
# These are PRIOR weights — spike-in calibration will refine them
DEFAULT_WEIGHTS = {
    'alignment_consensus': 0.0,     # Excluded: circular validation
    'local_assembly': 0.30,         # Highest: assembly confirms exact breakpoints (Liu Fig 3)
    'depth_signature': 0.25,        # Independent of alignment (alignment-free)
    'kmer_spectrum': 0.25,          # Independent of alignment (alignment-free)
    'breakpoint_junction': 0.20,    # Relies on alignment SA tags (moderate)
}

class TriangulationScorer:
    """
    Combines evidence from multiple layers into a T-score.
    
    The T-score represents the probability-weighted consensus
    across independent evidence sources. It is NOT a probability
    itself, but correlates with likelihood of a true variant.
    
    Calibration against synthetic benchmarks is required to
    map T-scores to estimated FDR.
    """
    
    def __init__(self, weights: Optional[Dict[str, float]] = None,
                 calibration: Optional[Dict[str, float]] = None,
                 include_consensus: bool = False):
        """
        Args:
            weights: Layer weights (default: DEFAULT_WEIGHTS)
            calibration: T-score to FDR mapping from benchmarks
            include_consensus: If True, include alignment_consensus in T-score.
                Default False to avoid circular validation (ICB output as input).
        """
        self.weights = weights or DEFAULT_WEIGHTS
        self.calibration = calibration or self._default_calibration()
        self.include_consensus = include_consensus
    
    @staticmethod
    def _default_calibration() -> Dict[str, float]:
        """Conservative default calibration. Replace with benchmark results."""
        return {
            '0.80': 0.05,   # T≥0.80 → est. FDR <5%
            '0.60': 0.15,   # T≥0.60 → est. FDR <15%
            '0.40': 0.35,   # T≥0.40 → est. FDR <35%
            '0.20': 0.60,   # T≥0.20 → est. FDR <60%
        }
    
    def score(self, sv_id: str, sv_type: str, sv_chrom: str,
              sv_start: int, sv_end: int, icb_support: int,
              layer_results: List[LayerResult]) -> TriangulationResult:
        """
        Compute triangulation score from evidence layers.
        
        Args:
            sv_id: SV identifier
            sv_type: DEL, INS, DUP, INV, BND
            sv_chrom: Chromosome/contig
            sv_start: Start position
            sv_end: End position
            icb_support: Number of callers (1-3)
            layer_results: Results from each evidence layer
        
        Returns:
            TriangulationResult with T-score and interpretation
        """
        # Filter out alignment_consensus unless explicitly included (circular validation)
        available_layers = [l for l in layer_results 
                          if l.available and (self.include_consensus or l.layer_name != 'alignment_consensus')]
        
        if not available_layers:
            # No evidence at all
            return TriangulationResult(
                sv_id=sv_id, sv_type=sv_type, sv_chrom=sv_chrom,
                sv_start=sv_start, sv_end=sv_end,
                sv_size=sv_end - sv_start, icb_support=icb_support,
                layers=layer_results,
                t_score=0.0, t_score_unweighted=0.0,
                layers_available=0, layers_agreeing=0, layers_contradicting=0,
                tier=TScoreTier.CONTRADICTED,
                estimated_fdr=1.0,
                interpretation="No evidence available for validation",
                score_uncertainty=0.0, completeness=0.0
            )
        
        # Compute weighted T-score
        weighted_sum = 0.0
        weight_total = 0.0
        unweighted_sum = 0.0
        
        agreeing = 0
        contradicting = 0
        
        for layer in available_layers:
            w = self.weights.get(layer.layer_name, 0.20)
            weighted_sum += layer.evidence_score * w
            weight_total += w
            unweighted_sum += layer.evidence_score
            
            if layer.evidence_score >= 0.6:
                agreeing += 1
            elif layer.evidence_score <= 0.2:
                contradicting += 1
        
        t_score = weighted_sum / weight_total if weight_total > 0 else 0.0
        t_score_unweighted = unweighted_sum / len(available_layers)
        
        # Score uncertainty (std dev of layer scores)
        scores = [l.evidence_score for l in available_layers]
        score_uncertainty = float(np.std(scores)) if len(scores) > 1 else 0.0
        
        # Completeness
        total_possible = len(layer_results)
        completeness = len(available_layers) / total_possible if total_possible > 0 else 0
        
        # Apply completeness penalty
        if completeness < 0.5:
            t_score *= completeness  # Penalize low completeness
        
        # Determine tier
        tier = TScoreTier.classify(t_score)
        
        # Estimate FDR from calibration
        estimated_fdr = self._estimate_fdr(t_score)
        
        # Generate interpretation
        interpretation = self._interpret(
            t_score, tier, available_layers, agreeing, contradicting
        )
        
        return TriangulationResult(
            sv_id=sv_id, sv_type=sv_type, sv_chrom=sv_chrom,
            sv_start=sv_start, sv_end=sv_end,
            sv_size=sv_end - sv_start, icb_support=icb_support,
            layers=layer_results,
            t_score=round(t_score, 4),
            t_score_unweighted=round(t_score_unweighted, 4),
            layers_available=len(available_layers),
            layers_agreeing=agreeing,
            layers_contradicting=contradicting,
            tier=tier,
            estimated_fdr=round(estimated_fdr, 4),
            interpretation=interpretation,
            score_uncertainty=round(score_uncertainty, 4),
            completeness=round(completeness, 2)
        )
    
    def _estimate_fdr(self, t_score: float) -> float:
        """Map T-score to estimated FDR using calibration."""
        thresholds = sorted([(float(k), v) for k, v in self.calibration.items()],
                           reverse=True)
        for threshold, fdr in thresholds:
            if t_score >= threshold:
                return fdr
        return 1.0
    
    def _interpret(self, t_score: float, tier: TScoreTier,
                   layers: List[LayerResult], agreeing: int,
                   contradicting: int) -> str:
        """Generate human-readable interpretation."""
        parts = []
        parts.append(f"T-score: {t_score:.3f} ({tier.name})")
        parts.append(f"Evidence: {agreeing} layers support, {contradicting} contradict")
        
        if tier == TScoreTier.TRIPLE_TRIANGULATED:
            parts.append("Suitable for functional follow-up (experimental validation still recommended)")
        elif tier == TScoreTier.DOUBLE_CONFIRMED:
            parts.append("Prioritized candidate; recommend targeted PCR validation")
        elif tier == TScoreTier.SINGLE_LINE:
            parts.append("Exploratory only; do not base conclusions solely on this SV")
        elif tier == TScoreTier.WEAK:
            parts.append("Low confidence; very likely artifact")
        else:
            parts.append("Probable false positive or reference assembly error")
        
        return " | ".join(parts)


# Convenience function
def compute_tscore(sv_id: str, sv_type: str, sv_chrom: str,
                   sv_start: int, sv_end: int, icb_support: int,
                   layer_results: List[LayerResult],
                   weights: Optional[Dict[str, float]] = None) -> TriangulationResult:
    """Compute T-score for a single SV."""
    scorer = TriangulationScorer(weights=weights)
    return scorer.score(sv_id, sv_type, sv_chrom, sv_start, sv_end,
                        icb_support, layer_results)


# Need numpy for std calculation
import numpy as np


if __name__ == '__main__':
    # Test with simulated layer results
    test_layers = [
        LayerResult("alignment_consensus", 1.0, "3/3 callers", True, 0.10),
        LayerResult("local_assembly", 0.95, "confirmed", True, 0.25),
        LayerResult("depth_signature", 0.80, "consistent", True, 0.20),
        LayerResult("kmer_spectrum", 0.90, "strong_support", True, 0.25),
        LayerResult("breakpoint_junction", 0.70, "confirmed", True, 0.20),
    ]
    
    scorer = TriangulationScorer()
    result = scorer.score(
        "test_DEL_001", "DEL", "chr1", 50000, 60000, 3, test_layers
    )
    
    print(f"\n{'='*60}")
    print(f"  VALID-SV: Triangulation Scorer — Test")
    print(f"{'='*60}")
    print(f"  SV: {result.sv_id}")
    print(f"  T-score: {result.t_score:.4f}")
    print(f"  Tier: {result.tier.name}")
    print(f"  Est. FDR: {result.estimated_fdr:.1%}")
    print(f"  Layers: {result.layers_agreeing} agree, {result.layers_contradicting} contradict")
    print(f"  Interpretation: {result.interpretation}")
    print(f"{'='*60}\n")

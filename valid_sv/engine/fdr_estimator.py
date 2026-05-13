#!/usr/bin/env python3
"""
FDR Estimator for Triangulation Scores
=======================================
Uses empirical distribution of T-scores to estimate
false discovery rates without a truth set.

Method: The lowest-scoring null_proportion of SVs are
assumed to be predominantly false positives. Their score
distribution is used to estimate FDR at each threshold.

CRITICAL CAVEAT: This is an ESTIMATE, not a measurement.
Calibrate with synthetic benchmarks for publication.

Author: VALID-SV / FUNGUS-SV
"""

import numpy as np
from typing import List, Dict, Optional
from dataclasses import dataclass


@dataclass
class FDREstimate:
    """Estimated FDR at various T-score thresholds."""
    total_svs: int
    method: str
    null_proportion: float
    thresholds: Dict[float, float]  # T-score threshold → estimated FDR
    caveat: str = "APPROXIMATE — calibrate with synthetic benchmarks"


def estimate_fdr(scores: List[float],
                 thresholds: Optional[List[float]] = None,
                 null_proportion: float = 0.30) -> FDREstimate:
    """
    Estimate FDR using empirical null distribution.
    
    Args:
        scores: List of T-scores for all SVs
        thresholds: T-score thresholds to evaluate
        null_proportion: Assumed max fraction of false positives
    
    Returns:
        FDREstimate with per-threshold FDR values
    
    Method:
        1. Sort scores ascending
        2. Bottom null_proportion = empirical null distribution
        3. For each threshold, FDR = (expected null above) / (total above)
    """
    if thresholds is None:
        thresholds = [0.2, 0.4, 0.6, 0.8]
    
    scores_array = np.array(sorted(scores))
    n = len(scores_array)
    
    if n < 10:
        return FDREstimate(
            total_svs=n, method="insufficient_data",
            null_proportion=null_proportion,
            thresholds={t: 1.0 for t in thresholds},
            caveat="Too few SVs for FDR estimation (need ≥10)"
        )
    
    # Empirical null: bottom null_proportion of scores
    n_null = max(5, int(n * null_proportion))
    null_scores = scores_array[:n_null]
    
    fdr_at_threshold = {}
    for t in thresholds:
        # Fraction of null SVs above threshold
        null_rate_above = np.mean(null_scores >= t) if len(null_scores) > 0 else 1.0
        
        # Expected false positives above threshold
        expected_fp = n * null_proportion * null_rate_above
        
        # Total SVs above threshold
        n_total_above = np.sum(scores_array >= t)
        
        # FDR
        fdr = expected_fp / n_total_above if n_total_above > 0 else 1.0
        fdr_at_threshold[t] = round(min(1.0, fdr), 4)
    
    return FDREstimate(
        total_svs=n,
        method="empirical_null",
        null_proportion=null_proportion,
        thresholds=fdr_at_threshold,
    )


def estimate_fdr_mixture(scores: List[float],
                          thresholds: Optional[List[float]] = None) -> FDREstimate:
    """
    Mixture model FDR estimate (requires scikit-learn).
    Falls back to empirical method if sklearn unavailable.
    """
    if thresholds is None:
        thresholds = [0.2, 0.4, 0.6, 0.8]
    
    try:
        from sklearn.mixture import GaussianMixture
        from scipy.stats import norm
        
        scores_array = np.array(scores).reshape(-1, 1)
        
        gmm = GaussianMixture(n_components=2, random_state=42, n_init=5)
        gmm.fit(scores_array)
        
        means = gmm.means_.flatten()
        vars_ = gmm.covariances_.flatten()
        weights = gmm.weights_
        
        true_idx = np.argmax(means)
        false_idx = 1 - true_idx
        
        fdr_at_threshold = {}
        for t in thresholds:
            false_density = weights[false_idx] * (1.0 - norm.cdf(t, means[false_idx], np.sqrt(vars_[false_idx])))
            true_density = weights[true_idx] * (1.0 - norm.cdf(t, means[true_idx], np.sqrt(vars_[true_idx])))
            total_density = false_density + true_density
            
            fdr = false_density / total_density if total_density > 0 else 1.0
            fdr_at_threshold[t] = round(min(1.0, fdr), 4)
        
        return FDREstimate(
            total_svs=len(scores),
            method="gaussian_mixture",
            null_proportion=round(weights[false_idx], 3),
            thresholds=fdr_at_threshold,
            caveat="Mixture model — validate with benchmarks"
        )
    
    except ImportError:
        return estimate_fdr(scores, thresholds)
    except Exception:
        return estimate_fdr(scores, thresholds)


if __name__ == '__main__':
    np.random.seed(42)
    
    # Simulate score distribution
    true_scores = np.random.beta(5, 1.5, 80) * 0.5 + 0.5
    false_scores = np.random.beta(1.5, 5, 20) * 0.6
    all_scores = np.concatenate([true_scores, false_scores])
    np.random.shuffle(all_scores)
    
    result = estimate_fdr(all_scores.tolist())
    print(f"\n{'='*60}")
    print(f"  VALID-SV: FDR Estimator — Test")
    print(f"{'='*60}")
    print(f"  Method: {result.method}")
    print(f"  Null proportion: {result.null_proportion:.1%}")
    print(f"  Total SVs: {result.total_svs}")
    for t, fdr in sorted(result.thresholds.items()):
        print(f"    T ≥ {t:.1f}: est. FDR = {fdr:.1%}")
    print(f"  {result.caveat}")
    print(f"{'='*60}\n")

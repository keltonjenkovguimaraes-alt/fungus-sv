"""
VALID-SV: Validation by Independent Layer Integration & Detection
==================================================================
A computational-only structural variant validation package for
non-model haploid organisms with PacBio HiFi data.

When no benchmark truth set exists, VALID-SV triangulates
multiple orthogonal lines of evidence with uncorrelated failure modes
to estimate SV confidence and false discovery rate.

Architecture inspired by:
  - NEEDLE (Ko & Brandizzi 2025): Separate prediction from validation
  - OrthoGarden (Turner et al. 2026): All-vs-all evidence inference replaces reference truth
  - DMRichR FAIRification (Salam et al. 2025): Make model-organism tools work for any annotated genome

Author: FUNGUS-SV team
Version: 0.1.0
Status: Development - NOT FOR PRODUCTION USE WITHOUT EXPERIMENTAL VALIDATION
"""

__version__ = "0.1.0"
__status__ = "Development"
__author__ = "FUNGUS-SV team"

# Re-export key functions for clean imports
from .evidence.layer_depth import DepthEvidence, analyze_depth_signature
from .evidence.layer_kmer import KmerEvidence, analyze_kmer_spectrum
from .evidence.layer_breakpoint import BreakpointEvidence, analyze_breakpoint_junctions
from .quality.triangulability import assess_triangulability
from .engine.scorer import TriangulationScorer, compute_tscore
from .engine.fdr_estimator import estimate_fdr

__all__ = [
    'DepthEvidence', 'analyze_depth_signature',
    'KmerEvidence', 'analyze_kmer_spectrum',
    'BreakpointEvidence', 'analyze_breakpoint_junctions',
    'assess_triangulability',
    'TriangulationScorer', 'compute_tscore',
    'estimate_fdr',
]

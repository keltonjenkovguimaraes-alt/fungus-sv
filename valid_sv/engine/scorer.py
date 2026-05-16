from typing import List, Optional
from dataclasses import dataclass
from enum import Enum

class TScoreTier(Enum):
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"

@dataclass
class LayerResult:
    layer_name: str
    score: float
    evidence: str
    available: bool
    weight: float
    details: str

@dataclass
class TriangulationResult:
    sv_id: str
    svtype: str
    chrom: str
    pos: int
    end: int
    support: int
    t_score: float
    confidence_tier: str
    layer_results: List[LayerResult]
    
    def to_dict(self):
        return {
            'sv_id': self.sv_id,
            'svtype': self.svtype,
            'chrom': self.chrom,
            'pos': self.pos,
            'end': self.end,
            'support': self.support,
            't_score': self.t_score,
            'confidence_tier': self.confidence_tier,
            'layer_results': [
                {
                    'layer_name': lr.layer_name,
                    'score': lr.score,
                    'evidence': lr.evidence,
                    'available': lr.available,
                    'weight': lr.weight,
                    'details': lr.details
                }
                for lr in self.layer_results
            ]
        }

class TriangulationScorer:
    def score(self, sv_id, svtype, chrom, pos, end, support, layer_results):
        # Simple scoring for test
        available_layers = [lr for lr in layer_results if lr.available]
        if available_layers:
            weighted_score = sum(lr.score * lr.weight for lr in available_layers)
            total_weight = sum(lr.weight for lr in available_layers)
            t_score = weighted_score / total_weight if total_weight > 0 else 0
        else:
            t_score = 0
        
        if t_score > 0.7:
            tier = "HIGH"
        elif t_score > 0.4:
            tier = "MEDIUM"
        else:
            tier = "LOW"
        
        return TriangulationResult(sv_id, svtype, chrom, pos, end, support, t_score, tier, layer_results)

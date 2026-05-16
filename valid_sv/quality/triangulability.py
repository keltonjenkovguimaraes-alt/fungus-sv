from enum import Enum
from types import SimpleNamespace

class TriangulabilityTier(Enum):
    FULLY_TRIANGULABLE = "FULLY_TRIANGULABLE"
    PARTIALLY_TRIANGULABLE = "PARTIALLY_TRIANGULABLE"
    LIMITED = "LIMITED"
    NOT_TRIANGULABLE = "NOT_TRIANGULABLE"

class TriangulabilityReport:
    def __init__(self, sv_id, tier, layers):
        self.sv_id = sv_id
        self.tier = tier
        self.layers = layers

def assess_triangulability(sv_id, svtype, size, has_reference, has_raw_reads, has_bam):
    layers = []
    layers.append(SimpleNamespace(layer_name="alignment_consensus", available=True))
    layers.append(SimpleNamespace(layer_name="local_assembly", available=has_reference and has_raw_reads))
    layers.append(SimpleNamespace(layer_name="depth_signature", available=has_bam and svtype in ['DEL', 'DUP']))
    layers.append(SimpleNamespace(layer_name="kmer_spectrum", available=has_raw_reads and svtype in ['DEL', 'INS']))
    layers.append(SimpleNamespace(layer_name="breakpoint_junction", available=has_bam))
    
    available_count = sum(1 for l in layers if l.available)
    if available_count >= 4:
        tier = TriangulabilityTier.FULLY_TRIANGULABLE
    elif available_count >= 2:
        tier = TriangulabilityTier.PARTIALLY_TRIANGULABLE
    else:
        tier = TriangulabilityTier.LIMITED
    
    return TriangulabilityReport(sv_id, tier, layers)

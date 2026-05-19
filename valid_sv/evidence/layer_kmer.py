#!/usr/bin/env python3
"""
Evidence Layer 4: k-mer Spectrum Analysis (v3)
================================================
Pre-built jellyfish database approach. The database is built ONCE
from all reads, then queried per SV. This avoids the per-SV
jellyfish overhead that was timing out on full datasets.

Requires: jellyfish >= 2.0 (conda install -c bioconda jellyfish)

Author: VALID-SV / FUNGUS-SV
"""

import subprocess
import tempfile
import os
import shutil
import numpy as np
from dataclasses import dataclass
from typing import List, Optional
from enum import Enum


class KmerVerdict(Enum):
    STRONG_SUPPORT = "strong_support"
    WEAK_SUPPORT = "weak_support"
    AMBIGUOUS = "ambiguous"
    CONTRADICTS = "contradicts"
    NOT_APPLICABLE = "not_applicable"
    INSUFFICIENT_DATA = "insufficient_data"


@dataclass
class KmerEvidence:
    sv_id: str
    sv_type: str
    verdict: KmerVerdict
    kmer_ratio: float
    ref_kmer_mean: float
    flank_kmer_mean: float
    evidence_score: float
    details: str = ""
    kmers_checked: int = 0


# Module-level cache for the jellyfish DB path
_JELLYFISH_DB = None


def build_kmer_database(fastq_path: str, output_dir: str = "/tmp/valid_sv_kmers",
                        k: int = 31, threads: int = 4) -> str:
    """
    Build jellyfish k-mer database from raw reads.
    Call ONCE per dataset, not per SV.
    
    Returns path to the .jf file.
    """
    global _JELLYFISH_DB
    
    os.makedirs(output_dir, exist_ok=True)
    jf_file = os.path.join(output_dir, "reads.jf")
    
    # Skip if already built
    if os.path.exists(jf_file) and os.path.getsize(jf_file) > 1000:
        _JELLYFISH_DB = jf_file
        return jf_file
    
    print(f"[k-mer] Building jellyfish database (one-time, ~5-10 min for full dataset)...")
    print(f"[k-mer] This will be cached at {jf_file}")
    
    if fastq_path.endswith('.gz'):
        cmd = (f"zcat {fastq_path} | jellyfish count -m {k} -s 100M -t {threads} "
               f"-o {jf_file} /dev/stdin")
        shell = True
    else:
        cmd = ['jellyfish', 'count', '-m', str(k), '-s', '100M', '-t', str(threads),
               '-o', jf_file, fastq_path]
        shell = False
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, shell=shell)
    
    if result.returncode != 0:
        raise RuntimeError(f"Jellyfish database build failed: {result.stderr}")
    
    _JELLYFISH_DB = jf_file
    print(f"[k-mer] Database built: {os.path.getsize(jf_file)/1e6:.1f} MB")
    return jf_file


def set_database_path(jf_path: str):
    """Set path to pre-built jellyfish database."""
    global _JELLYFISH_DB
    if os.path.exists(jf_path):
        _JELLYFISH_DB = jf_path
    else:
        raise FileNotFoundError(f"Jellyfish DB not found: {jf_path}")


def extract_kmers_from_reference(reference_path: str, chrom: str,
                                  start: int, end: int, k: int = 31,
                                  stride: int = 10) -> List[str]:
    """Extract k-mers from a reference genome region."""
    region = f"{chrom}:{start}-{end}"
    
    result = subprocess.run(
        ['samtools', 'faidx', reference_path, region],
        capture_output=True, text=True, timeout=30
    )
    
    if result.returncode != 0:
        return []
    
    lines = result.stdout.strip().split('\n')
    seq = ''.join(lines[1:]).upper()
    
    kmers = []
    for i in range(0, len(seq) - k + 1, stride):
        kmer = seq[i:i + k]
        if 'N' not in kmer:
            kmers.append(kmer)
    
    return kmers


def query_kmer_counts(jf_file: str, kmers: List[str]) -> List[int]:
    """Query pre-built jellyfish DB for k-mer counts."""
    if not kmers:
        return []
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        for i, kmer in enumerate(kmers):
            f.write(f">k{i}\n{kmer}\n")
        query_file = f.name
    
    try:
        result = subprocess.run(
            ['jellyfish', 'query', jf_file, '-s', query_file],
            capture_output=True, text=True, timeout=30
        )
    finally:
        os.unlink(query_file)
    
    if result.returncode != 0:
        return [-1] * len(kmers)
    
    counts = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                counts.append(int(parts[1]))
            except ValueError:
                counts.append(-1)
        else:
            counts.append(0)
    
    return counts


def analyze_kmer_spectrum(fastq_path: str, reference_path: str,
                           sv_id: str, sv_type: str, chrom: str,
                           start: int, end: int, k: int = 31,
                           flank_size: int = 1000,
                           jf_db: Optional[str] = None) -> KmerEvidence:
    """
    Analyze k-mer presence/absence for SV validation.
    Uses pre-built or cached jellyfish database.
    """
    global _JELLYFISH_DB
    
    if sv_type not in ('DEL', 'INS'):
        return KmerEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=KmerVerdict.NOT_APPLICABLE,
            kmer_ratio=1.0, ref_kmer_mean=0, flank_kmer_mean=0,
            evidence_score=0.0,
            details=f"k-mer analysis not applicable for {sv_type}"
        )
    
    # Determine DB path
    db_path = jf_db or _JELLYFISH_DB
    if not db_path or not os.path.exists(db_path):
        return KmerEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=KmerVerdict.INSUFFICIENT_DATA,
            kmer_ratio=0, ref_kmer_mean=0, flank_kmer_mean=0,
            evidence_score=0.0,
            details="No k-mer database. Run build_kmer_database() first."
        )
    
    if not shutil.which('jellyfish'):
        return KmerEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=KmerVerdict.INSUFFICIENT_DATA,
            kmer_ratio=0, ref_kmer_mean=0, flank_kmer_mean=0,
            evidence_score=0.0,
            details="jellyfish not installed"
        )
    
    # Extract kmers
    ref_kmers = extract_kmers_from_reference(reference_path, chrom, start, end, k)
    
    flank_start = max(1, start - flank_size)
    flank_end = start
    flank_kmers = extract_kmers_from_reference(reference_path, chrom, flank_start, flank_end, k)
    
    if not ref_kmers or not flank_kmers:
        return KmerEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=KmerVerdict.INSUFFICIENT_DATA,
            kmer_ratio=0, ref_kmer_mean=0, flank_kmer_mean=0,
            evidence_score=0.0,
            details="Could not extract k-mers from reference"
        )
    
    # Query counts (sample up to 100 kmers for speed)
    ref_sample = ref_kmers[:min(100, len(ref_kmers))]
    flank_sample = flank_kmers[:min(100, len(flank_kmers))]
    
    ref_counts = query_kmer_counts(db_path, ref_sample)
    flank_counts = query_kmer_counts(db_path, flank_sample)
    
    ref_counts = [c for c in ref_counts if c >= 0]
    flank_counts = [c for c in flank_counts if c >= 0]
    
    if not ref_counts or not flank_counts:
        return KmerEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=KmerVerdict.INSUFFICIENT_DATA,
            kmer_ratio=0, ref_kmer_mean=0, flank_kmer_mean=0,
            evidence_score=0.0,
            details="No valid k-mer counts obtained"
        )
    
    ref_median = np.median(ref_counts)
    flank_median = np.median(flank_counts)
    
    if flank_median == 0:
        return KmerEvidence(
            sv_id=sv_id, sv_type=sv_type,
            verdict=KmerVerdict.AMBIGUOUS,
            kmer_ratio=0, ref_kmer_mean=ref_median, flank_kmer_mean=flank_median,
            evidence_score=0.0,
            details="Control k-mers absent from reads"
        )
    
    kmer_ratio = ref_median / flank_median
    
    if sv_type == 'DEL':
        if kmer_ratio < 0.10:
            verdict = KmerVerdict.STRONG_SUPPORT
            score = 1.0
            details = f"Deleted kmers depleted {1-kmer_ratio:.0%} (ratio={kmer_ratio:.4f})"
        elif kmer_ratio < 0.30:
            verdict = KmerVerdict.WEAK_SUPPORT
            score = 0.7
            details = f"Deleted kmers partially depleted (ratio={kmer_ratio:.4f})"
        elif kmer_ratio > 0.70:
            verdict = KmerVerdict.CONTRADICTS
            score = 0.0
            details = f"Deleted kmers present at normal level (ratio={kmer_ratio:.4f}) — possible FP"
        else:
            verdict = KmerVerdict.AMBIGUOUS
            score = 0.3
            details = f"Intermediate depletion (ratio={kmer_ratio:.4f})"
    
    elif sv_type == 'INS':
        if kmer_ratio > 0.50:
            verdict = KmerVerdict.STRONG_SUPPORT
            score = min(1.0, kmer_ratio)
            details = f"Insert kmers detected (ratio={kmer_ratio:.4f})"
        elif kmer_ratio > 0.20:
            verdict = KmerVerdict.WEAK_SUPPORT
            score = 0.5
            details = f"Insert kmers partially detected (ratio={kmer_ratio:.4f})"
        elif kmer_ratio < 0.05:
            verdict = KmerVerdict.CONTRADICTS
            score = 0.0
            details = "Insert kmers absent — contradicts insertion"
        else:
            verdict = KmerVerdict.AMBIGUOUS
            score = 0.3
            details = f"Unclear (ratio={kmer_ratio:.4f})"
    
    else:
        verdict = KmerVerdict.NOT_APPLICABLE
        score = 0.0
        details = f"Unknown type: {sv_type}"
    
    return KmerEvidence(
        sv_id=sv_id, sv_type=sv_type,
        verdict=verdict,
        kmer_ratio=round(kmer_ratio, 4),
        ref_kmer_mean=round(ref_median, 2),
        flank_kmer_mean=round(flank_median, 2),
        evidence_score=round(score, 4),
        details=details,
        kmers_checked=len(ref_counts) + len(flank_counts)
    )

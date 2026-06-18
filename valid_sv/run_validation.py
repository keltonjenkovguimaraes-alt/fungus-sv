#!/usr/bin/env python3
"""
VALID-SV: Main Validation Pipeline
====================================
Entry point for the triangulation-based SV validation.

Architecture:
    Layer 0: ICB consensus (REPORTED, not scored)
    Layer 1: Local Assembly Refinement (LAR) - standalone
    Layer 2: Read-Depth Signature
    Layer 3: k-mer Spectrum Analysis
    Layer 4: Breakpoint Junction Analysis
    Layer 5: Ploidy Confirmation (fast bcftools pileup)
    Layer 6: Genomic Context (hard filter)
    -> Triangulation Engine -> T-score + Confidence Estimate
"""

import sys, os, argparse, json, yaml, glob, re
from pathlib import Path
from typing import List, Dict, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))

from valid_sv.evidence.layer_depth import analyze_depth_signature, DepthEvidence
from valid_sv.evidence.layer_kmer import analyze_kmer_spectrum, KmerEvidence
from valid_sv.evidence.layer_breakpoint import analyze_breakpoint_junctions, BreakpointEvidence
from valid_sv.evidence.layer_ploidy import analyze_ploidy, PloidyEvidence
from valid_sv.evidence.layer_genomic_context import analyze_genomic_context, FilterVerdict
from valid_sv.evidence.layer_lar import run_lar, LAREvidence, LARVerdict
from valid_sv.quality.triangulability import assess_triangulability, TriangulabilityReport
from valid_sv.engine.scorer import (
    TriangulationScorer, LayerResult, TriangulationResult, TScoreTier
)
from valid_sv.engine.fdr_estimator import estimate_fdr
from valid_sv.reporting.report_card import generate_report_card, generate_summary_table


def parse_consensus_vcf(vcf_path):
    svs = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 8: continue
            chrom, pos, sv_id, info = parts[0], int(parts[1]), parts[2], parts[7]
            svtype_match = re.search(r'SVTYPE=(\w+)', info)
            svtype = svtype_match.group(1) if svtype_match else 'UNK'
            end_match = re.search(r'END=(\d+)', info)
            svlen_match = re.search(r'SVLEN=(\d+)', info)
            if end_match: end = int(end_match.group(1))
            elif svlen_match: end = pos + abs(int(svlen_match.group(1)))
            else: end = pos
            support_match = re.search(r'SUPPORT=(\d+)', info)
            if not support_match:
                support_match = re.search(r'NUMCALLERS=(\d+)', info)
            support = int(support_match.group(1)) if support_match else 1
            svs.append({'id': sv_id, 'chrom': chrom, 'pos': pos, 'end': end,
                        'svtype': svtype, 'size': end - pos, 'support': support})
    return svs


def run_validation_pipeline(consensus_vcf, bam_path, reference_path, fastq_path,
                            output_dir, min_support=1, max_svs=None, skip_kmer=False,
                            jellyfish_db=None, ablation=False, threads=4):
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("  VALID-SV: Triangulation-Based SV Validation")
    print("=" * 70)
    print(f"  BAM: {bam_path}")
    print(f"  Reference: {reference_path}")
    print("=" * 70)
    print()

    # Annotation TSV
    _annotation_tsv = None
    _bam_basename = os.path.basename(bam_path)
    for _tsv_path in glob.glob('data/yeast/*_sv_annotations.tsv'):
        _tsv_strain = os.path.basename(_tsv_path).replace('_sv_annotations.tsv', '')
        if _tsv_strain in _bam_basename or _tsv_strain in reference_path:
            _annotation_tsv = _tsv_path; break
    if not _annotation_tsv:
        _annotation_tsv = 'data/yeast/S288C_sv_annotations.tsv'

    # Parse SVs
    all_svs = parse_consensus_vcf(consensus_vcf)
    print(f"  Loaded {len(all_svs)} SVs")
    svs_to_validate = [sv for sv in all_svs if sv['support'] >= min_support]
    print(f"  Validating {len(svs_to_validate)} SVs with SUPPORT >= {min_support}")
    if max_svs:
        svs_to_validate = svs_to_validate[:max_svs]

    # k-mer database
    kmer_db = None
    if jellyfish_db and os.path.exists(jellyfish_db):
        from valid_sv.evidence.layer_kmer import set_database_path
        kmer_db = jellyfish_db; set_database_path(kmer_db)
    elif fastq_path and not skip_kmer:
        from valid_sv.evidence.layer_kmer import build_kmer_database, set_database_path
        try:
            kmer_db = build_kmer_database(fastq_path, output_dir=os.path.join(output_dir, "kmer_db"))
            set_database_path(kmer_db)
        except Exception as e:
            print(f"  WARNING: k-mer DB build failed: {e}"); skip_kmer = True

    # Load calibrated weights
    _config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'config.yaml')
    with open(_config_path) as f:
        _config = yaml.safe_load(f)
    CALIBRATED_WEIGHTS = _config.get('weights', {
        'local_assembly': 0.20, 'depth_signature': 0.35,
        'kmer_spectrum': 0.15, 'breakpoint_junction': 0.30,
        'ploidy_confirmation': 0.0,
    })

    # ====== RUN PLOIDY CHECK ONCE BEFORE THE LOOP ======
    print("\n  Checking ploidy...")
    try:
        ploidy_result = analyze_ploidy(bam_path, reference_path)
        print(f"    Het rate: {ploidy_result.het_rate:.3f}, Haploid: {ploidy_result.is_haploid}")
    except Exception as e:
        ploidy_result = PloidyEvidence(0, 0, 0, 0.0, True, 1.0, f"Ploidy check failed: {e}")
        print(f"    WARNING: {e}")

    # ====== RUN EVIDENCE LAYERS ======
    print("\n  Running evidence layers...")
    scorer = TriangulationScorer(weights=CALIBRATED_WEIGHTS)
    results = []

    for i, sv in enumerate(svs_to_validate):
        if i % 50 == 0:
            print(f"    Processing {i}/{len(svs_to_validate)}...")

        layer_results = []

        # Layer 0: ICB consensus
        icb_score = {1: 0.33, 2: 0.67, 3: 1.0}.get(sv['support'], 0.33)
        layer_results.append(LayerResult("alignment_consensus", icb_score,
            f"{sv['support']}/3 callers", True, 0.0, f"ICB support: {sv['support']} callers"))

        # Layer 1: LAR
        layer_results.append(LayerResult("local_assembly", 0.0, "not_run", False,
            CALIBRATED_WEIGHTS['local_assembly'], "LAR must be run separately"))

        # Layer 2: Depth
        if sv['svtype'] in ('DEL', 'DUP'):
            try:
                depth_result = analyze_depth_signature(bam_path, sv['id'], sv['svtype'],
                    sv['chrom'], int(sv['pos']), int(sv['end']))
                layer_results.append(LayerResult("depth_signature", depth_result.evidence_score,
                    depth_result.verdict.value, True, CALIBRATED_WEIGHTS['depth_signature'],
                    depth_result.details))
            except Exception as e:
                layer_results.append(LayerResult("depth_signature", 0.0, "error", False,
                    CALIBRATED_WEIGHTS['depth_signature'], f"Depth failed: {e}"))
        else:
            layer_results.append(LayerResult("depth_signature", 0.0, "not_applicable", False,
                CALIBRATED_WEIGHTS['depth_signature'], f"Not applicable for {sv['svtype']}"))

        # Layer 3: k-mer
        if not skip_kmer and fastq_path and sv['svtype'] in ('DEL', 'INS'):
            try:
                kmer_result = analyze_kmer_spectrum(fastq_path, reference_path,
                    sv['id'], sv['svtype'], sv['chrom'], int(sv['pos']), int(sv['end']), jf_db=kmer_db)
                layer_results.append(LayerResult("kmer_spectrum", kmer_result.evidence_score,
                    kmer_result.verdict.value, True, CALIBRATED_WEIGHTS['kmer_spectrum'],
                    kmer_result.details))
            except Exception as e:
                layer_results.append(LayerResult("kmer_spectrum", 0.0, "error", False,
                    CALIBRATED_WEIGHTS['kmer_spectrum'], f"k-mer failed: {e}"))
        else:
            layer_results.append(LayerResult("kmer_spectrum", 0.0, "unavailable", False,
                CALIBRATED_WEIGHTS['kmer_spectrum'], "k-mer layer not available"))

        # Layer 4: Breakpoint
        try:
            bp_result = analyze_breakpoint_junctions(bam_path, sv['id'], sv['svtype'],
                sv['chrom'], int(sv['pos']), int(sv['end']))
            layer_results.append(LayerResult("breakpoint_junction", bp_result.evidence_score,
                bp_result.verdict.value, True, CALIBRATED_WEIGHTS['breakpoint_junction'],
                bp_result.details))
        except Exception as e:
            layer_results.append(LayerResult("breakpoint_junction", 0.0, "error", False,
                CALIBRATED_WEIGHTS['breakpoint_junction'], f"Breakpoint failed: {e}"))

        # Layer 5: Ploidy (cached from pre-check)
        layer_results.append(LayerResult("ploidy_confirmation", ploidy_result.evidence_score,
            f"het_rate={ploidy_result.het_rate:.3f}", True,
            CALIBRATED_WEIGHTS['ploidy_confirmation'], ploidy_result.details))

        # Layer 6: Genomic context
        if os.path.exists(_annotation_tsv):
            try:
                genomic_result = analyze_genomic_context(sv['id'], sv['svtype'], _annotation_tsv)
                if genomic_result.verdict == FilterVerdict.PASS: gc_score = 1.0
                elif genomic_result.verdict == FilterVerdict.FLAG: gc_score = 0.5
                else: gc_score = 0.0
                layer_results.append(LayerResult("genomic_context", gc_score,
                    genomic_result.verdict.value, True, 0.0, genomic_result.details))
            except Exception as e:
                layer_results.append(LayerResult("genomic_context", 0.0, "error", False, 0.0,
                    f"Genomic context failed: {e}"))

        # Score
        result = scorer.score(sv['id'], sv['svtype'], sv['chrom'],
            sv['pos'], sv['end'], sv['support'], layer_results)
        results.append(result)

    print(f"    Completed {len(results)} SVs\n")

    # FDR estimation
    all_tscores = [r.t_score for r in results]
    try:
        fdr_estimate = estimate_fdr(all_tscores)
    except Exception:
        fdr_estimate = None

    # Reports
    reports_dir = os.path.join(output_dir, 'reports')
    os.makedirs(reports_dir, exist_ok=True)
    for result in results:
        report = generate_report_card(result)
        with open(os.path.join(reports_dir, f"{result.sv_id}.txt"), 'w') as f:
            f.write(report)

    summary = generate_summary_table(results)
    with open(os.path.join(output_dir, 'validation_summary.txt'), 'w') as f:
        f.write(summary)
    print(summary)

    json_output = {'pipeline': 'VALID-SV v0.9.4', 'n_svs_validated': len(results),
                   'results': [r.to_dict() for r in results]}
    with open(os.path.join(output_dir, 'validation_results.json'), 'w') as f:
        json.dump(json_output, f, indent=2)

    print("\n" + "=" * 70)
    print("  VALIDATION COMPLETE")
    print("=" * 70 + "\n")
    return json_output


def main():
    parser = argparse.ArgumentParser(description='VALID-SV: Triangulation-based SV validation')
    parser.add_argument('--consensus-vcf', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--fastq', default=None)
    parser.add_argument('--output', default='results/validation')
    parser.add_argument('--min-support', type=int, default=1)
    parser.add_argument('--max-svs', type=int, default=None)
    parser.add_argument('--skip-kmer', action='store_true')
    parser.add_argument('--jellyfish-db', default=None)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()

    run_validation_pipeline(args.consensus_vcf, args.bam, args.reference, args.fastq,
        args.output, args.min_support, args.max_svs, args.skip_kmer,
        args.jellyfish_db, threads=args.threads)


if __name__ == '__main__':
    main()

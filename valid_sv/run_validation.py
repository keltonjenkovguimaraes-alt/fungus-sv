#!/usr/bin/env python3
"""
VALID-SV: Main Validation Pipeline
====================================
Entry point for the triangulation-based SV validation.

Usage:
    python -m valid_sv.run_validation \
        --consensus-vcf results/variants/consensus/consensus.vcf \
        --bam results/alignment/sample.sorted.bam \
        --reference data/reference/reference.fasta \
        --fastq data/raw/sample.fastq.gz \
        --output results/validation/

Architecture:
    Prediction phase (existing FUNGUS-SV):
        ICB consensus → Candidate SV list

    Validation phase (VALID-SV):
        Layer 1: ICB multi-caller agreement (REPORTED, not scored - circular)
        Layer 2: Local Assembly Refinement (LAR) — must be run separately
        Layer 3: Read-Depth Signature
        Layer 4: k-mer Spectrum Analysis
        Layer 5: Breakpoint Junction Analysis
        Layer 6: Ploidy Confirmation (SNV het rate)
        → Triangulation Engine → T-score + Confidence Estimate

Author: VALID-SV / FUNGUS-SV
Status: Development — NOT FOR PRODUCTION USE
"""

import sys
import os
import argparse
import json
from pathlib import Path
from typing import List, Dict, Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from valid_sv.evidence.layer_depth import analyze_depth_signature, DepthEvidence
from valid_sv.evidence.layer_kmer import analyze_kmer_spectrum, KmerEvidence
from valid_sv.evidence.layer_breakpoint import analyze_breakpoint_junctions, BreakpointEvidence
from valid_sv.evidence.layer_ploidy import analyze_ploidy, run_longshot, PloidyEvidence
from valid_sv.evidence.layer_lar import run_lar, LAREvidence, LARVerdict
from valid_sv.quality.triangulability import assess_triangulability, TriangulabilityReport
from valid_sv.engine.scorer import (
    TriangulationScorer, LayerResult, TriangulationResult, TScoreTier
)
from valid_sv.engine.fdr_estimator import estimate_fdr
from valid_sv.reporting.report_card import generate_report_card, generate_summary_table


def parse_consensus_vcf(vcf_path: str) -> List[dict]:
    """Parse SV calls from the ICB consensus VCF."""
    svs = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            sv_id = parts[2]
            info = parts[7]
            # Extract SV type
            import re
            svtype_match = re.search(r'SVTYPE=(\w+)', info)
            svtype = svtype_match.group(1) if svtype_match else 'UNK'

            # Extract end position
            end_match = re.search(r'END=(\d+)', info)
            svlen_match = re.search(r'SVLEN=(\d+)', info)
            if end_match:
                end = int(end_match.group(1))
            elif svlen_match:
                end = pos + abs(int(svlen_match.group(1)))
            else:
                end = pos

            # Extract support
            support_match = re.search(r'SUPPORT=(\d+)', info)
            support = int(support_match.group(1)) if support_match else 1

            svs.append({
                'id': sv_id,
                'chrom': chrom,
                'pos': pos,
                'end': end,
                'svtype': svtype,
                'size': end - pos,
                'support': support,
            })

    return svs


def run_validation_pipeline(consensus_vcf: str, bam_path: str,
                            reference_path: str, fastq_path: Optional[str],
                            output_dir: str, min_support: int = 1,
                            max_svs: Optional[int] = None,
                            skip_kmer: bool = False,
                            jellyfish_db: Optional[str] = None,
                            ablation: bool = False,
                            threads: int = 4) -> dict:
    """
    Run complete validation pipeline on consensus SVs.

    Args:
        consensus_vcf: Path to ICB consensus VCF
        bam_path: Path to aligned BAM
        reference_path: Path to reference FASTA
        fastq_path: Path to raw reads (for k-mer layer)
        output_dir: Output directory for reports
        min_support: Minimum ICB support to validate (1, 2, or 3)
        max_svs: Maximum number of SVs to validate (for testing)
        skip_kmer: Skip k-mer layer if jellyfish unavailable

    Returns:
        Dictionary with summary statistics
    """
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("  VALID-SV: Triangulation-Based SV Validation")
    print("=" * 70)
    print(f"  Consensus VCF: {consensus_vcf}")
    print(f"  BAM: {bam_path}")
    print(f"  Reference: {reference_path}")
    print(f"  Raw reads: {fastq_path or 'NOT PROVIDED (k-mer layer disabled)'}")
    print(f"  Output: {output_dir}")
    print(f"  Min ICB support: {min_support}")
    print("=" * 70)
    print()
    
    # Parse SVs
    all_svs = parse_consensus_vcf(consensus_vcf)
    print(f"  Loaded {len(all_svs)} SVs from consensus VCF")

    # Filter by support
    svs_to_validate = [sv for sv in all_svs if sv['support'] >= min_support]
    print(f"  Validating {len(svs_to_validate)} SVs with SUPPORT ≥ {min_support}")

    if max_svs:
        svs_to_validate = svs_to_validate[:max_svs]
        print(f"  (Limited to {max_svs} for testing)")

    # Assess triangulability first
    print("\n  Assessing triangulability...")
    triangulability_reports = {}
    for sv in svs_to_validate:
        report = assess_triangulability(
            sv['id'], sv['svtype'], sv['size'],
            has_reference=bool(reference_path),
            has_raw_reads=bool(fastq_path),
            has_bam=bool(bam_path)
        )
        triangulability_reports[sv['id']] = report

    fully = sum(1 for r in triangulability_reports.values()
                if r.tier.name == 'FULLY_TRIANGULABLE')
    partially = sum(1 for r in triangulability_reports.values()
                   if r.tier.name == 'PARTIALLY_TRIANGULABLE')
    limited = sum(1 for r in triangulability_reports.values()
                 if r.tier.name in ('LIMITED', 'NOT_TRIANGULABLE'))
    print(f"    Fully triangulable: {fully}")
    print(f"    Partially triangulable: {partially}")
    print(f"    Limited/not triangulable: {limited}")
    
    # Build k-mer database once (if fastq provided)
    kmer_db = None
    if jellyfish_db and os.path.exists(jellyfish_db):
        from valid_sv.evidence.layer_kmer import set_database_path
        kmer_db = jellyfish_db
        set_database_path(kmer_db)
        print(f"  Using pre-built k-mer DB: {kmer_db}")
    elif fastq_path and not skip_kmer:
        from valid_sv.evidence.layer_kmer import build_kmer_database, set_database_path
        try:
            print("\n  Building k-mer database (one-time)...")
            kmer_db = build_kmer_database(fastq_path, output_dir=os.path.join(output_dir, "kmer_db"))
            set_database_path(kmer_db)
            print(f"  k-mer DB ready: {kmer_db}")
        except Exception as e:
            print(f"  WARNING: k-mer DB build failed: {e}")
            print(f"  k-mer layer will be skipped")
            skip_kmer = True

    # Run evidence layers
    print("\n  Running evidence layers...")

    scorer = TriangulationScorer()
    results = []
    
    for i, sv in enumerate(svs_to_validate):
        if i % 10 == 0 and i > 0:
            print(f"    Processed {i}/{len(svs_to_validate)} SVs...")

        layer_results = []

        # Layer 1: Alignment consensus (pre-computed from VCF)
        icb_score = {1: 0.33, 2: 0.67, 3: 1.0}.get(sv['support'], 0.33)
        layer_results.append(LayerResult(
            "alignment_consensus", icb_score,
            f"{sv['support']}/3 callers", True, 0.0,
            f"ICB support: {sv['support']} callers (reported but NOT used in T-score)"
        ))

        # Layer 2: Local assembly (LAR)
        if args.lar:
            try:
                lar_result = run_lar(
                    bam_path, reference_path,
                    sv['id'], sv['svtype'],
                    sv['chrom'], sv['pos'], sv['end']
                )
                layer_results.append(LayerResult(
                    "local_assembly", lar_result.evidence_score,
                    lar_result.verdict.value, lar_result.evidence_score > 0,
                    0.20, lar_result.details
                ))
            except Exception as e:
                layer_results.append(LayerResult(
                    "local_assembly", 0.0, "error", False, 0.20,
                    f"LAR failed: {str(e)}"
                ))
        else:
            layer_results.append(LayerResult(
                "local_assembly", 0.0, "not_run", False, 0.20,
                "LAR not requested (use --lar flag)"
            ))

        # Layer 3: Depth signature
        depth_available = any(l.layer_name == 'depth_signature' and l.available
                              for l in triang.layers) if triang else False
        if depth_available and sv['svtype'] in ('DEL', 'DUP'):
            try:
                depth_result = analyze_depth_signature(
                    bam_path, sv['id'], sv['svtype'],
                    sv['chrom'], sv['pos'], sv['end']
                )
                layer_results.append(LayerResult(
                    "depth_signature", depth_result.evidence_score,
                    depth_result.verdict.value, True, 0.25,
                    depth_result.details
                ))
            except Exception as e:
                layer_results.append(LayerResult(
                    "depth_signature", 0.0, "error", False, 0.25,
                    f"Depth analysis failed: {str(e)}"
                ))
        else:
            layer_results.append(LayerResult(
                "depth_signature", 0.0, "not_applicable", False, 0.25,
                f"Not applicable for {sv['svtype']}"
            ))

        # Layer 4: k-mer spectrum
        kmer_available = (any(l.layer_name == 'kmer_spectrum' and l.available
                             for l in triang.layers) if triang else False)
        kmer_available = kmer_available and not skip_kmer and fastq_path

        if kmer_available and sv['svtype'] in ('DEL', 'INS'):
            try:
                kmer_result = analyze_kmer_spectrum(
                    fastq_path, reference_path,
                    sv['id'], sv['svtype'],
                    sv['chrom'], sv['pos'], sv['end'],
                    jf_db=kmer_db
                )
                layer_results.append(LayerResult(
                    "kmer_spectrum", kmer_result.evidence_score,
                    kmer_result.verdict.value, True, 0.25,
                    kmer_result.details
                ))
            except Exception as e:
                layer_results.append(LayerResult(
                    "kmer_spectrum", 0.0, "error", False, 0.25,
                    f"k-mer analysis failed: {str(e)}"
                ))
        else:
            layer_results.append(LayerResult(
                "kmer_spectrum", 0.0,
                "not_applicable" if sv['svtype'] not in ('DEL', 'INS') else "unavailable",
                False, 0.25,
                "k-mer layer not available (missing FASTQ or jellyfish)"
            ))

        # Layer 5: Breakpoint junction
        bp_available = any(l.layer_name == 'breakpoint_junction' and l.available
                          for l in triang.layers) if triang else False
        if bp_available:
            try:
                bp_result = analyze_breakpoint_junctions(
                    bam_path, sv['id'], sv['svtype'],
                    sv['chrom'], sv['pos'], sv['end']
                )
                layer_results.append(LayerResult(
                    "breakpoint_junction", bp_result.evidence_score,
                    bp_result.verdict.value, True, 0.25,
                    bp_result.details
                ))
            except Exception as e:
                layer_results.append(LayerResult(
                    "breakpoint_junction", 0.0, "error", False, 0.25,
                    f"Breakpoint analysis failed: {str(e)}"
                ))
        else:
            layer_results.append(LayerResult(
                "breakpoint_junction", 0.0, "unavailable", False, 0.25,
                "No BAM available"
            ))

        # Layer 6: Ploidy confirmation (SNV het rate)
        try:
            ploidy_vcf = os.path.join(output_dir, 'longshot_snvs.vcf')
            if os.path.exists(ploidy_vcf):
                ploidy_result = analyze_ploidy(ploidy_vcf)
            else:
                run_longshot(bam_path, reference_path, ploidy_vcf)
                ploidy_result = analyze_ploidy(ploidy_vcf)

            ploidy_score = ploidy_result.evidence_score if ploidy_result.is_haploid else 0.3
            layer_results.append(LayerResult(
                "ploidy_confirmation", ploidy_score,
                f"het_rate={ploidy_result.het_rate:.3f}", True, 0.15,
                ploidy_result.details
            ))
        except Exception as e:
            layer_results.append(LayerResult(
                "ploidy_confirmation", 0.0, "error", False, 0.00,
                f"Ploidy analysis failed: {str(e)}"
            ))
        
        # Score
        result = scorer.score(
            sv['id'], sv['svtype'], sv['chrom'],
            sv['pos'], sv['end'], sv['support'],
            layer_results
        )
        results.append(result)

    print(f"    Completed {len(results)} SVs\n")

    # Estimate FDR
    all_tscores = [r.t_score for r in results]

    print("  Estimating FDR from T-score distribution...")
    try:
        fdr_estimate = estimate_fdr(all_tscores)
        print(f"    True component mean: {fdr_estimate.true_component_mean:.3f}")
        print(f"    False component mean: {fdr_estimate.false_component_mean:.3f}")
        print(f"    Estimated true proportion: {fdr_estimate.true_component_weight:.1%}")
    except Exception:
        fdr_estimate = None
        print("    Using simple estimator (install scikit-learn for mixture model)")

    # Generate report cards
    print("\n  Generating report cards...")

    reports_dir = os.path.join(output_dir, 'reports')
    os.makedirs(reports_dir, exist_ok=True)

    # Individual reports
    for result in results:
        report = generate_report_card(result)
        report_path = os.path.join(reports_dir, f"{result.sv_id}.txt")
        with open(report_path, 'w') as f:
            f.write(report)

    # Summary table
    summary = generate_summary_table(results)
    summary_path = os.path.join(output_dir, 'validation_summary.txt')
    with open(summary_path, 'w') as f:
        f.write(summary)

    print(summary)
    print(f"\n  Individual reports: {reports_dir}/")
    print(f"  Summary: {summary_path}")

    # JSON output for programmatic use
    json_output = {
        'pipeline': 'VALID-SV v0.1.0',
        'status': 'DEVELOPMENT — estimates are approximate',
        'n_svs_validated': len(results),
        'fdr_estimate': {
            'method': 'mixture_model' if fdr_estimate else 'simple_empirical',
            'thresholds': fdr_estimate.thresholds if fdr_estimate else
                         estimate_fdr(all_tscores).thresholds,
        } if all_tscores else {},
        'results': [r.to_dict() for r in results],
    }

    json_path = os.path.join(output_dir, 'validation_results.json')
    with open(json_path, 'w') as f:
        json.dump(json_output, f, indent=2)
    print(f"  JSON output: {json_path}")

    print("\n" + "=" * 70)
    print("  VALIDATION COMPLETE")
    print("=" * 70)
    print("\n  ⚠️  IMPORTANT CAVEATS:")
    print("  1. T-scores and FDR estimates are APPROXIMATE")
    print("  2. Calibrate with synthetic benchmarks before publication")
    print("  3. Experimental validation required for biological conclusions")
    print("  4. This is a hypothesis-generation tool, not a truth machine")
    print()

    return json_output


def main():
    parser = argparse.ArgumentParser(
        description='VALID-SV: Triangulation-based SV validation',
        epilog='Development version — estimates are approximate.'
    )

    parser.add_argument('--consensus-vcf', required=True,
                       help='ICB consensus VCF file')
    parser.add_argument('--bam', required=True,
                       help='Aligned BAM file (sorted + indexed)')
    parser.add_argument('--reference', required=True,
                       help='Reference genome FASTA (indexed)')
    parser.add_argument('--fastq', default=None,
                       help='Raw PacBio HiFi reads (for k-mer layer)')
    parser.add_argument('--output', default='results/validation',
                       help='Output directory')
    parser.add_argument('--min-support', type=int, default=1,
                       help='Minimum ICB support (1-3)')
    parser.add_argument('--max-svs', type=int, default=None,
                       help='Max SVs to validate (for testing)')
    parser.add_argument('--skip-kmer', action='store_true',
                       help='Skip k-mer layer (if jellyfish unavailable)')
    parser.add_argument('--ablation', action='store_true',
                       help='Run each validation layer individually for ablation study')
    parser.add_argument('--jellyfish-db', default=None,
                       help='Pre-built jellyfish database (.jf)')
    parser.add_argument('--threads', type=int, default=4,
                       help='Threads for parallel steps')

    args = parser.parse_args()

    run_validation_pipeline(
        consensus_vcf=args.consensus_vcf,
        bam_path=args.bam,
        reference_path=args.reference,
        fastq_path=args.fastq,
        output_dir=args.output,
        min_support=args.min_support,
        max_svs=args.max_svs,
        skip_kmer=args.skip_kmer,
        ablation=args.ablation,
        jellyfish_db=args.jellyfish_db,
        threads=args.threads,
    )


if __name__ == '__main__':
    main()

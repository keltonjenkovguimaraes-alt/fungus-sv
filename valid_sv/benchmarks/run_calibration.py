#!/usr/bin/env python3
"""Complete calibration pipeline for FUNGUS-SV T-score thresholds.

Runs: spike_in → simulate reads → align → ICB → validate → evaluate.
Outputs precision/recall/F1 at different T-score thresholds.
"""

import subprocess, os, sys, json, argparse, re
from pathlib import Path


def run_cmd(cmd, desc="", timeout=3600):
    """Run a shell command with feedback."""
    if desc:
        print(f"\n  [{desc}]")
    print(f"  $ {' '.join(cmd) if isinstance(cmd, list) else cmd[:100]}")
    result = subprocess.run(cmd, capture_output=True, text=True, 
                          timeout=timeout, shell=isinstance(cmd, str))
    if result.returncode != 0:
        print(f"  WARNING: {result.stderr[:300]}")
    return result


def parse_truth(truth_vcf):
    """Parse the spike-in truth VCF."""
    truth = {}
    with open(truth_vcf) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            m = re.search(r'SVTYPE=(\w+)', parts[7])
            svtype = m.group(1) if m else 'UNK'
            truth[parts[2]] = {
                'chrom': parts[0], 'pos': int(parts[1]),
                'type': svtype
            }
    return truth


def evaluate(consensus_vcf, truth_vcf, validation_json, output_dir):
    """Compare pipeline output to truth set and compute metrics."""
    truth = parse_truth(truth_vcf)
    
    # Parse consensus SVs
    consensus = []
    with open(consensus_vcf) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            m = re.search(r'SVTYPE=(\w+)', parts[7])
            svtype = m.group(1) if m else 'UNK'
            consensus.append({
                'id': parts[2], 'chrom': parts[0],
                'pos': int(parts[1]), 'type': svtype
            })
    
    # Parse validation results
    with open(validation_json) as f:
        vdata = json.load(f)
    
    tscores = {r['sv_id']: r['t_score'] for r in vdata['results']}
    
    # Match consensus to truth (≤2kb distance, same type)
    matched = []
    for csv in consensus:
        t = tscores.get(csv['id'], 0)
        for tid, tv in truth.items():
            if (tv['chrom'] == csv['chrom'] and tv['type'] == csv['type'] 
                and abs(tv['pos'] - csv['pos']) <= 2000):
                matched.append({'id': csv['id'], 't_score': t, 'truth_id': tid})
                break
    
    n_truth = len(truth)
    n_consensus = len(consensus)
    n_tp = len(matched)
    
    # Compute metrics at different T-score thresholds
    thresholds = [0.2, 0.4, 0.6, 0.8]
    results = {'n_truth': n_truth, 'n_consensus': n_consensus, 'n_tp': n_tp}
    
    for thresh in thresholds:
        passing = [m for m in matched if m['t_score'] >= thresh]
        tp = len(passing)
        fp = sum(1 for csv in consensus 
                if tscores.get(csv['id'], 0) >= thresh 
                and not any(abs(tv['pos'] - csv['pos']) <= 2000 
                           and tv['chrom'] == csv['chrom'] 
                           and tv['type'] == csv['type'] 
                           for tv in truth.values()))
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / n_truth if n_truth > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        results[f'T≥{thresh}'] = {
            'tp': tp, 'fp': fp,
            'precision': round(precision, 4),
            'recall': round(recall, 4),
            'f1': round(f1, 4)
        }
    
    # Save results
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'calibration_results.json'), 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"  CALIBRATION RESULTS")
    print(f"{'='*60}")
    print(f"  Truth SVs: {n_truth}")
    print(f"  Consensus SVs: {n_consensus}")
    print(f"  True Positives (matched): {n_tp}")
    print(f"  Overall Precision: {n_tp/n_consensus*100:.1f}%" if n_consensus else "  N/A")
    print(f"  Overall Recall: {n_tp/n_truth*100:.1f}%" if n_truth else "  N/A")
    print(f"\n  T-score → FDR mapping:")
    for thresh in thresholds:
        r = results[f'T≥{thresh}']
        fdr = 1 - r['precision'] if r['precision'] > 0 else 1.0
        print(f"    T≥{thresh}: Precision={r['precision']:.3f}, Recall={r['recall']:.3f}, "
              f"F1={r['f1']:.3f}, est.FDR={fdr:.3f}")
    print(f"{'='*60}\n")
    
    return results


def main():
    parser = argparse.ArgumentParser(description='FUNGUS-SV Calibration Pipeline')
    parser.add_argument('--reference', required=True)
    parser.add_argument('--output', default='results/calibration')
    parser.add_argument('--n-svs', type=int, default=50)
    parser.add_argument('--coverage', type=int, default=58)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    
    print("=" * 60)
    print("  FUNGUS-SV: Calibration Pipeline")
    print("=" * 60)
    
    # Step 1: Generate spike-in SVs + simulate reads
    print("\n[1/4] Generating spike-in data...")
    run_cmd([
        'python3', 'valid_sv/benchmarks/spike_in.py',
        '--reference', args.reference,
        '--output-dir', args.output,
        '--n-svs', str(args.n_svs),
        '--coverage', str(args.coverage)
    ], "spike_in")
    
    # Step 2: Align reads
    print("\n[2/4] Aligning reads...")
    ref = f'{args.output}/modified_reference.fasta'
    reads = f'{args.output}/simulated_reads.fastq.gz'
    bam = f'{args.output}/alignment/calib.sorted.bam'
    
    if os.path.exists(reads):
        os.makedirs(f'{args.output}/alignment', exist_ok=True)
        run_cmd([
            'minimap2', '-t', str(args.threads), '-ax', 'map-hifi',
            '-R', '@RG\\tID:calib\\tSM:calib_sample', ref, reads
        ], "minimap2")
        # Note: piping to samtools sort would need shell=True
        os.system(f"minimap2 -t {args.threads} -ax map-hifi "
                 f"-R '@RG\\tID:calib\\tSM:calib_sample' {ref} {reads} "
                 f"| samtools sort -@ 4 -o {bam} - 2>/dev/null")
        os.system(f"samtools index {bam} 2>/dev/null")
    
    # Step 3: Run ICB
    print("\n[3/4] Running ICB consensus...")
    consensus_vcf = f'{args.output}/consensus_svs.vcf'
    run_cmd([
        'python3', 'fungus_sv/core/icb.py',
        '--bam', bam,
        '--reference', args.reference,
        '--output', args.output,
        '--callers', 'sniffles2', 'cutesv', 'svim',
        '--min-callers', '2',
        '--threads', str(args.threads)
    ], "ICB consensus")
    
    # Fix output path
    actual_consensus = f'{args.output}/consensus_svs.vcf'
    
    # Step 4: Run validation
    print("\n[4/4] Running validation...")
    val_dir = f'{args.output}/validation'
    truth_vcf = 'results/benchmarks/synthetic_truth.vcf'
    
    run_cmd([
        'python3', '-m', 'valid_sv.run_validation',
        '--consensus-vcf', actual_consensus,
        '--bam', bam,
        '--reference', args.reference,
        '--fastq', reads,
        '--output', val_dir,
        '--threads', str(args.threads)
    ], "validation")
    
    # Step 5: Evaluate
    val_json = f'{val_dir}/validation_results.json'
    if os.path.exists(actual_consensus) and os.path.exists(val_json):
        evaluate(actual_consensus, truth_vcf, val_json, args.output)
    else:
        print("\n  Could not find consensus VCF or validation JSON.")
        print(f"  Expected: {actual_consensus}")
        print(f"  Expected: {val_json}")


if __name__ == '__main__':
    main()

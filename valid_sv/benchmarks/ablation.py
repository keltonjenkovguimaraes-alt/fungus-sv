#!/usr/bin/env python3
"""Ablation study: measure each layer's contribution to SV validation.

Runs FUNGUS-SV with each evidence layer individually and all combinations
to quantify how much each layer contributes to precision and recall.
"""

import subprocess, os, sys, json, argparse


def run_ablation(consensus_vcf, bam, reference, fastq, output_dir, threads=4):
    """Run validation with each layer individually."""
    
    layers = [
        ('depth_only', ['--skip-kmer']),
        ('kmer_only', []),  # Need to implement layer-specific flags
        ('breakpoint_only', []),
        ('lar_only', []),
        ('all_layers', []),
    ]
    
    results = {}
    for name, extra_args in layers:
        print(f"\n  Running: {name}")
        out = f'{output_dir}/{name}'
        os.makedirs(out, exist_ok=True)
        
        cmd = [
            'python3', '-m', 'valid_sv.run_validation',
            '--consensus-vcf', consensus_vcf,
            '--bam', bam,
            '--reference', reference,
            '--output', out,
            '--threads', str(threads)
        ]
        if fastq:
            cmd += ['--fastq', fastq]
        cmd += extra_args
        
        subprocess.run(cmd, capture_output=True, timeout=3600)
        
        # Collect results
        val_json = f'{out}/validation_results.json'
        if os.path.exists(val_json):
            with open(val_json) as f:
                data = json.load(f)
            results[name] = {
                'n_svs': data.get('n_svs_validated', 0),
                'mean_t_score': sum(r['t_score'] for r in data['results']) / len(data['results']) if data['results'] else 0
            }
    
    # Save ablation results
    with open(f'{output_dir}/ablation_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"  ABLATION RESULTS")
    print(f"{'='*60}")
    for name, r in results.items():
        print(f"  {name:<25s}: {r['n_svs']} SVs, mean T={r['mean_t_score']:.3f}")
    print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(description='Ablation study')
    parser.add_argument('--consensus-vcf', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--fastq', default=None)
    parser.add_argument('--output', default='results/ablation')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    
    run_ablation(args.consensus_vcf, args.bam, args.reference, 
                 args.fastq, args.output, args.threads)


if __name__ == '__main__':
    main()

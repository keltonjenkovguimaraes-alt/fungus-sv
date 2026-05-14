#!/usr/bin/env python3
"""Spike-in Calibration for VALID-SV T-score thresholds."""

import subprocess, os, sys, random, json, argparse
from pathlib import Path


def create_spike_truth(reference_path: str, n_del: int = 20,
                       n_inv: int = 5) -> list:
    """Generate synthetic SV truth set at random positions."""
    contigs = {}
    current_contig = None
    current_len = 0
    
    with open(reference_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_contig:
                    contigs[current_contig] = current_len
                current_contig = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_contig:
            contigs[current_contig] = current_len
    
    spike_svs = []
    
    for i in range(n_del):
        contig = random.choice(list(contigs.keys()))
        contig_len = contigs[contig]
        size = random.randint(100, 10000)
        start = random.randint(1000, contig_len - size - 1000)
        spike_svs.append({
            'id': f'SPIKE_DEL_{i+1:03d}', 'type': 'DEL',
            'chrom': contig, 'start': start, 'end': start + size,
            'size': size, 'true_positive': True
        })
    
    for i in range(n_inv):
        contig = random.choice(list(contigs.keys()))
        contig_len = contigs[contig]
        size = random.randint(500, 5000)
        start = random.randint(1000, contig_len - size - 1000)
        spike_svs.append({
            'id': f'SPIKE_INV_{i+1:03d}', 'type': 'INV',
            'chrom': contig, 'start': start, 'end': start + size,
            'size': size, 'true_positive': True
        })
    
    os.makedirs('results/calibration', exist_ok=True)
    truth_path = 'results/calibration/spike_truth.json'
    with open(truth_path, 'w') as f:
        json.dump(spike_svs, f, indent=2)
    
    print(f"Created {len(spike_svs)} spike-in SVs: {n_del} DEL + {n_inv} INV")
    return spike_svs


def main():
    parser = argparse.ArgumentParser(description='Spike-in calibration')
    parser.add_argument('--reference', required=True)
    parser.add_argument('--n-del', type=int, default=20)
    parser.add_argument('--n-inv', type=int, default=5)
    parser.add_argument('--output', default='results/calibration')
    args = parser.parse_args()
    
    spike_svs = create_spike_truth(args.reference, args.n_del, args.n_inv)
    print(f"Truth set: {args.output}/spike_truth.json")
    print("\nTo complete calibration:")
    print("1. Simulate reads with spike-in SVs using pbsim3")
    print("2. Run FUNGUS-SV pipeline on simulated reads")
    print("3. Compare T-scores against known truth positions")


if __name__ == '__main__':
    main()

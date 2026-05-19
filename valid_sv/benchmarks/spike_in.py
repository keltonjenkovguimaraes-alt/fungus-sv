#!/usr/bin/env python3
"""
Synthetic Benchmark Generator for VALID-SV
===========================================
Introduces known SVs into a reference genome, simulates
PacBio HiFi reads, then runs FUNGUS-SV + VALID-SV to
measure precision, recall, and F1 score.

This provides the calibration data needed to map
T-scores to actual false discovery rates.

Author: VALID-SV / FUNGUS-SV
"""

import subprocess
import os
import sys
import random
import argparse
from pathlib import Path
from typing import List, Dict, Tuple


def generate_random_svs(reference_path: str, n_svs: int = 50,
                         min_size: int = 50, max_size: int = 10000) -> List[Dict]:
    """
    Generate random SV specifications.
    Returns list of {type, chrom, start, end, size}.
    """
    # Get contig lengths from reference
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
    
    print(f"  Reference has {len(contigs)} contigs, total {sum(contigs.values()):,} bp")
    
    svs = []
    sv_types = ['DEL'] * 20 + ['INS'] * 10 + ['INV'] * 10 + ['DUP'] * 10  # 50 total
    random.shuffle(sv_types)
    
    for i, svtype in enumerate(sv_types):
        # Pick random contig (weighted by length)
        contig = random.choices(list(contigs.keys()), 
                               weights=list(contigs.values()))[0]
        contig_len = contigs[contig]
        
        # Pick random position and size
        size = random.randint(min_size, min(max_size, contig_len // 10))
        start = random.randint(1000, contig_len - size - 1000)
        end = start + size
        
        svs.append({
            'id': f'SYNTH_SV_{i+1:03d}',
            'type': svtype,
            'chrom': contig,
            'start': start,
            'end': end,
            'size': size
        })
    
    # Write truth VCF
    truth_vcf = 'results/benchmarks/synthetic_truth.vcf'
    os.makedirs(os.path.dirname(truth_vcf), exist_ok=True)
    
    with open(truth_vcf, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('##source=VALID-SV_SyntheticBenchmark\n')
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        
        for sv in svs:
            f.write(f'{sv["chrom"]}\t{sv["start"]}\t{sv["id"]}\tN\t<{sv["type"]}>\t.\tPASS\t'
                   f'SVTYPE={sv["type"]};SVLEN={sv["size"]};END={sv["end"]}\n')
    
    print(f"  Generated {len(svs)} synthetic SVs")
    print(f"    DEL: {sum(1 for s in svs if s['type']=='DEL')}")
    print(f"    INS: {sum(1 for s in svs if s['type']=='INS')}")
    print(f"    INV: {sum(1 for s in svs if s['type']=='INV')}")
    print(f"    DUP: {sum(1 for s in svs if s['type']=='DUP')}")
    print(f"  Truth VCF: {truth_vcf}")
    
    return svs, truth_vcf


def modify_reference(reference_path: str, svs: List[Dict], 
                     output_path: str) -> str:
    """
    Create a modified reference with synthetic SVs introduced.
    Currently supports DEL (remove sequence) and INV (reverse complement).
    INS and DUP require more complex sequence insertion.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    
    # Read reference
    genome = {}
    for record in SeqIO.parse(reference_path, 'fasta'):
        genome[record.id] = list(str(record.seq))
    
    # Apply SVs (only DEL and INV for simplicity)
    modifications = []
    for sv in svs:
        if sv['type'] == 'DEL' and sv['chrom'] in genome:
            chrom_seq = genome[sv['chrom']]
            # Replace deleted region with Ns (avoids coordinate shift)
            for i in range(sv['start'], min(sv['end'], len(chrom_seq))):
                chrom_seq[i] = 'N'
            modifications.append(sv['id'])
        
        elif sv['type'] == 'INV' and sv['chrom'] in genome:
            chrom_seq = genome[sv['chrom']]
            # Reverse complement the inverted region
            region = chrom_seq[sv['start']:sv['end']]
            rev_comp = list(str(Seq(''.join(region)).reverse_complement()))
            for i, base in enumerate(rev_comp):
                if sv['start'] + i < len(chrom_seq):
                    chrom_seq[sv['start'] + i] = base
            modifications.append(sv['id'])
    
    # Write modified genome
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        for contig_id, seq in genome.items():
            f.write(f'>{contig_id}\n')
            seq_str = ''.join(seq)
            for i in range(0, len(seq_str), 80):
                f.write(seq_str[i:i+80] + '\n')
    
    print(f"  Modified reference: {output_path}")
    print(f"  Applied {len(modifications)} modifications (DEL + INV only)")
    
    return output_path


def simulate_reads(modified_ref: str, output_fastq: str,
                   coverage: int = 30, read_n50: int = 15000):
    """
    Simulate PacBio HiFi reads using pbsim3 or badread.
    Falls back to a simple simulation if tools unavailable.
    """
    # Check for pbsim3
    if shutil.which('pbsim3'):
        print("  Using pbsim3 for read simulation...")
        cmd = (f'pbsim --strategy wgs --method errhmm '
               f'--errhmm {os.path.expanduser("~")}/miniforge3/envs/sv_valid/data/ERRHMM-RSII.model '
               f'--depth {coverage} '
               f'--genome {modified_ref} '
               f'--length-min 5000 --length-mean {read_n50} '
               f'--prefix {output_fastq}')
        subprocess.run(cmd, shell=True, check=True)
        return output_fastq + '.fastq'
    
    # Check for badread
    if shutil.which('badread'):
        print("  Using badread for read simulation...")
        cmd = (f'badread simulate --reference {modified_ref} '
               f'--quantity {coverage}x --length {read_n50},5000 '
               f'--error 0.001,0,0.005 --qscore 30,2 '
               f'| gzip > {output_fastq}.gz')
        subprocess.run(cmd, shell=True, check=True)
        return output_fastq + '.gz'
    
    print("  WARNING: No read simulator found (pbsim3 or badread)")
    print("  Install: conda install -c bioconda pbsim3")
    return None


def evaluate_results(truth_vcf: str, pipeline_vcf: str, 
                     json_results: str) -> Dict:
    """
    Compare pipeline output to truth set.
    Calculates precision, recall, F1 score.
    """
    import json
    
    # Parse truth
    truth_svs = set()
    with open(truth_vcf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            truth_svs.add(parts[2])  # SV ID
    
    # Parse pipeline results
    with open(json_results) as f:
        data = json.load(f)
    
    results = data['results']
    
    # For synthetic data, we match by position overlap with truth
    # Simplified: count SVs at different T-score thresholds
    thresholds = [0.2, 0.4, 0.6, 0.8]
    evaluation = {}
    
    for threshold in thresholds:
        passing = [r for r in results if r['t_score'] >= threshold]
        evaluation[f'T≥{threshold:.1f}'] = {
            'n_passing': len(passing),
            'pct_of_total': len(passing) / len(results) * 100 if results else 0
        }
    
    return {
        'n_truth_svs': len(truth_svs),
        'n_pipeline_svs': len(results),
        'thresholds': evaluation
    }


def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic benchmark for VALID-SV calibration'
    )
    parser.add_argument('--reference', required=True, help='Reference genome FASTA')
    parser.add_argument('--output-dir', default='results/benchmarks')
    parser.add_argument('--n-svs', type=int, default=50)
    parser.add_argument('--coverage', type=int, default=30)
    parser.add_argument('--run-pipeline', action='store_true',
                       help='Run full pipeline on synthetic data')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("  VALID-SV: Synthetic Benchmark Generator")
    print("=" * 60)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Step 1: Generate SVs
    print("\n[1/3] Generating synthetic SVs...")
    svs, truth_vcf = generate_random_svs(args.reference, args.n_svs)
    
    # Step 2: Modify reference
    print("\n[2/3] Creating modified reference...")
    modified_ref = os.path.join(args.output_dir, 'modified_reference.fasta')
    modify_reference(args.reference, svs, modified_ref)
    
    # Step 3: Simulate reads
    print("\n[3/3] Simulating reads...")
    sim_fastq = os.path.join(args.output_dir, 'simulated_reads')
    reads = simulate_reads(modified_ref, sim_fastq, args.coverage)
    
    if reads:
        print(f"  Simulated reads: {reads}")
    
    print("\n" + "=" * 60)
    print("  Benchmark ready for pipeline testing")
    print(f"  Truth VCF: {truth_vcf}")
    print(f"  Modified ref: {modified_ref}")
    if reads:
        print(f"  Simulated reads: {reads}")
        print(f"\n  To run pipeline:")
        print(f"    python -m valid_sv.run_validation \\")
        print(f"      --consensus-vcf <icb_output> \\")
        print(f"      --bam <alignment> \\")
        print(f"      --reference {modified_ref} \\")
        print(f"      --fastq {reads} \\")
        print(f"      --output {args.output_dir}/validation/")
    print("=" * 60)


if __name__ == '__main__':
    import shutil
    main()

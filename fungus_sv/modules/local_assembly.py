#!/usr/bin/env python3
"""
LAR: Local Assembly Refinement
===============================
Extracts reads spanning SV breakpoints and performs local
de novo assembly to precisely resolve breakpoint coordinates.

Uses Flye for assembly, then aligns contigs back to reference
to identify exact breakpoints.

Author: FUNGUS-SV team
"""

import subprocess
import os
import sys
import tempfile
import argparse
from pathlib import Path


def extract_region_reads(bam_path, chrom, start, end, output_fastq, flank=2000):
    """
    Extract reads spanning an SV region with flanks.
    For deletions, the deleted region itself has no reads in haploid,
    so we expand to include flanking sequence.
    """
    region_start = max(1, start - flank)
    region_end = end + flank
    region = f"{chrom}:{region_start}-{region_end}"
    
    cmd = f"samtools view -b {bam_path} {region} | samtools fastq -o {output_fastq} -"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    
    # Count reads
    if os.path.exists(output_fastq):
        count = sum(1 for _ in open(output_fastq)) // 4
        return count
    return 0


def local_assemble(reads_fastq, output_prefix, threads=4, min_reads=10):
    """
    Perform local de novo assembly using Flye.
    Falls back to wtdbg2 for very small regions.
    """
    # Count reads
    read_count = sum(1 for _ in open(reads_fastq)) // 4
    
    if read_count < min_reads:
        return None, f"Insufficient reads ({read_count}, need ≥{min_reads})"
    
    try:
        result = subprocess.run([
            'flye', '--pacbio-hifi', reads_fastq,
            '--out-dir', output_prefix,
            '--threads', str(threads),
            '--min-overlap', '1000'
        ], capture_output=True, text=True, timeout=120)
        
        assembly_fasta = os.path.join(output_prefix, 'assembly.fasta')
        if os.path.exists(assembly_fasta) and os.path.getsize(assembly_fasta) > 0:
            return assembly_fasta, f"Assembled from {read_count} reads"
        else:
            return None, f"Flye produced no assembly from {read_count} reads"
    
    except subprocess.TimeoutExpired:
        return None, "Local assembly timed out"
    except FileNotFoundError:
        return None, "Flye not installed (conda install -c bioconda flye)"


def align_assembly_to_reference(assembly_fasta, reference_fasta, output_paf):
    """Align assembled contig to reference to find breakpoints."""
    result = subprocess.run([
        'minimap2', '-cx', 'asm5', reference_fasta, assembly_fasta
    ], capture_output=True, text=True, timeout=30)
    
    with open(output_paf, 'w') as f:
        f.write(result.stdout)
    
    return result.stdout


def parse_paf_for_breakpoints(paf_content, sv_type):
    """
    Parse PAF alignment to extract refined breakpoints.
    
    Returns dict with refined coordinates and identity.
    """
    lines = [l for l in paf_content.strip().split('\n') if l]
    if not lines:
        return None
    
    # Take the longest alignment
    best = None
    best_len = 0
    for line in lines:
        parts = line.split('\t')
        if len(parts) < 12:
            continue
        query_len = int(parts[1])
        if query_len > best_len:
            best_len = query_len
            best = parts
    
    if best is None:
        return None
    
    # PAF format: query_name, query_len, query_start, query_end,
    #             strand, target_name, target_len, target_start, target_end,
    #             matches, block_len, mapq
    query_name = best[0]
    query_len = int(best[1])
    q_start = int(best[2])
    q_end = int(best[3])
    strand = best[4]
    target_name = best[5]
    target_len = int(best[6])
    t_start = int(best[7])
    t_end = int(best[8])
    matches = int(best[9])
    block_len = int(best[10])
    mapq = int(best[11])
    
    identity = matches / block_len if block_len > 0 else 0
    
    return {
        'contig_name': query_name,
        'contig_len': query_len,
        'ref_contig': target_name,
        'ref_start': t_start,
        'ref_end': t_end,
        'strand': strand,
        'identity': round(identity, 4),
        'mapq': mapq,
        'aligned_len': block_len
    }


def refine_sv(vcf_line, bam_path, reference_path, output_dir, flank=2000, threads=4):
    """
    Refine a single SV call using local assembly.
    
    Args:
        vcf_line: VCF record line
        bam_path: Path to aligned BAM
        reference_path: Path to reference FASTA
        output_dir: Directory for assembly outputs
        flank: Flanking region size
        threads: CPU threads
    
    Returns:
        dict with refinement results
    """
    import re
    
    parts = vcf_line.strip().split('\t')
    if len(parts) < 8:
        return {'status': 'error', 'message': 'Invalid VCF line'}
    
    chrom = parts[0]
    pos = int(parts[1])
    sv_id = parts[2]
    info = parts[7]
    
    svtype_match = re.search(r'SVTYPE=(\w+)', info)
    svtype = svtype_match.group(1) if svtype_match else 'UNK'
    
    end_match = re.search(r'END=(\d+)', info)
    if end_match:
        end = int(end_match.group(1))
    else:
        svlen_match = re.search(r'SVLEN=(-?\d+)', info)
        end = pos + abs(int(svlen_match.group(1))) if svlen_match else pos
    
    # Create output directory for this SV
    sv_dir = os.path.join(output_dir, sv_id)
    os.makedirs(sv_dir, exist_ok=True)
    
    # Extract reads
    reads_fq = os.path.join(sv_dir, 'reads.fastq')
    read_count = extract_region_reads(bam_path, chrom, pos, end, reads_fq, flank)
    
    if read_count < 10:
        return {
            'status': 'insufficient_reads',
            'sv_id': sv_id,
            'read_count': read_count,
            'message': f'Only {read_count} reads in region'
        }
    
    # Local assembly
    assembly_dir = os.path.join(sv_dir, 'assembly')
    assembly_fasta, message = local_assemble(reads_fq, assembly_dir, threads)
    
    if assembly_fasta is None:
        return {
            'status': 'assembly_failed',
            'sv_id': sv_id,
            'read_count': read_count,
            'message': message
        }
    
    # Align to reference
    paf_path = os.path.join(sv_dir, 'alignment.paf')
    paf_content = align_assembly_to_reference(assembly_fasta, reference_path, paf_path)
    
    # Parse breakpoints
    bp_info = parse_paf_for_breakpoints(paf_content, svtype)
    
    if bp_info is None:
        return {
            'status': 'alignment_failed',
            'sv_id': sv_id,
            'read_count': read_count,
            'message': 'Could not parse alignment'
        }
    
    # Calculate size difference from ICB call
    original_size = end - pos
    refined_size = bp_info['ref_end'] - bp_info['ref_start']
    size_ratio = refined_size / original_size if original_size > 0 else 1.0
    
    return {
        'status': 'success',
        'sv_id': sv_id,
        'sv_type': svtype,
        'read_count': read_count,
        'original_size': original_size,
        'refined_size': refined_size,
        'size_ratio': round(size_ratio, 2),
        'identity': bp_info['identity'],
        'contig_len': bp_info['contig_len'],
        'refined_start': bp_info['ref_start'],
        'refined_end': bp_info['ref_end'],
        'message': f'Refined: {original_size}bp → {refined_size}bp (identity={bp_info["identity"]:.2%})'
    }


def main():
    parser = argparse.ArgumentParser(
        description='LAR: Local Assembly Refinement for SV breakpoints'
    )
    parser.add_argument('--consensus', required=True, help='ICB consensus VCF')
    parser.add_argument('--bam', required=True, help='Aligned BAM file')
    parser.add_argument('--reference', required=True, help='Reference FASTA')
    parser.add_argument('--output', required=True, help='Output refined VCF')
    parser.add_argument('--flank', type=int, default=2000, help='Flanking region size')
    parser.add_argument('--threads', type=int, default=4, help='CPU threads')
    parser.add_argument('--max-svs', type=int, default=None, help='Max SVs to refine')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("  FUNGUS-SV: Local Assembly Refinement (LAR)")
    print("=" * 60)
    
    # Read consensus VCF
    with open(args.consensus) as f:
        vcf_lines = [l for l in f if not l.startswith('#')]
    
    print(f"  Loaded {len(vcf_lines)} SVs from consensus VCF")
    
    if args.max_svs:
        vcf_lines = vcf_lines[:args.max_svs]
        print(f"  (Limited to {args.max_svs})")
    
    # Create output directory
    lar_dir = os.path.join(os.path.dirname(args.output), 'lar_assemblies')
    os.makedirs(lar_dir, exist_ok=True)
    
    results = []
    for i, line in enumerate(vcf_lines):
        if i % 10 == 0:
            print(f"  Processing {i+1}/{len(vcf_lines)}...")
        
        result = refine_sv(line, args.bam, args.reference, lar_dir,
                          args.flank, args.threads)
        results.append(result)
    
    # Summarize
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] != 'success']
    
    print(f"\n  Results:")
    print(f"    Successfully refined: {len(successful)}")
    print(f"    Failed/insufficient: {len(failed)}")
    
    if successful:
        identities = [r['identity'] for r in successful]
        size_changes = [abs(1 - r['size_ratio']) for r in successful]
        print(f"    Mean identity: {sum(identities)/len(identities):.2%}")
        print(f"    Mean size change: {sum(size_changes)/len(size_changes):.1%}")
        
        # Flag large corrections
        big_changes = [r for r in successful if abs(1 - r['size_ratio']) > 0.5]
        if big_changes:
            print(f"    Large corrections (>50% size change): {len(big_changes)}")
            for r in big_changes[:5]:
                print(f"      {r['sv_id']}: {r['original_size']}bp → {r['refined_size']}bp")
    
    print("=" * 60)


if __name__ == '__main__':
    main()

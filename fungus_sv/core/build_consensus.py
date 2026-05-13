#!/usr/bin/env python3
"""
ICB Consensus Builder
======================
Part of FUNGUS-SV: Builds high-confidence SV set from multiple callers.

Algorithm:
1. Parse VCFs from all callers
2. Cluster overlapping SVs using reciprocal overlap
3. Score each cluster by caller agreement
4. Output high-confidence consensus VCF
"""

import sys
import re
from collections import defaultdict

def parse_vcf(vcf_file, caller_name):
    """Extract SVs from a VCF file."""
    svs = []
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, sv_id, ref, alt, qual, sv_filter, info = parts[:8]
            
            # Extract SV type and length
            svtype = re.search(r'SVTYPE=(\w+)', info)
            svlen = re.search(r'SVLEN=(-?\d+)', info)
            end_match = re.search(r';END=(\d+)', info)
            
            if svtype:
                sv = {
                    'chrom': chrom,
                    'pos': int(pos),
                    'end': int(end_match.group(1)) if end_match else int(pos) + abs(int(svlen.group(1))) if svlen else int(pos),
                    'svtype': svtype.group(1),
                    'svlen': abs(int(svlen.group(1))) if svlen else 0,
                    'caller': caller_name,
                    'qual': float(qual) if qual != '.' else 0,
                    'line': line.strip()
                }
                svs.append(sv)
    return svs


def overlap(sv1, sv2, min_overlap=0.5):
    """Calculate reciprocal overlap between two SVs."""
    if sv1['chrom'] != sv2['chrom']:
        return 0.0
    if sv1['svtype'] != sv2['svtype']:
        return 0.0
    
    start1, end1 = sv1['pos'], sv1['end']
    start2, end2 = sv2['pos'], sv2['end']
    
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start >= overlap_end:
        return 0.0
    
    overlap_len = overlap_end - overlap_start
    len1 = end1 - start1
    len2 = end2 - start2
    
    reciprocal = (overlap_len / len1) * (overlap_len / len2)
    return reciprocal


def cluster_svs(all_svs, min_overlap=0.5):
    """Cluster overlapping SVs from different callers."""
    clusters = []
    used = set()
    
    for i, sv1 in enumerate(all_svs):
        if i in used:
            continue
        cluster = [sv1]
        used.add(i)
        
        for j, sv2 in enumerate(all_svs):
            if j in used:
                continue
            # Check if sv2 overlaps with any SV in the cluster
            for member in cluster:
                if overlap(member, sv2, min_overlap) > 0:
                    cluster.append(sv2)
                    used.add(j)
                    break
        
        clusters.append(cluster)
    
    return clusters


def score_cluster(cluster, min_callers=2):
    """Score a cluster based on caller agreement."""
    callers = set(sv['caller'] for sv in cluster)
    num_callers = len(callers)
    
    if num_callers < min_callers:
        return {'support': num_callers, 'callers': callers, 'consensus': False}
    
    # Calculate average position and length
    avg_pos = sum(sv['pos'] for sv in cluster) / len(cluster)
    avg_end = sum(sv['end'] for sv in cluster) / len(cluster)
    svtype = cluster[0]['svtype']
    
    return {
        'support': num_callers,
        'callers': list(callers),
        'consensus': num_callers >= min_callers,
        'chrom': cluster[0]['chrom'],
        'pos': int(avg_pos),
        'end': int(avg_end),
        'svlen': int(avg_end - avg_pos),
        'svtype': svtype,
        'num_calls': len(cluster)
    }


def write_consensus_vcf(consensus_svs, output_file, reference):
    """Write consensus SVs to VCF format."""
    with open(output_file, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('##source=FUNGUS-SV_ICB\n')
        f.write('##reference={}\n'.format(reference))
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        f.write('##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of supporting callers">\n')
        f.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Supporting callers">\n')
        f.write('##INFO=<ID=NUM_CALLS,Number=1,Type=Integer,Description="Number of total calls in cluster">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        
        for i, sv in enumerate(consensus_svs):
            f.write('{chrom}\t{pos}\tFUNGUS_SV_{i}\tN\t<{svtype}>\t.\tPASS\t'
                    'SVTYPE={svtype};SVLEN={svlen};END={end};'
                    'SUPPORT={support};CALLERS={callers};NUM_CALLS={num_calls}\n'.format(
                chrom=sv['chrom'], pos=sv['pos'], i=i+1,
                svtype=sv['svtype'], svlen=sv['svlen'], end=sv['end'],
                support=sv['support'], callers=','.join(sv['callers']),
                num_calls=sv['num_calls']
            ))


def main():
    import argparse
    parser = argparse.ArgumentParser(description='ICB Consensus Builder')
    parser.add_argument('--vcfs', nargs='+', required=True, help='Input VCF files')
    parser.add_argument('--caller-names', nargs='+', required=True, help='Caller names')
    parser.add_argument('--output', required=True, help='Output consensus VCF')
    parser.add_argument('--reference', default='reference.fasta', help='Reference genome')
    parser.add_argument('--min-callers', type=int, default=2, help='Minimum callers for consensus')
    parser.add_argument('--min-overlap', type=float, default=0.5, help='Minimum reciprocal overlap')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("  FUNGUS-SV: Iterative Consensus Builder (ICB)")
    print("=" * 60)
    print(f"  Input VCFs: {len(args.vcfs)}")
    print(f"  Callers: {', '.join(args.caller_names)}")
    print(f"  Min callers: {args.min_callers}")
    print(f"  Min overlap: {args.min_overlap}")
    print("=" * 60)
    
    # Step 1: Parse all VCFs
    all_svs = []
    for vcf, name in zip(args.vcfs, args.caller_names):
        svs = parse_vcf(vcf, name)
        all_svs.extend(svs)
        print(f"\n  {name}: {len(svs)} SVs")
    
    print(f"\n  Total raw calls: {len(all_svs)}")
    
    # Step 2: Cluster by overlap
    clusters = cluster_svs(all_svs, args.min_overlap)
    print(f"  Total clusters: {len(clusters)}")
    
    # Step 3: Score clusters
    consensus = []
    high_conf = 0
    medium_conf = 0
    low_conf = 0
    
    for cluster in clusters:
        score = score_cluster(cluster, args.min_callers)
        if score['consensus']:
            consensus.append(score)
            if score['support'] == 3:
                high_conf += 1
            elif score['support'] == 2:
                medium_conf += 1
        else:
            low_conf += 1
    
    print(f"\n  High confidence (3 callers): {high_conf}")
    print(f"  Medium confidence (2 callers): {medium_conf}")
    print(f"  Low confidence (1 caller): {low_conf}")
    print(f"  Total consensus: {len(consensus)}")
    
    # Step 4: Write output
    write_consensus_vcf(consensus, args.output, args.reference)
    print(f"\n  Consensus VCF written to: {args.output}")
    print("=" * 60)


if __name__ == '__main__':
    main()

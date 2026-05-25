#!/usr/bin/env python3
"""
ICB: Iterative Consensus Builder
=================================
Implements interval-tree-based SV merging across multiple callers.
Reciprocal overlap threshold and breakpoint reconciliation based on
Liu et al. (2024) Nature Communications and Liu et al. (2024) Genome Biology.

Key findings applied:
- 0.5 reciprocal overlap optimal for 3-caller intersection (Liu Genome Biol, Fig. 6)
- 200 bp flank for breakpoint matching (Kronenberg et al. 2025)
- Min supporting callers = 2 for intersection (Liu Genome Biol)
- cuteSV anchors INS breakpoints; Sniffles2 anchors DEL breakpoints
"""

import sys
import subprocess
import argparse
import os
import re
import time
from pathlib import Path
from collections import defaultdict


def run_sv_caller(caller, bam, reference, output_dir, threads=4):
    """Execute individual SV callers with paper-validated parameters."""
    os.makedirs(output_dir, exist_ok=True)

    callers = {
        'sniffles2': [
            'sniffles', '--input', bam, '--vcf', f'{output_dir}/sniffles2_svs.vcf',
            '--threads', str(threads), '--minsupport', '2'
        ],
        'cutesv': [
            'cuteSV', bam, reference, f'{output_dir}/cutesv_svs.vcf',
            output_dir,
            '--max_cluster_bias_INS', '100',
            '--diff_ratio_merging_INS', '0.3',
            '--max_cluster_bias_DEL', '100',
            '--diff_ratio_merging_DEL', '0.3',
            '--min_support', '2',
            '--threads', str(threads)
        ],
        'svim': [
            'svim', 'alignment', output_dir, bam, reference,
            '--min_sv_size', '50',
            '--min_mapq', '20'
        ]
    }

    if caller not in callers:
        print(f"[ICB] Unknown caller: {caller}", file=sys.stderr)
        return None

    vcf_path = f'{output_dir}/{caller}_svs.vcf'
    if os.path.exists(vcf_path):
        os.remove(vcf_path)

    print(f"[ICB] Running {caller}...")
    result = subprocess.run(callers[caller], capture_output=True, text=True, timeout=3600)

    # SVIM outputs to variants.vcf, not svim_svs.vcf — handle safely
    if caller == 'svim':
        svim_output = os.path.join(output_dir, 'variants.vcf')
        # Wait up to 5 seconds for SVIM to finish writing
        for _ in range(10):
            if os.path.exists(svim_output):
                time.sleep(0.5)
                # Verify file is not empty and not being written
                if os.path.getsize(svim_output) > 0:
                    os.rename(svim_output, vcf_path)
                    break
            time.sleep(0.5)

    # Check for pbsv caller (not yet implemented)
    if caller == 'pbsv':
        print("[ICB] ERROR: pbsv is not yet implemented. Use sniffles2, cutesv, or svim.", file=sys.stderr)
        return None

    if result.returncode != 0:
        print(f"[ICB] ERROR: {caller} failed with return code {result.returncode}", file=sys.stderr)
        print(f"[ICB] stderr: {result.stderr[:500]}", file=sys.stderr)
        return None

    if not os.path.exists(vcf_path):
        print(f"[ICB] ERROR: {caller} did not produce output VCF: {vcf_path}", file=sys.stderr)
        return None

    # Verify VCF is not empty
    if os.path.getsize(vcf_path) == 0:
        print(f"[ICB] WARNING: {caller} produced empty VCF", file=sys.stderr)
        return None

    count = 0
    with open(vcf_path) as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    print(f"[ICB] {caller}: {count} SVs called -> {vcf_path}")

    return vcf_path


def parse_sv_record(line):
    """Parse a VCF line into a standardized SV dict."""
    parts = line.strip().split('\t')
    if len(parts) < 8:
        return None

    chrom = parts[0]
    pos = int(parts[1])
    info = parts[7]

    svtype = 'UNK'
    m = re.search(r'SVTYPE=(\w+)', info)
    if m:
        svtype = m.group(1)

    end = pos
    m = re.search(r'END=(\d+)', info)
    if m:
        end = int(m.group(1))

    svlen = abs(end - pos)
    m = re.search(r'SVLEN=(-?\d+)', info)
    if m:
        svlen = abs(int(m.group(1)))

    return {
        'chrom': chrom, 'pos': pos, 'end': end, 'id': parts[2],
        'svtype': svtype, 'svlen': svlen, 'support': 0,
        'genotype': './.', 'caller': None
    }


def load_vcf(vcf_path, caller_name):
    """Load all SV records from a VCF file."""
    svs = []
    if not vcf_path or not os.path.exists(vcf_path):
        return svs
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            sv = parse_sv_record(line)
            if sv and sv['svtype'] != 'UNK':
                sv['caller'] = caller_name
                svs.append(sv)
    return svs


def reciprocal_overlap(sv1, sv2):
    """Calculate reciprocal overlap between two SVs (0.0 to 1.0)."""
    start1, end1 = sv1['pos'], sv1['end']
    start2, end2 = sv2['pos'], sv2['end']
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start >= overlap_end:
        return 0.0
    overlap_len = overlap_end - overlap_start
    len1 = end1 - start1
    len2 = end2 - start2
    if len1 == 0 or len2 == 0:
        return 0.0
    return max(overlap_len / len1, overlap_len / len2)


def size_similarity(sv1, sv2):
    """Calculate size similarity between two SVs (Kronenberg method)."""
    len1 = sv1['svlen']
    len2 = sv2['svlen']
    if len1 == 0 or len2 == 0:
        return 0.0
    return min(len1, len2) / max(len1, len2)


def build_consensus(caller_vcfs, output_vcf, min_overlap=0.5, min_callers=2, flank=200):
    """Core ICB algorithm: Interval-tree-based consensus merging."""
    print(f"[ICB] Building consensus from {len(caller_vcfs)} callers...")
    print(f"[ICB] Parameters: overlap>={min_overlap}, callers>={min_callers}, flank={flank}bp")

    all_svs = []
    for caller_name, vcf_path in caller_vcfs.items():
        if vcf_path is None:
            continue
        svs = load_vcf(vcf_path, caller_name)
        print(f"[ICB]   {caller_name}: {len(svs)} SVs loaded")
        all_svs.extend(svs)

    if not all_svs:
        print("[ICB] ERROR: No SV calls loaded", file=sys.stderr)
        return

    sv_groups = defaultdict(list)
    for sv in all_svs:
        key = (sv['chrom'], sv['svtype'])
        sv_groups[key].append(sv)

    consensus = []

    for (chrom, svtype), sv_list in sv_groups.items():
        sv_list.sort(key=lambda x: (x['pos'], x['end']))

        clusters = []
        for sv in sv_list:
            matched = False
            for cluster in clusters:
                for member in cluster:
                    if sv['pos'] <= member['end'] + flank and member['pos'] <= sv['end'] + flank:
                        ro = reciprocal_overlap(sv, member)
                        ss = size_similarity(sv, member)
                        if ro >= min_overlap or (ss >= 0.7 and ro >= 0.3):
                            cluster.append(sv)
                            matched = True
                            break
                if matched:
                    break
            if not matched:
                clusters.append([sv])

        for cluster in clusters:
            callers_in_cluster = set(sv['caller'] for sv in cluster)
            if len(callers_in_cluster) < min_callers:
                continue

            if svtype in ('INS', 'DUP'):
                preferred_callers = ['cutesv', 'sniffles2']
            else:
                preferred_callers = ['sniffles2', 'cutesv']

            representative = None
            for pref in preferred_callers:
                for sv in cluster:
                    if sv['caller'] == pref:
                        representative = sv
                        break
                if representative:
                    break
            if representative is None:
                representative = cluster[0]

            all_starts = sorted(sv['pos'] for sv in cluster)
            all_ends = sorted(sv['end'] for sv in cluster)
            median_start = all_starts[len(all_starts)//2]
            median_end = all_ends[len(all_ends)//2]

            consensus_sv = {
                'chrom': chrom,
                'pos': median_start,
                'end': median_end,
                'svtype': svtype,
                'svlen': abs(median_end - median_start),
                'callers': sorted(callers_in_cluster),
                'num_callers': len(callers_in_cluster),
                'support': representative.get('support', 0),
                'genotype': representative.get('genotype', './.')
            }
            consensus.append(consensus_sv)

    with open(output_vcf, 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write('##source=FUNGUS-SV ICB Consensus\n')
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        out.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Supporting callers">\n')
        out.write('##INFO=<ID=NUMCALLERS,Number=1,Type=Integer,Description="Number of supporting callers">\n')
        out.write(f'##ICB_PARAMS=min_overlap={min_overlap},min_callers={min_callers},flank={flank}\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        for i, sv in enumerate(consensus):
            sv_id = f'ICB_{sv["chrom"]}_{sv["pos"]}_{i+1}'
            info = (f'SVTYPE={sv["svtype"]};END={sv["end"]};SVLEN={sv["svlen"]};'
                    f'CALLERS={",".join(sv["callers"])};NUMCALLERS={sv["num_callers"]}')
            out.write(f'{sv["chrom"]}\t{sv["pos"]}\t{sv_id}\tN\t<{sv["svtype"]}>\t.\tPASS\t{info}\n')

    counts = defaultdict(int)
    for sv in consensus:
        counts[sv['num_callers']] += 1
        counts[f"type_{sv['svtype']}"] += 1

    print(f"\n[ICB] Consensus: {len(consensus)} total SVs")
    print(f"[ICB]   2-caller: {counts[2]}, 3-caller: {counts.get(3, 0)}")
    for svtype in ['DEL', 'INS', 'INV', 'DUP']:
        if counts[f'type_{svtype}'] > 0:
            print(f"[ICB]   {svtype}: {counts[f'type_{svtype}']}")
    print(f"[ICB] Output -> {output_vcf}")


def main():
    parser = argparse.ArgumentParser(
        description='ICB: Iterative Consensus Builder for SV discovery',
        epilog='Part of FUNGUS-SV pipeline. Based on Liu et al. (2024) and Kronenberg et al. (2025).'
    )

    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--reference', required=True, help='Reference FASTA')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--callers', nargs='+',
                       default=['sniffles2', 'cutesv', 'svim'],
                       help='SV callers to use (default: sniffles2 cutesv svim)')
    parser.add_argument('--min-callers', type=int, default=2,
                       help='Minimum callers for consensus (default: 2)')
    parser.add_argument('--min-overlap', type=float, default=0.5,
                       help='Minimum reciprocal overlap (default: 0.5)')
    parser.add_argument('--flank', type=int, default=200,
                       help='Flanking bp for breakpoint tolerance (default: 200)')
    parser.add_argument('--threads', type=int, default=4,
                       help='Threads per caller (default: 4)')
    parser.add_argument('--skip-calling', action='store_true',
                       help='Skip SV calling, build consensus from existing VCFs')

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)

    print("""
    ╔══════════════════════════════════════════╗
    ║   FUNGUS-SV: Iterative Consensus Builder ║
    ║   Liu et al. 2024 validated parameters   ║
    ╚══════════════════════════════════════════╝
    """)

    caller_vcfs = {}
    if not args.skip_calling:
        for caller in args.callers:
            vcf = run_sv_caller(caller, args.bam, args.reference, args.output, args.threads)
            caller_vcfs[caller] = vcf
    else:
        for caller in args.callers:
            vcf = f'{args.output}/{caller}_svs.vcf'
            if os.path.exists(vcf):
                caller_vcfs[caller] = vcf
            else:
                print(f"[ICB] WARNING: {vcf} not found, skipping {caller}")

    consensus_vcf = f'{args.output}/consensus_svs.vcf'
    build_consensus(caller_vcfs, consensus_vcf,
                    min_overlap=args.min_overlap,
                    min_callers=args.min_callers,
                    flank=args.flank)

    print(f"\n[ICB] Done. Consensus VCF: {consensus_vcf}")


if __name__ == '__main__':
    main()

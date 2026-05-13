#!/usr/bin/env python3
"""
FUNGUS-SV Annotator
====================
Custom SV annotation for non-model organisms.
Annotates SVs with overlapping genes from GFF.
"""

import re
from collections import Counter

def parse_gff(gff_file):
    """Parse GFF and extract gene features with products."""
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feat_type, start, end = parts[0], parts[1], parts[2], int(parts[3]), int(parts[4])
            score, strand, phase, attrs = parts[5], parts[6], parts[7], parts[8]
            
            if feat_type in ['gene', 'CDS', 'exon', 'mRNA', 'tRNA', 'rRNA']:
                gene_id = re.search(r'ID=([^;]+)', attrs)
                gene_name = re.search(r'Name=([^;]+)', attrs)
                locus_tag = re.search(r'locus_tag=([^;]+)', attrs)
                product = re.search(r'product=([^;]+)', attrs)
                
                genes.append({
                    'chrom': chrom, 'start': start, 'end': end,
                    'type': feat_type, 'strand': strand,
                    'id': gene_id.group(1) if gene_id else 'unknown',
                    'name': gene_name.group(1) if gene_name else (locus_tag.group(1) if locus_tag else 'unknown'),
                    'product': product.group(1) if product else '?'
                })
    return genes


def annotate_sv(sv, genes, min_overlap=10):
    """Find genes overlapping an SV."""
    results = []
    
    for gene in genes:
        if gene['chrom'] != sv['chrom']:
            continue
        
        overlap_start = max(sv['pos'], gene['start'])
        overlap_end = min(sv['end'], gene['end'])
        
        if overlap_start >= overlap_end:
            continue
        
        overlap_len = overlap_end - overlap_start
        gene_len = gene['end'] - gene['start']
        overlap_pct = (overlap_len / gene_len) * 100 if gene_len > 0 else 0
        
        if overlap_pct < min_overlap:
            continue
        
        # Predict impact
        if gene['type'] in ['CDS', 'exon']:
            impact = 'HIGH' if overlap_pct > 50 else 'MODERATE'
        elif gene['type'] in ['mRNA', 'tRNA', 'rRNA']:
            impact = 'HIGH' if overlap_pct > 50 else 'MODERATE'
        else:
            impact = 'HIGH' if overlap_pct > 80 else ('MODERATE' if overlap_pct > 20 else 'LOW')
        
        results.append({**gene, 'overlap_pct': round(overlap_pct, 1), 'impact': impact})
    
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description='FUNGUS-SV Annotator')
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--gff', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-overlap', type=float, default=10)
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("  FUNGUS-SV: Custom SV Annotator")
    print("=" * 70)
    
    genes = parse_gff(args.gff)
    print(f"\n  Loaded {len(genes)} gene features from GFF")
    
    types = Counter(g['type'] for g in genes)
    for t, c in types.most_common():
        print(f"    {t}: {c}")
    
    sv_count, annotated, high, moderate, low = 0, 0, 0, 0, 0
    
    with open(args.vcf) as vcf, open(args.output, 'w') as out:
        out.write("SV_ID\tType\tChrom\tStart\tEnd\tLen\tSupport\t"
                 "Gene_ID\tGene_Name\tFeature\tProduct\tOverlap%\tImpact\tStrand\n")
        
        for line in vcf:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            chrom, pos, sv_id = parts[0], int(parts[1]), parts[2]
            info = parts[7]
            
            svtype = re.search(r'SVTYPE=(\w+)', info)
            svtype = svtype.group(1) if svtype else 'UNK'
            
            end_m = re.search(r'END=(\d+)', info)
            svlen_m = re.search(r'SVLEN=(\d+)', info)
            
            if end_m:
                end = int(end_m.group(1))
            elif svlen_m:
                end = pos + int(svlen_m.group(1))
            else:
                end = pos
            
            svlen = end - pos
            support = re.search(r'SUPPORT=(\d+)', info)
            support = support.group(1) if support else '0'
            
            sv_count += 1
            sv = {'chrom': chrom, 'pos': pos, 'end': end}
            
            anns = annotate_sv(sv, genes, args.min_overlap)
            
            if anns:
                for a in anns:
                    out.write(f"{sv_id}\t{svtype}\t{chrom}\t{pos}\t{end}\t{svlen}\t{support}\t"
                             f"{a['id']}\t{a['name']}\t{a['type']}\t{a['product']}\t"
                             f"{a['overlap_pct']}\t{a['impact']}\t{a['strand']}\n")
                    annotated += 1
                    if a['impact'] == 'HIGH':
                        high += 1
                    elif a['impact'] == 'MODERATE':
                        moderate += 1
                    else:
                        low += 1
    
    print(f"\n  Total SVs analyzed: {sv_count}")
    print(f"  SVs with gene overlap ≥{args.min_overlap}%: {annotated}")
    print(f"    HIGH impact: {high}")
    print(f"    MODERATE impact: {moderate}")
    print(f"    LOW impact: {low}")
    print(f"\n  Output written to: {args.output}")
    print("=" * 70)

if __name__ == '__main__':
    main()

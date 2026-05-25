#!/usr/bin/env python3
"""Annotate SVs with overlapping genes from GFF or GBFF file."""
import sys
import re
import argparse

def parse_gff(gff_path):
    """Parse GFF3 and return gene regions."""
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'gene': continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]
            gene_id = re.search(r'ID=([^;]+)', attrs)
            gene_name = re.search(r'Name=([^;]+)', attrs)
            locus = re.search(r'locus_tag=([^;]+)', attrs)
            name = gene_name.group(1) if gene_name else (locus.group(1) if locus else (gene_id.group(1) if gene_id else 'unknown'))
            genes.append({'chrom': chrom, 'start': start, 'end': end, 'name': name})
    return genes

def parse_gbff(gbff_path):
    """Parse GenBank flat file (.gbff) and return gene regions."""
    genes = []
    current_chrom = None
    current_gene = None
    
    with open(gbff_path) as f:
        for line in f:
            # Detect LOCUS line for chromosome name
            if line.startswith('LOCUS'):
                parts = line.split()
                if len(parts) >= 2:
                    current_chrom = parts[1]
                continue
            
            # Detect gene feature
            if line.startswith('     gene            '):
                if current_gene and current_chrom:
                    genes.append(current_gene)
                coords = line.strip().replace('gene', '').strip()
                # Handle complement() and join()
                coords = coords.replace('complement(', '').replace('join(', '').replace(')', '')
                if '..' in coords:
                    parts = coords.split('..')
                    try:
                        current_gene = {
                            'chrom': current_chrom,
                            'start': int(parts[0]),
                            'end': int(parts[1]),
                            'name': 'unknown'
                        }
                    except:
                        current_gene = None
                else:
                    current_gene = None
                continue
            
            # Extract locus_tag or gene name
            if current_gene:
                m = re.search(r'/locus_tag="([^"]+)"', line)
                if m:
                    current_gene['name'] = m.group(1)
                m = re.search(r'/gene="([^"]+)"', line)
                if m:
                    current_gene['name'] = m.group(1)
    
    # Don't forget last gene
    if current_gene and current_chrom:
        genes.append(current_gene)
    
    return genes

def annotate_svs(vcf_path, genes, output_path):
    """Find genes overlapping each SV."""
    svs = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            info = parts[7]
            m = re.search(r'SVTYPE=(\w+)', info)
            svtype = m.group(1) if m else 'UNK'
            m2 = re.search(r'END=(\d+)', info)
            end = int(m2.group(1)) if m2 else pos
            svs.append({'chrom': chrom, 'pos': pos, 'end': end, 'type': svtype, 'id': parts[2]})
    
    with open(output_path, 'w') as out:
        out.write("SV_ID\tChrom\tStart\tEnd\tType\tSize_bp\tGenes_Affected\n")
        for sv in svs:
            affected = []
            for gene in genes:
                # Match chromosome (try both LR813585 and chromosome 1 formats)
                if gene['chrom'] in sv['chrom'] or sv['chrom'] in gene['chrom']:
                    if max(sv['pos'], gene['start']) <= min(sv['end'], gene['end']):
                        affected.append(gene['name'])
            
            sv_size = abs(sv['end'] - sv['pos'])
            genes_str = ', '.join(affected[:10]) if affected else '(intergenic)'
            if len(affected) > 10: genes_str += f' ... +{len(affected)-10} more'
            out.write(f"{sv['id']}\t{sv['chrom']}\t{sv['pos']}\t{sv['end']}\t{sv['type']}\t{sv_size}\t{genes_str}\n")
    
    total = len(svs)
    genic = sum(1 for sv in svs if any(
        max(sv['pos'], g['start']) <= min(sv['end'], g['end'])
        for g in genes if g['chrom'] in sv['chrom'] or sv['chrom'] in g['chrom']
    ))
    intergenic = total - genic
    print(f"\nAnnotation Summary:")
    print(f"  Total SVs: {total}")
    print(f"  Genic: {genic} ({genic/total*100:.1f}%)" if total else "  No SVs")
    print(f"  Intergenic: {intergenic} ({intergenic/total*100:.1f}%)" if total else "")
    print(f"  Output: {output_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--gff', required=True)
    parser.add_argument('--output', default='sv_annotations.tsv')
    args = parser.parse_args()
    
    print(f"Loading: {args.gff}")
    
    if args.gff.endswith('.gbff') or args.gff.endswith('.gbk'):
        genes = parse_gbff(args.gff)
    else:
        genes = parse_gff(args.gff)
    
    print(f"  {len(genes)} genes loaded")
    annotate_svs(args.vcf, genes, args.output)

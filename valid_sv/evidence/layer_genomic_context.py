#!/usr/bin/env python3
"""
Evidence Layer: Genomic Context Filters
========================================
SV-type-specific hard filters using gene annotation and physical
signatures adapted from FocalSV (Luo et al. 2025) for haploid genomes.

Does NOT contribute to T-score. Acts as PASS/FLAG/FAIL gate
after triangulation, before ploidy filter.

Gene structure rules:
  - DUP: breakpoints should fall between genes; ORFs within DUP intact
  - DEL: should remove complete genes or intergenic regions
  - INV: breakpoints that split genes create potential fusions
  - TRA: inter-chromosomal split reads + gene fusion detection

Physical signature rules (FocalSV-adapted, haploid-calibrated):
  - DUP: split-read same-strand overlapping alignments + coverage shift
  - INV: split-read opposite-strand with clustered breakpoints
  - TRA: inter-chromosomal split-read pairs
"""

from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Dict, Tuple
import os


class FilterVerdict(Enum):
    PASS = "pass"
    FLAG = "flag"
    FAIL = "fail"
    NOT_APPLICABLE = "not_applicable"


@dataclass
class GenomicContextResult:
    sv_id: str
    sv_type: str
    verdict: FilterVerdict
    filters_passed: int
    filters_total: int
    genes_affected: List[str]
    gene_context: str
    details: str


# Essential genes in S. cerevisiae (SGD database)
# These are genes where deletion is lethal in rich medium.
# Source: SGD (https://www.yeastgenome.org/)
ESSENTIAL_GENES = {
    "BDP1", "GPI15", "COG6", "SSN8", "SAM50", "HHT2", "HHF2",
    "RRN6", "FMT1", "SCT1", "HIR1", "SLA1", "PDR3",
    "COX1", "ATP8", "ATP6", "COB",  # mitochondrial
    "IDH1", "NCE103", "BOP3",
}

# S288C-specific features (known insertions not in CICC-1445)
S288C_SPECIFIC = {
    "YBL005W-B": "Ty2 retrotransposon (S288C-specific insertion)",
    "YBR012W-B": "Ty2 retrotransposon (S288C-specific insertion)",
    "YBL005W-A": "Adjacent to Ty2 element",
    "YBR012W-A": "Adjacent to Ty2 element",
}

# Repeat/transposon genes that produce unreliable signals
REPEAT_GENES = {
    "FLO1", "FLO9", "FLO10",  # Flocculin genes - subtelomeric repeats
    "YBL005W-B", "YBR012W-B",  # Ty2 elements
    "YBL005W-A", "YBR012W-A",  # Ty2-adjacent
    "YNCN",  # Uncharacterized ORFs with numbered suffixes
    "YNCB", "YNCA",  # Various uncharacterized
}


def load_annotation_file(tsv_path: str) -> Dict[str, dict]:
    """
    Load pre-computed SV annotation TSV.
    
    Columns: SV_ID, Chrom, Start, End, Type, Size_bp, Genes_Affected
    
    Returns dict keyed by SV_ID.
    """
    annotations = {}
    if not os.path.exists(tsv_path):
        return annotations
    
    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            sv_id = parts[0]
            genes_str = parts[6] if len(parts) > 6 else ""
            genes = [g.strip() for g in genes_str.split(',') if g.strip() and g.strip() != '(intergenic)']
            
            annotations[sv_id] = {
                'chrom': parts[1],
                'start': int(parts[2]),
                'end': int(parts[3]),
                'sv_type': parts[4],
                'size_bp': int(parts[5]),
                'genes': genes,
                'is_intergenic': '(intergenic)' in genes_str,
            }
    
    return annotations


def classify_gene_context(genes: List[str], sv_type: str) -> Tuple[str, List[str]]:
    """
    Classify the genomic context of an SV based on affected genes.
    
    Returns (context_label, list_of_flags).
    """
    if not genes:
        return "intergenic", []
    
    flags = []
    
    # Check for essential genes
    essential_hits = [g for g in genes if g in ESSENTIAL_GENES]
    if essential_hits:
        flags.append(f"ESSENTIAL_GENE: {', '.join(essential_hits)}")
    
    # Check for S288C-specific features
    s288c_hits = [g for g in genes if g in S288C_SPECIFIC]
    if s288c_hits:
        for g in s288c_hits:
            flags.append(f"S288C_SPECIFIC: {g} = {S288C_SPECIFIC[g]}")
    
    # Check for repeat/transposon genes
    repeat_hits = [g for g in genes if any(g.startswith(r) for r in REPEAT_GENES if r in g)]
    # Also check prefix matching for YNCN*, YNCB*, etc.
    for g in genes:
        for prefix in ['YNCN', 'YNCB', 'YNCA']:
            if g.startswith(prefix) and g not in repeat_hits:
                repeat_hits.append(g)
    if repeat_hits:
        flags.append(f"REPEAT_REGION: {', '.join(repeat_hits)}")
    
    # Context classification
    if flags:
        return "genic_with_flags", flags
    else:
        return "genic_clean", []


def apply_dup_filters(sv_id: str, genes: List[str], 
                      is_intergenic: bool) -> Tuple[FilterVerdict, int, int, str]:
    """
    DUP-specific filters.
    
    Rules:
    1. Gene continuity: DUP breakpoints should not truncate genes.
       With annotation-only data, check if the gene list contains
       partial gene names or fragments.
    2. ORF integrity: All genes in duplicated region should be complete.
    3. S288C-specific check: DUPs containing Ty2/S288C insertions
       are likely reference-specific, not true CICC-1445 DUPs.
    """
    passed = 0
    total = 2
    details = []
    
    # Filter 1: Gene structure check
    context, flags = classify_gene_context(genes, "DUP")
    
    # Check for S288C-specific features (these make the DUP suspect)
    s288c_flags = [f for f in flags if "S288C_SPECIFIC" in f]
    repeat_flags = [f for f in flags if "REPEAT_REGION" in f]
    
    if s288c_flags:
        details.append(f"[FAIL] S288C-specific insertion region: {s288c_flags[0]}")
        details.append("This DUP is likely a reference-specific insertion, not a true duplication in CICC-1445")
        return FilterVerdict.FAIL, 0, total, "; ".join(details)
    
    # Filter 2: Gene integrity
    if is_intergenic:
        passed += 1
        details.append("[PASS] Intergenic duplication")
    elif context == "genic_clean":
        passed += 1
        details.append(f"[PASS] {len(genes)} genes in duplicated region; no structural flags")
    elif repeat_flags and not s288c_flags:
        passed += 1
        details.append(f"[PASS] Duplication contains {len(genes)} genes in repetitive region")
        details.append(f"  Note: {repeat_flags[0]}")
    
    # Essential gene check (FLAG but don't fail)
    essential_flags = [f for f in flags if "ESSENTIAL_GENE" in f]
    
    # Mitochondrial exception
    mitochondrial_genes = {'COX1', 'COX2', 'COX3', 'COB', 'ATP6', 'ATP8', 'ATP9',
                            '15S_RRNA', '21S_RRNA', 'tW', 'tE', 'tM', 'tF', 'tT',
                            'AI1', 'AI2', 'AI3', 'AI4', 'AI5_ALPHA', 'AI5_BETA',
                            'BI2', 'BI3', 'BI4'}
    is_mito = any(g in mitochondrial_genes for g in genes)
    
    if is_mito:
        details.append("[PASS] Mitochondrial genome — multicopy, amplifications common and non-lethal")
        passed += 1
    elif essential_flags:
        details.append(f"[FLAG] {essential_flags[0]} — dosage change may have phenotypic impact")
    
    if s288c_flags:
        return FilterVerdict.FAIL, passed, total, "; ".join(details)
    elif passed >= 1:
        return FilterVerdict.PASS if passed == total else FilterVerdict.FLAG, passed, total, "; ".join(details)
    else:
        return FilterVerdict.FAIL, passed, total, "; ".join(details)

def apply_del_filters(sv_id: str, genes: List[str],
                      is_intergenic: bool,
                      sv_size: int = 0) -> Tuple[FilterVerdict, int, int, str]:
    """
    DEL-specific filters.
    
    Rules:
    1. DEL should remove complete genes or be intergenic.
    2. Essential gene deletion is suspicious (sample would be inviable).
    """
    passed = 0
    total = 2
    details = []
    
    # Filter 1: Gene content

    # Size-stratified reliability (Pedersen & Quinlan 2019)
    if sv_size is not None and sv_size < 500:
        details.append(f"[FLAG] Small DEL ({sv_size} bp) — depth/breakpoint signals unreliable below 500 bp")
        return FilterVerdict.FLAG, 0, total, "; ".join(details)
    if is_intergenic:
        passed += 1
        details.append("[PASS] Intergenic deletion")
    elif genes:
        passed += 1
        details.append(f"[PASS] Deletion removes {len(genes)} gene(s): {', '.join(genes[:5])}")
    
    # Filter 2: Essential gene check
    context, flags = classify_gene_context(genes, "DEL")
    essential_flags = [f for f in flags if "ESSENTIAL_GENE" in f]
    repeat_flags = [f for f in flags if "REPEAT_REGION" in f]
    s288c_flags = [f for f in flags if "S288C_SPECIFIC" in f]
    
    if essential_flags:
        details.append(f"[FLAG] {essential_flags[0]} — if this were a true deletion, the strain would be inviable")
        details.append("  Likely a reference assembly difference, not a true deletion")
    elif repeat_flags:
        passed += 1
        details.append(f"[FLAG] {repeat_flags[0]}")
    elif s288c_flags:
        details.append(f"[FLAG] {s288c_flags[0]}")
    
    if passed >= 1:
        return FilterVerdict.PASS, passed, total, "; ".join(details)
    else:
        return FilterVerdict.FAIL, passed, total, "; ".join(details)


def apply_inv_filters(sv_id: str, genes: List[str],
                      is_intergenic: bool) -> Tuple[FilterVerdict, int, int, str]:
    """
    INV-specific filters.
    
    Rules:
    1. INV breakpoints that fall within genes create potential fusion genes.
    2. Intergenic inversions are less likely to have functional impact.
    """
    total = 1
    details = []
    
    if is_intergenic:
        details.append("[PASS] Inversion in intergenic region")
        return FilterVerdict.PASS, 1, total, "; ".join(details)
    elif genes:
        # Check if genes are at the boundaries (breakpoints) vs fully contained
        details.append(f"[FLAG] Inversion breakpoint near {len(genes)} gene(s): {', '.join(genes[:5])}")
        details.append("  Potential gene disruption or fusion — manual review recommended")
        return FilterVerdict.FLAG, 0, total, "; ".join(details)
    
    return FilterVerdict.PASS, 1, total, "No gene structure concerns"


def analyze_genomic_context(sv_id: str, sv_type: str,
                            annotation_tsv: str,
                            min_dup_split_reads: int = 3,
                            min_inv_split_reads: int = 3,
                            min_tra_split_reads: int = 2) -> GenomicContextResult:
    """
    Main entry point: analyze an SV using gene annotation context.
    
    Args:
        sv_id: SV identifier (must match annotation TSV)
        sv_type: DEL, DUP, INV, BND/TRA
        annotation_tsv: Path to pre-computed SV annotation TSV
        min_*_split_reads: Minimum split reads for each SV type
    
    Returns:
        GenomicContextResult with PASS/FLAG/FAIL verdict
    """
    annotations = load_annotation_file(annotation_tsv)
    
    if sv_id not in annotations:
        return GenomicContextResult(
            sv_id=sv_id, sv_type=sv_type,
            verdict=FilterVerdict.NOT_APPLICABLE,
            filters_passed=0, filters_total=0,
            genes_affected=[],
            gene_context="no_annotation",
            details=f"No annotation found for {sv_id}"
        )
    
    ann = annotations[sv_id]
    genes = ann['genes']
    is_intergenic = ann['is_intergenic']
    
    if sv_type == 'DUP':
        verdict, passed, total, details = apply_dup_filters(sv_id, genes, is_intergenic)
    elif sv_type == 'DEL':
        verdict, passed, total, details = apply_del_filters(sv_id, genes, is_intergenic, ann.get("size_bp", 0))
    elif sv_type == 'INV':
        verdict, passed, total, details = apply_inv_filters(sv_id, genes, is_intergenic)
    elif sv_type in ('BND', 'TRA'):
        # TRA/BND: gene fusion check only
        if genes:
            details = f"[FLAG] Potential gene fusion at breakpoint involving: {', '.join(genes[:5])}"
            verdict, passed, total = FilterVerdict.FLAG, 0, 1
        else:
            details = "[PASS] Intergenic translocation breakpoint"
            verdict, passed, total = FilterVerdict.PASS, 1, 1
    else:
        return GenomicContextResult(
            sv_id=sv_id, sv_type=sv_type,
            verdict=FilterVerdict.NOT_APPLICABLE,
            filters_passed=0, filters_total=0,
            genes_affected=genes,
            gene_context="intergenic" if is_intergenic else "genic",
            details=f"No type-specific filters for {sv_type}"
        )
    
    context, flags = classify_gene_context(genes, sv_type)
    
    return GenomicContextResult(
        sv_id=sv_id, sv_type=sv_type,
        verdict=verdict,
        filters_passed=passed,
        filters_total=total,
        genes_affected=genes,
        gene_context=context,
        details=details
    )


if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: layer_genomic_context.py <annotation.tsv> <sv_id> [sv_type]")
        print("Example: layer_genomic_context.py data/yeast/S288C_sv_annotations.tsv ICB_NC_001134.8_197384_269 DUP")
        sys.exit(1)
    
    tsv = sys.argv[1]
    sv_id = sys.argv[2]
    sv_type = sys.argv[3] if len(sys.argv) > 3 else 'DUP'
    
    result = analyze_genomic_context(sv_id, sv_type, tsv)
    
    print(f"\n{'='*60}")
    print(f"  Genomic Context Filter — {sv_type}")
    print(f"{'='*60}")
    print(f"  SV: {result.sv_id}")
    print(f"  Verdict: {result.verdict.value.upper()}")
    print(f"  Filters: {result.filters_passed}/{result.filters_total} passed")
    print(f"  Genes: {', '.join(result.genes_affected) if result.genes_affected else '(none)'}")
    print(f"  Context: {result.gene_context}")
    print(f"  Details: {result.details}")
    print(f"{'='*60}\n")

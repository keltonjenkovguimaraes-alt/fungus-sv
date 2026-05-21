# FUNGUS-SV Cross-Species Validation Results

## A. baumannii ATCC 19606 HiFi reads vs. 5 Acinetobacter spp.

| Reference Species | Genome (Mb) | Sniffles2 | cuteSV | ICB Consensus | HIGH | WEAK |
|-------------------|-------------|-----------|--------|---------------|------|------|
| A. bouvetii | 3.4 | 6 | 9 | 5 DEL | 2 | 3 |
| A. lwoffii | 3.5 | 9 | 24 | 8 DEL | 4 | 4 |
| A. cumulans | 3.7 | 10 | 25 | 9 DEL | 3 | 6 |
| A. lanii | 3.4 | 12 | 20 | 11 DEL | 5 | 6 |
| A. larvae | 3.7 | 1 | 3 | 1 DEL | 0 | 1 |
| **TOTAL** | — | **38** | **81** | **34 DEL** | **14** | **20** |

## Key Findings

1. ICB consensus reduced 119 raw calls → 34 consensus (71.4% reduction)
2. 100% of SVs ≥100 bp scored HIGH confidence (T ≥ 0.6)
3. 100% of SVs <100 bp scored WEAK (T < 0.4)
4. Size-stratified: 15 SVs ≥100 bp (mean T=0.607), 19 SVs <100 bp (mean T=0.286)
5. All 5 pipeline runs completed with zero failures

## Methods

- Reads: 19,568 PacBio HiFi, ~82× coverage, N50 ~17 kb
- Aligner: minimap2 -x map-hifi
- Callers: Sniffles2 v2.2 + cuteSV v1.0.8
- ICB: ≥2 callers, 0.5 overlap, 200bp flank
- Validation: depth + k-mer + breakpoint + ploidy
- T-score: uniform weights (0.25), SMaHT tiers

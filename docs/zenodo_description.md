# FUNGUS-SV: Structural Variant Discovery and Triangulation-Based Validation for Haploid Fungal Genomes

This dataset contains an interactive HTML report of structural variants (SVs) 
detected between the query strain CICC-1445 and five reference strains of 
*Saccharomyces cerevisiae* (S288C, BJ4, IMX2600, Makgeolli, SX2).

**Interactive report:** [fungus-sv.netlify.app](https://fungus-sv.netlify.app)

---

## Methods

- **Sequencing:** PacBio HiFi reads from CICC-1445 (SRR18210299, 274,915 reads, ~20 kb N50)
- **Alignment:** minimap2 (map-hifi preset)
- **SV Callers:** Sniffles2 + cuteSV + SVIM (ICB consensus, ≥2 caller agreement)
- **Validation:** 5-layer triangulation (local assembly, read depth, k-mer spectrum, breakpoint junction, ploidy confirmation)
- **Orthogonal validation:** Read depth analysis (samtools depth, DHFFC), split-read junction detection (samtools view), Samplot visualization
- **Reference genomes:** S288C (NCBI R64-3-1), BJ4, IMX2600, Makgeolli, SX2
- **Gene annotation:** NCBI RefSeq GFF (R64-3-1), 6,459 genes

---

## Results Summary

| Reference | Total SVs | DEL | DUP | INV |
|-----------|-----------|-----|-----|-----|
| S288C     | 277       | 248 | 18  | 11  |
| BJ4       | 165       | 140 | 13  | 12  |
| IMX2600   | 303       | 276 | 19  | 8   |
| Makgeolli | 250       | 225 | 17  | 8   |
| SX2       | 290       | 261 | 18  | 11  |

---

## Orthogonal Validation (CICC-1445 vs S288C)

### Genome-Wide Depth Validation (DHFFC)
| Category | Count | % | Threshold |
|----------|-------|---|-----------|
| Confirmed | 115 | 43.2% | DEL: DHFFC < 0.3, DUP: DHFFC > 2.0 |
| Weak | 67 | 25.2% | DEL: 0.3–0.7, DUP: 1.3–2.0 |
| Not DEL/DUP | 84 | 31.6% | Contradicts called SV type |

### Inversion Split-Read Validation
- **11/11 inversions (100%) confirmed** by split-read junction evidence
- Split reads per breakpoint: 32–1,091
- All scored CONTRADICTED (T=0.167) by triangulation despite confirmation — systematic under-scoring of inversions documented

### Key Validated Loci
| Locus | Type | Size | DHFFC | Result |
|-------|------|------|-------|--------|
| rDNA (chrXII) | DEL | 3,657 bp | 0.03 | Confirmed |
| chrII 62 kb | INV | 62,196 bp | N/A | Confirmed (split reads) |
| FLO1 | DEL | 554 bp | 1.24 | Not simple DEL |
| YBL005W-B (Ty2) | DEL | 5,921 bp | 0.30 | Confirmed |
| YBR012W-B (Ty2) | DEL | 5,922 bp | 0.23 | Confirmed |

---

## Calibrated Parameters (Haploid Fungi)

| Parameter | Original (Human Diploid) | Calibrated (Haploid Fungi) | Source |
|-----------|-------------------------|---------------------------|--------|
| DHFFC DEL | < 0.7 | **< 0.3** | Pedersen & Quinlan (2019); this study |
| DHFFC DUP | > 1.3 | **> 2.0** | Pedersen & Quinlan (2019); this study |
| Depth weight | 0.25 | **0.35** | Liu et al. (2024); this study |
| Breakpoint weight | 0.20 | **0.30** | Liu et al. (2024); this study |
| Assembly weight | 0.30 | **0.20** | This study |
| k-mer weight | 0.25 | **0.15** | This study |

---

## Contents

- `FUNGUS_SV_report.html` — Interactive genome track visualization with gene annotations and validation results
- `S288C_vs_CICC1445_ideogram.png` — High-resolution ideogram with affected gene labels
- `samplot_images/` — Samplot visualizations of key validated loci (rDNA DEL, chrII INV, FLO1)

---

## Live Report

**https://fungus-sv.netlify.app**

## Source Code

**https://github.com/keltonjenkovguimaraes-alt/fungus-sv**

---

## License

MIT License — Copyright (c) 2026 Kelton Guimarães et al.

# FUNGUS-SV Manual Curation Guide
## S288C vs CICC-1445 — CONTRADICTED Calls

### How to Curate Each Image

For each Samplot image, look at:
1. **Read depth (gray background)** — Does it drop/increase where expected?
2. **Split reads (colored lines)** — Are there reads spanning the breakpoint?
3. **Discordant pairs** — Do paired reads support the SV type?

### Scoring
- **YES** = Clear evidence, real SV
- **MAYBE** = Some evidence but uncertain
- **NO** = No supporting evidence, likely false positive

### Images Generated (17)

#### Inversions (11 images) — All previously confirmed by split reads
- ICB_NC_001134.8_197380_22 (62 kb, chrII)
- ICB_NC_001136.10_437149_68 (370 bp, chrIV)
- ICB_NC_001137.3_443710_84 (5.6 kb, chrV)
- ICB_NC_001139.9_535210_119 (6.4 kb, chrVII)
- ICB_NC_001142.9_206574_164 (168 kb, chrX)
- ICB_NC_001142.9_710189_165 (16.6 kb, chrX)
- ICB_NC_001145.3_362919_205 (70.5 kb, chrXIII)
- ICB_NC_001146.8_562033_222 (6 kb, chrXIV)
- ICB_NC_001147.6_594821_240 (6 kb, chrXV)
- ICB_NC_001224.1_51576_268 (16.5 kb, chrM)
- ICB_NC_001138.5_143849_271 (23.4 kb, chrVI)

#### Duplications (3 images) — Previously confirmed/likely by depth
- ICB_NC_001134.8_801640_270 (2.9 kb, DHFFC=2.85)
- ICB_NC_001144.5_451418_275 (17.5 kb, rDNA, DHFFC=82.3)
- ICB_NC_001133.9_205662_11 (913 bp, FLO1, DHFFC=1.50)

#### Deletions (3 images) — Confirmed by depth but CONTRADICTED by v2
- ICB_NC_001134.8_221032_19 (5.9 kb, YBL005W-B Ty2)
- ICB_NC_001134.8_259571_20 (5.9 kb, YBR012W-B Ty2)
- ICB_NC_001134.8_427856_21 (1.9 kb)


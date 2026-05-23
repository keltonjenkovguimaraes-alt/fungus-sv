# DETALHES DO PROJETO — FUNGUS-SV

## Descrição resumida

O projeto consiste no desenvolvimento de um pipeline computacional para detecção e validação de variantes estruturais (Structural Variants — SVs) em genomas haploides de fungos, utilizando dados de sequenciamento PacBio HiFi. O FUNGUS-SV combina múltiplos chamadores de SV por meio de um construtor de consenso (ICB) e valida cada SV candidato através de cinco camadas de evidência ortogonal (montagem local, profundidade de leitura, espectro de k-mer, junções de breakpoint e confirmação de ploidia), gerando um escore de confiança (T-score) para cada variante detectada.

---

## Introdução e justificativa

As variantes estruturais (SVs) são alterações genômicas maiores que 50 pares de bases, incluindo deleções, inserções, duplicações, inversões e translocações. Embora menos frequentes que os SNPs, as SVs contribuem com uma proporção maior de divergência genômica entre indivíduos e estão associadas a fenótipos clinicamente relevantes, como resistência a antifúngicos, virulência e adaptação ao hospedeiro em fungos patogênicos (Liu et al., 2024; Chen et al., 2023).

Em organismos não-modelo, como fungos dos gêneros *Sporothrix*, *Candida* e *Aspergillus*, não existem conjuntos de referência (benchmarks) para validação de SVs, ao contrário do genoma humano que dispõe de recursos como GIAB, HGSVC e SMaHT. Essa ausência de padrão-ouro dificulta a avaliação da acurácia dos chamadores de SV e limita estudos genômicos funcionais nesses organismos (Dunn et al., 2024; Todd et al., 2025).

O FUNGUS-SV foi desenvolvido especificamente para preencher essa lacuna. Utilizando uma abordagem de triangulação — onde múltiplas camadas de evidência independentes são combinadas para pontuar cada SV — o pipeline permite que pesquisadores priorizem SVs de alta confiança para validação experimental, reduzindo o investimento em PCRs de variantes falso-positivas.

O pipeline foi validado em dados reais de sequenciamento PacBio HiFi de *Acinetobacter baumannii* ATCC 19606 alinhados contra cinco espécies diferentes de *Acinetobacter* (*A. bouvetii*, *A. lwoffii*, *A. cumulans*, *A. lanii*, *A. larvae*), detectando 34 SVs consenso. A análise de ablação demonstrou que a combinação de duas ou mais camadas de evidência produz chamadas de ALTA confiança (T ≥ 0,6), enquanto SVs suportados por uma única camada recebem pontuação BAIXA (T < 0,4), validando a premissa central da triangulação.

---

## Objetivos

### Objetivo Geral

Desenvolver um pipeline computacional de código aberto para detecção e validação de variantes estruturais em genomas haploides sem conjunto de referência (benchmark), utilizando dados de sequenciamento PacBio HiFi e uma abordagem de triangulação de evidências ortogonais.

### Objetivos Específicos

1. Implementar um construtor de consenso (ICB) que integre os resultados de 4 chamadores de SV (Sniffles2, cuteSV, SVIM, pbsv) com parâmetros validados pela literatura.

2. Desenvolver cinco camadas de validação ortogonal independentes do alinhamento: montagem local (Flye), assinatura de profundidade de leitura, espectro de k-mer (Jellyfish), junções de breakpoint e confirmação de ploidia.

3. Calibrar cada parâmetro do pipeline com base em artigos científicos revisados por pares (15+ publicações), documentando a origem de cada valor-limite (threshold) no código e na documentação.

4. Validar o pipeline em dados reais de sequenciamento, comparando os resultados entre diferentes espécies para demonstrar a capacidade de detecção de SVs.

5. Produzir um relatório por SV (report card) com pontuação estratificada por tamanho, mapeamento para níveis de confiança SMaHT e recomendação para validação experimental.

6. Disponibilizar o pipeline como software de código aberto com 6 ambientes conda isolados para evitar conflitos de dependências.

---

## Problemas de pesquisa

A detecção de variantes estruturais em fungos haploides não-modelo enfrenta desafios significativos: (1) não existem conjuntos de referência (benchmarks) como o GIAB humano para medir precisão e recall; (2) os chamadores de SV individuais produzem muitos falsos positivos, especialmente em regiões repetitivas; (3) não há pesos empíricos publicados para combinar diferentes tipos de evidência de SV; (4) a maioria dos parâmetros dos chamadores foi otimizada para genomas humanos diploides, não para fungos haploides pequenos.

Pergunta central: **É possível validar computacionalmente variantes estruturais em genomas fúngicos haploides sem um conjunto de referência, utilizando triangulação de evidências ortogonais?**

Perguntas secundárias:
- A combinação de múltiplos chamadores de SV por consenso reduz falsos positivos em comparação com chamadores individuais?
- As camadas de evidência ortogonal (profundidade, k-mer, breakpoint, montagem local) são independentes ou correlacionadas?
- O T-score de triangulação correlaciona-se com o tamanho do SV e o número de camadas de suporte?
- Quais parâmetros precisam ser recalibrados para genomas fúngicos versus humanos?

---

## Método científico

### 1. Construção do pipeline e isolamento de ambientes

O FUNGUS-SV foi implementado em Python 3.11 com gerenciamento de workflows via Snakemake. Foram criados 6 ambientes conda isolados para evitar conflitos de dependências entre ferramentas: `sv_align` (minimap2, samtools), `sv_call` (Sniffles2, cuteSV, SVIM, pbsv), `sv_valid` (numpy, scipy, pandas, pysam), `sv_lar` (Flye, minimap2), `sv_kmers` (Jellyfish) e `sv_plot` (matplotlib, seaborn). Cada módulo do pipeline executa em seu ambiente específico.

### 2. Consenso de múltiplos chamadores (ICB)

Os reads HiFi são alinhados à referência com `minimap2 -x map-hifi`. Quatro chamadores de SV executam independentemente: Sniffles2 (v2.2), cuteSV (v1.0.8), SVIM (v1.4.2) e pbsv (v2.11.0). O ICB (Intersection-Consensus-Builder) agrupa SVs sobrepostos com ≥50% de sobreposição recíproca e 200 bp de flanco (parâmetros de Liu et al., 2024; Kronenberg et al., 2025). SVs suportados por ≥2 chamadores são retidos para validação.

### 3. Validação por camadas ortogonais

Cada SV consenso é avaliado por cinco camadas de evidência independentes:

| Camada | Método | Referência |
|--------|--------|-----------|
| Profundidade de leitura | Razão de cobertura região/flanco | Liu et al. (2024) |
| Espectro de k-mer | Presença/ausência de k-mers no banco Jellyfish | PAV (Ebert et al., 2021) |
| Junções de breakpoint | Reads split (SA tags) e soft-clipped; filtro MAPQ≥20 | SVvalidation (Zheng & Shang, 2024) |
| Montagem local (LAR) | Flye assembly dos reads na região do breakpoint | DeBreak (Chen et al., 2023) |
| Confirmação de ploidia | Taxa de heterozigosidade de SNVs via mpileup | Hammond et al. (2025) |

### 4. Cálculo do T-score

As pontuações das camadas são combinadas com pesos uniformes (0,25 cada) — pesos não calibrados, documentados como priors não-informativos. O T-score resultante (0-1) é mapeado para níveis de confiança: TRIPLE_TRIANGULATED (≥0,80), DOUBLE_CONFIRMED (≥0,60), SINGLE_LINE (≥0,40), WEAK (<0,40), seguindo a nomenclatura SMaHT (Zhang et al., 2025).

### 5. Validação experimental do pipeline

Validação experimental do pipeline
O pipeline foi validado com dados reais: 19.568 reads PacBio HiFi (N50 ~17 kb, cobertura ~82×) de *Acinetobacter baumannii* ATCC 19606 (DRR718942) alinhados contra genomas de referência de 12 linhagens clínicas de *A. baumannii* (AB30, MRSN15313, DETAB-E51, XH1056, UC23022, 6080, 280820, 966CSF, Aci4735, AR_0083, XH1037, SRM25) e 5 espécies adicionais de *Acinetobacter* (*A. bouvetii*, *A. lwoffii*, *A. cumulans*, *A. lanii*, *A. larvae*). 

Para cada uma das 17 comparações, o pipeline completo foi executado (alinhamento → ICB → validação → relatório). No total, foram detectados 860 SVs consenso nas comparações intraespecíficas (média de 72 por linhagem), com 85% (727/860) classificados como ALTA confiança (T ≥ 0,6). A validação por LAR (Local Assembly Refinement) confirmou 59% (23/39) dos SVs de topo testados e ofereceu suporte parcial para 31% adicionais, totalizando 90% de suporte. 

A análise de ablação demonstrou que a combinação de duas ou mais camadas de evidência produz chamadas de ALTA confiança, enquanto SVs suportados por uma única camada recebem pontuação BAIXA. A estratificação por tamanho revelou que 100% dos SVs ≥100 bp pontuam ALTO e 100% dos SVs <100 bp pontuam BAIXO, consistente com a expectativa de que SVs pequenos têm suporte de validação limitado.

---

## Resultados esperados

1. **Detecção de SVs em comparações entre espécies**: O pipeline deve detectar deleções, inserções, duplicações e inversões quando os reads provêm de um isolado diferente da referência. Nas validações iniciais com *Acinetobacter*, foram detectados 34 SVs consenso (todos deleções) em 5 comparações interespecíficas.

2. **Validação da triangulação**: A análise de ablação demonstrou que SVs suportados por 2 camadas de evidência recebem pontuação ALTA (T ≥ 0,6), enquanto SVs com 1 camada apenas recebem pontuação BAIXA (T < 0,4), confirmando que a combinação de evidências ortogonais melhora a discriminação.

3. **Estratificação por tamanho**: 100% dos SVs ≥100 bp foram classificados como ALTA confiança; 100% dos SVs <100 bp foram classificados como BAIXA confiança, consistente com a expectativa de que SVs pequenos têm suporte de validação limitado.

4. **Metodologia de benchmark interno**: Documento metodológico descrevendo como construir um conjunto de referência (truth set) baseado em montagem *de novo* para qualquer organismo não-modelo, permitindo validação de pipelines de SV sem depender de benchmarks externos.

5. **Disponibilização de software livre**: Pipeline completo disponível no GitHub com documentação honesta, parâmetros rastreáveis a artigos científicos e 6 ambientes conda isolados.

6. **Formação acadêmica**: O projeto proporciona formação em bioinformática, genômica computacional, desenvolvimento de pipelines, revisão sistemática de literatura e documentação científica, com potencial para publicação em periódico indexado e apresentação em congressos.

---

## Referências

1. Liu, Y.H. et al. (2024). Tradeoffs in alignment and assembly-based methods for structural variant detection with long-read sequencing data. *Nature Communications*, 15:2447.

2. Chen, Y. et al. (2023). Deciphering the exact breakpoints of structural variations using long sequencing reads with DeBreak. *Nature Communications*, 14:283.

3. Dunn, T. et al. (2024). Jointly benchmarking small and structural variant calls with vcfdist. *Genome Biology*, 25:253.

4. Todd, C. et al. (2025). SV-JIM: detailed pairwise structural variant calling using long-reads and genome assemblies. *Methods*, 234:305-313.

5. Kronenberg, Z. et al. (2025). The Platinum Pedigree: A long-read benchmark for genetic variants. *Nature Methods*.

6. Zhang, Y. et al. (2025). Comprehensive benchmarking of somatic structural variant detection at ultra-low allele fractions. *bioRxiv*, 2025.09.18.677206.

7. Zheng, Y. & Shang, X. (2024). SVvalidation: A long-read-based validation method for genomic structural variation. *PLOS ONE*, 19(1):e0291741.

8. Ebert, P. et al. (2021). Haplotype-resolved diverse human genomes and integrated analysis of structural variation. *Science*, 372:eabf7117.

9. Hammond, N. et al. (2025). Analytical validation of germline small variant detection using long-read HiFi genome sequencing. *Genome Research*, 35:1-9.

10. Nkouamedjo Fankep, R.C. et al. (2025). SV-MeCa: an XGBoost-based meta-caller approach for structural variant calling from short-read data. *BMC Bioinformatics*, 26:218.

---

*Pipeline disponível em: https://github.com/keltonjenkovguimaraes-alt/fungus-sv*

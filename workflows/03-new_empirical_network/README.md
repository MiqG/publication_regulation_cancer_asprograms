## New Empirical Network - Extend SF-exon networks determined empirically

Workflows extending splicing factor-exon networks obtained from experiments charaterizing exon inclusion changes upon perturbing SFs.

Outputs in `results/new_empirical_network/files` and `results/new_empirical_network/figures`.

## Outline
- `01-network_inference.smk`: create splicing factor exon networks from Rogalska2024 and add them to existing networks generated in AngladaGirotto2024
- `02-network_evaluation.smk`: evalutate predictive power of former and extended SF-exon networks.
- `03-cancer_splicing_programs.smk`: redefine cancer splicing programs (oncogenic-like and tumor suppressor-like) considering the additional splicing factors included in extended networks. This reanalyzed TCGA datasets across 14 cohorts.
- `04-carcinogenesis.smk`: reanalyze Danielsson *et al.* to check that extended cancer splicing programs show coordinated regulation along carcinogenesis experiment.


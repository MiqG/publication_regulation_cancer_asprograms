## SF programs in differentiation - activity analysis of cancer splicing programs during developmental tissue differentiation

Workflow analyzing Cardoso-Moreira *et al.* through our approach to assess the activity of splicing factors in cancer splicing programs during developmental differentiation, from embryo to adult, of human tissues. 

Outputs in `results/sf_programs_in_differentiation/files` and `results/sf_programs_in_differentiation/figures`.

## Outline
- `01-network_inference.smk`: create splicing factor exon networks from Rogalska2024 and add them to existing networks generated in AngladaGirotto2024
- `02-network_evaluation.smk`: evalutate predictive power of former and extended SF-exon networks.
- `03-cancer_splicing_programs.smk`: redefine cancer splicing programs (oncogenic-like and tumor suppressor-like) considering the additional splicing factors included in extended networks. This reanalyzed TCGA datasets across 14 cohorts.
- `04-carcinogenesis.smk`: reanalyze Danielsson *et al.* to check that extended cancer splicing programs show coordinated regulation along carcinogenesis experiment.
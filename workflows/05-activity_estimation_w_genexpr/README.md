# Activity estimation with gene expression - Adjusting splicing factor activities from gene expression signatures to resemble splicing factor activities from exon inclusion

Workflows making SF-gene networks to estimate SF activity from gene expression signatures and using shallow ANNs to adjust gene-based activies to resemble those obtained with SF-exon networks.

Outputs in `results/activity_estimation_w_genexpr/files` and `results/activity_estimation_w_genexpr/figures`.

## Outline
- `01-empirical_genexpr_networks.smk`: creates splicing factor-gene networks combining SF perturbation data from:
    - bulk transcriptomics: Rogalska2024 (analyzed in this study) and AngladaGirotto2024
    - single-cell transcriptomics: Perturb-seq
    
- `02-genexpr_networks_evaluation.smk`: performs leave-one-out network evaluation to benchmark the predictive power of exon-based and gene-based networks at predicting which splicing factor had been perturbed in held-out experiments (not used to create splicing factor networks for estimating splicing factor activities).

- `03-activity_modeling.smk`: computes exon inclusion and gene expression signatures of CCLE data to train different shallow ANNs to adjust splicing factor activities estimated from gene expression signatures to match those estimated from exon-inclusion signatures.

- `04-eval_carcinogenesis.smk`: estimate splicing factor activities using exon inclusion signatures, gene expression signatures and different adjustments in a bulk (Danielsson2019) and single-cell (Hodis2022) carcinogenesis dataset.

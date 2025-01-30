# Prepare supplementary tables and intermediate files for submission

Once all workflows have been run, this workflow (`main.smk`) generates the supplementary tables for the publication. Additionally, it places intermediary files from our analyses into a single folder.

Outputs in `results/prepare_submission/files`.

- Supplementary Tables:
    - Supplementary Table 1: redefined cancer splicing programs
    
- Intermediate files:
    ```shell
    prepare_submission/files/
    ├── 20250122-intermediate_files.zip
    ├── intermediate_files
    │   ├── datasets
    │   │   ├── activity
    │   │   │   ├── CardosoMoreira2019-EX.tsv.gz
    │   │   │   ├── Danielsson2013-bulkgenexpr-adjusted_fclayer.tsv.gz
    │   │   │   ├── Danielsson2013-EX.tsv.gz
    │   │   │   ├── Hodis2022-invitro_eng_melanoc-bulkgenexpr-adjusted_fclayer.tsv.gz
    │   │   │   ├── ReplogleWeissman2022_rpe1-bulkgenexpr-adjusted_fclayer.tsv.gz
    │   │   │   └── Urbanski2022-EX.tsv.gz
    │   │   ├── event_psi
    │   │   │   ├── Rogalska2024-ALTA.tsv.gz
    │   │   │   ├── Rogalska2024-ALTD.tsv.gz
    │   │   │   ├── Rogalska2024-EX.tsv.gz
    │   │   │   ├── Rogalska2024-INT.tsv.gz
    │   │   │   ├── Urbanski2022-ALTA.tsv.gz
    │   │   │   ├── Urbanski2022-ALTD.tsv.gz
    │   │   │   ├── Urbanski2022-EX.tsv.gz
    │   │   │   └── Urbanski2022-INT.tsv.gz
    │   │   ├── genexpr_tpm
    │   │   │   ├── Rogalska2024.tsv.gz
    │   │   │   └── Urbanski2022.tsv.gz
    │   │   ├── metadata
    │   │   │   ├── Hodis2022-invitro_eng_melanoc.tsv.gz
    │   │   │   ├── ReplogleWeissman2022_rpe1.tsv.gz
    │   │   │   ├── Rogalska2024.tsv.gz
    │   │   │   └── Urbanski2022.tsv.gz
    │   │   └── signatures
    │   │       ├── CardosoMoreira2019-EX.tsv.gz
    │   │       ├── Danielsson2013-EX.tsv.gz
    │   │       ├── Danielsson2013-genexpr_tpm.tsv.gz
    │   │       ├── Hodis2022-invitro_eng_melanoc-genexpr_cpm.tsv.gz
    │   │       ├── ReplogleWeissman2022_rpe1-genexpr_cpm.tsv.gz
    │   │       ├── Rogalska2024-EX.tsv.gz
    │   │       ├── Urbanski2022-EX.tsv.gz
    │   │       └── Urbanski2022-genexpr_tpm.tsv.gz
    │   ├── models
    │   │   ├── from_bulkgenexpr_to_EX
    │   │   │   ├── ewlayer
    │   │   │   │   ├── weights-0.pth
    │   │   │   │   ├── weights-1.pth
    │   │   │   │   ├── weights-2.pth
    │   │   │   │   ├── weights-3.pth
    │   │   │   │   └── weights-4.pth
    │   │   │   └── fclayer
    │   │   │       ├── weights-0.pth
    │   │   │       ├── weights-1.pth
    │   │   │       ├── weights-2.pth
    │   │   │       ├── weights-3.pth
    │   │   │       └── weights-4.pth
    │   │   ├── from_bulkscgenexpr_to_EX
    │   │   │   ├── ewlayer
    │   │   │   │   ├── weights-0.pth
    │   │   │   │   ├── weights-1.pth
    │   │   │   │   ├── weights-2.pth
    │   │   │   │   ├── weights-3.pth
    │   │   │   │   └── weights-4.pth
    │   │   │   └── fclayer
    │   │   │       ├── weights-0.pth
    │   │   │       ├── weights-1.pth
    │   │   │       ├── weights-2.pth
    │   │   │       ├── weights-3.pth
    │   │   │       └── weights-4.pth
    │   │   └── from_scgenexpr_to_EX
    │   │       ├── ewlayer
    │   │       │   ├── weights-0.pth
    │   │       │   ├── weights-1.pth
    │   │       │   ├── weights-2.pth
    │   │       │   ├── weights-3.pth
    │   │       │   └── weights-4.pth
    │   │       └── fclayer
    │   │           ├── weights-0.pth
    │   │           ├── weights-1.pth
    │   │           ├── weights-2.pth
    │   │           ├── weights-3.pth
    │   │           └── weights-4.pth
    │   └── networks
    │       ├── experimentally_derived_regulons_pruned-bulkgenexpr
    │       │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │       │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │       │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │       │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │       │   ├── ENCOREKD-benchmark-delta_psi.tsv.gz
    │       │   ├── ENCOREKO-benchmark-delta_psi.tsv.gz
    │       │   └── Rogalska2024-HELA_CERVIX-delta_psi.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-bulkscgenexpr
    │       │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │       │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │       │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │       │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │       │   ├── ENCOREKD-benchmark-delta_psi.tsv.gz
    │       │   ├── ENCOREKO-benchmark-delta_psi.tsv.gz
    │       │   ├── ReplogleWeissman2022-K562_essential-log2fc_genexpr.tsv.gz
    │       │   ├── ReplogleWeissman2022-K562_gwps-log2fc_genexpr.tsv.gz
    │       │   ├── ReplogleWeissman2022-rpe1-log2fc_genexpr.tsv.gz
    │       │   └── Rogalska2024-HELA_CERVIX-delta_psi.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-scgenexpr
    │       │   ├── ReplogleWeissman2022-K562_essential-log2fc_genexpr.tsv.gz
    │       │   ├── ReplogleWeissman2022-K562_gwps-log2fc_genexpr.tsv.gz
    │       │   └── ReplogleWeissman2022-rpe1-log2fc_genexpr.tsv.gz
    │       └── experimentally_derived_regulons_pruned_w_viper_networks-EX
    │           ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │           ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │           ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │           ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │           ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │           ├── ENCOREKD-K562-delta_psi.tsv.gz
    │           ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │           ├── ENCOREKO-K562-delta_psi.tsv.gz
    │           └── Rogalska2024-HELA_CERVIX-delta_psi.tsv.gz
    └── supplementary_tables
        └── suptab01_cancer_splicing_programs.tsv.gz

    23 directories, 89 files
    ```

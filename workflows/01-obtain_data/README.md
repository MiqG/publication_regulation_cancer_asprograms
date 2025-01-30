# Obtain Data - Workflows to download (and align) raw data

## Outline
1. `download_general.smk`
    - downloads data from other publications and databases
    - output directories:
        - `data/raw/HGNC`
        - `data/raw/Harmonizome`
        - `data/raw/MSigDB`
        - `data/raw/STRINGDB`
        - `data/raw/viper_splicing_intermediate_files`: contains preprocesssed TCGA, Danielsson, network inference, and Cardoso-Moreira datasets.
        
2. `download_sf_networks_AngladaGirotto2024.smk`
    - downloads and aligns RNA-seq reads with vast-tools producing a dataset-level exon inclusion table (PSI) and gene expression table (TPM)
    - outputs in `data/raw/viper_splicing_networks`
    
3. `download_and_align_Rogalska2024.smk`
    - downloads and aligns RNA-seq reads with vast-tools producing a dataset-level exon inclusion table (PSI) and gene expression table (TPM)
    - outputs in `data/raw/articles/Rogalska2024`
    
4. `download_perturbseq.smk`
    - downloads preprocessed AnnData count files from [Replogle2022](https://doi.org/10.1016/j.cell.2022.05.013).
    - outputs in `data/raw/scPerturb`
    
5. `download_and_align_Urbanski2022.smk`
    - downloads and aligns RNA-seq reads with vast-tools producing a dataset-level exon inclusion table (PSI) and gene expression table (TPM)
    - outputs in `data/raw/articles/Urbanski2022`
    
    
## Downloaded manually
- COSMIC
    Download manually from https://cancer.sanger.ac.uk/census. Expected output directory at:
    ```
    publication_spotter/data/raw/COSMIC/cancer_gene_census.tsv
    ```
    
- Melanoma tumorigenesis from [Hodis2022](https://doi.org/10.1126%2Fscience.abi8175):
    - downloaded by hand on 2024-05-14 due to login restrictions
    - outputs in `data/raw/articles/Hodis2022`
    - urls:
        - https://singlecell.broadinstitute.org/single_cell/data/public/SCP1334/engineered-melanocytes?filename=invitro_eng_melanoc_logTP10K.txt.gz
        - https://singlecell.broadinstitute.org/single_cell/data/public/SCP1334/engineered-melanocytes?filename=invivo_melanocytes_downsampled_logTP10K.txt.gz
        - https://singlecell.broadinstitute.org/single_cell/data/public/SCP1334/engineered-melanocytes?filename=invitro_invivo_all_metadatafile_mod_withCelltypes.csv

import os

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
DATA_DIR = os.path.join(ROOT,'data','raw','articles')
SUPPORT_DIR = os.path.join(ROOT,'support')

##### RULES #####
rule all:
    input:
        # Hodis2022 (by hand)
        # Becker2021
        os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_stromal_220213.rds"),
        os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_normal_epithelial_220213.rds"),
        os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_immune_220213.rds"),
        os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_all_epithelial_220213.rds"), # contains cancer cells as well
        os.path.join(DATA_DIR,"Becker2021","stromal_celltypes_rna.tsv"),
        os.path.join(DATA_DIR,"Becker2021","normal_celltypes_rna.tsv"), # corresponds to 'epithelial_celltypes_rna.tsv'
        os.path.join(DATA_DIR,"Becker2021","immune_celltypes_rna.tsv"),
        # Boiarsky2022
        os.path.join(DATA_DIR,"Boiarsky2022","GSE193531_umi-count-matrix.csv.gz"),
        os.path.join(DATA_DIR,"Boiarsky2022","GSE193531_cell-level-metadata.csv.gz")
        
        
rule download_Becker2021:
    params:
        stromal_seurat = "https://drive.google.com/file/d/135CVx9DZl7U3dIoXpLrx0p6AowBwdBtF/view?usp=drive_link",
        normal_seurat = "https://drive.google.com/file/d/12wzDPipillWBI-EmvGFqNSdu349N2WsB/view?usp=drive_link",
        immune_seurat = "https://drive.google.com/file/d/130Vb-OYXBC_oDSEAhl8s6I33S8rJS9I2/view?usp=drive_link",
        merged_seurat = "https://drive.google.com/file/d/12szgIYXLBwpLUXkK92A_TuLqZfQQAFAi/view?usp=drive_link",
        stromal_metadata = "https://drive.google.com/file/d/1HTeUdOJYAngj1QDOMnx1RwwWZXDojpvg/view?usp=drive_link",
        normal_metadata = "https://drive.google.com/file/d/17uo3Xs2BY8k4_psomuKIc7YDKiYcEnCD/view?usp=drive_link",
        immune_metadata = "https://drive.google.com/file/d/1hul6CFhrYGRf0fg3CNJ3BszQEv73aCW-/view?usp=drive_link",
    output:
        stromal_seurat = os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_stromal_220213.rds"),
        normal_seurat = os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_normal_epithelial_220213.rds"),
        immune_seurat = os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_immune_220213.rds"),
        merged_seurat = os.path.join(DATA_DIR,"Becker2021","Final_scHTAN_colon_all_epithelial_220213.rds"),
        stromal_metadata = os.path.join(DATA_DIR,"Becker2021","stromal_celltypes_rna.tsv"),
        normal_metadata = os.path.join(DATA_DIR,"Becker2021","normal_celltypes_rna.tsv"), # corresponds to 'epithelial_celltypes_rna.tsv'
        immune_metadata = os.path.join(DATA_DIR,"Becker2021","immune_celltypes_rna.tsv")
    shell:
        """
        set -eo pipefail
        
        gdown --user-agent="Chrome" --no-check-certificate --fuzzy {params.stromal_seurat} -O {output.stromal_seurat}
        gdown --user-agent="Chrome" --no-check-certificate --fuzzy {params.normal_seurat} -O {output.normal_seurat}
        gdown --user-agent="Chrome" --no-check-certificate --fuzzy {params.immune_seurat} -O {output.immune_seurat}
        gdown --user-agent="Chrome" --no-check-certificate --fuzzy {params.merged_seurat} -O {output.merged_seurat}
        
        gdown --user-agent="Chrome" --no-check-certificate --fuzzy {params.stromal_metadata} -O {output.stromal_metadata}
        gdown --user-agent="Chrome" --no-check-certificate --fuzzy {params.normal_metadata} -O {output.normal_metadata}
        gdown --user-agent="Chrome" --no-check-certificate --fuzzy {params.immune_metadata} -O {output.immune_metadata}
        
        echo "Done!"
        """
        
        
rule download_Boiarsky2022:
    params:
        genexpr = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193531&format=file&file=GSE193531%5Fumi%2Dcount%2Dmatrix%2Ecsv%2Egz",
        metadata = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193531&format=file&file=GSE193531%5Fcell%2Dlevel%2Dmetadata%2Ecsv%2Egz"
    output:
        genexpr = os.path.join(DATA_DIR,"Boiarsky2022","GSE193531_umi-count-matrix.csv.gz"),
        metadata = os.path.join(DATA_DIR,"Boiarsky2022","GSE193531_cell-level-metadata.csv.gz")
    shell:
        """
        set -eo pipefail
        
        wget --user-agent="Chrome" --no-check-certificate "{params.genexpr}" -O {output.genexpr}
        wget --user-agent="Chrome" --no-check-certificate "{params.metadata}" -O {output.metadata}
        
        echo "Done!"
        """

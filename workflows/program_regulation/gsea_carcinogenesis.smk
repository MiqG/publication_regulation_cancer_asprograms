import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
NETWORKS_DIR = os.path.join(ROOT,"results","network_inference")
RESULTS_DIR = os.path.join(ROOT,"results","program_regulation")

SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### RULES #####
rule all:
    input:
        os.path.join(RESULTS_DIR,"figures","gsea_carcinogenesis")
        
        
rule figures_gsea_carcinogenesis:
    input:
        carcinogenesis_bulk_activity = os.path.join(NETWORKS_DIR,"figures","eval_tumorigenesis","figdata","eval_tumorigenesis","protein_activity.tsv.gz"),
        carcinogenesis_bulk_hallmarks = os.path.join(NETWORKS_DIR,"files","gsea","tumorigenesis-genexpr-hallmarks.tsv.gz"),
        carcinogenesis_singlecell_activity = os.path.join(NETWORKS_DIR,"figures","eval_tumorigenesis_singlecell-Hodis2022-invitro_eng_melanoc","figdata","eval_tumorigenesis_singlecell","protein_activity.tsv.gz"),
        carcinogenesis_singlecell_hallmarks = os.path.join(NETWORKS_DIR,"files","gsea","Hodis2022-invitro_eng_melanoc-hallmarks.tsv.gz"),
        pertseq_rpe1_activity = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_rpe1-genexpr.tsv.gz"),
        pertseq_rpe1_hallmarks = os.path.join(RESULTS_DIR,"files","gsea","ReplogleWeissman2022_rpe1-hallmarks.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","gsea_carcinogenesis"))
    shell:
        """
        Rscript scripts/figures_gsea_carcinogenesis.R \
                    --carcinogenesis_bulk_activity_file={input.carcinogenesis_bulk_activity} \
                    --carcinogenesis_bulk_hallmarks_file={input.carcinogenesis_bulk_hallmarks} \
                    --carcinogenesis_singlecell_activity_file={input.carcinogenesis_singlecell_activity} \
                    --carcinogenesis_singlecell_hallmarks_file={input.carcinogenesis_singlecell_hallmarks} \
                    --pertseq_rpe1_activity_file={input.pertseq_rpe1_activity} \
                    --pertseq_rpe1_hallmarks_file={input.pertseq_rpe1_hallmarks} \
                    --figs_dir={output}
        """
        

    
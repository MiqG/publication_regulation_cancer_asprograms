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
VIPER_SPLICING_DIR = os.path.join(ROOT,"../../repositories/viper_splicing")

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}
ONTOLOGIES = ["hallmarks","reactome"]

##### RULES #####
rule all:
    input:
        # protein activity Urbanski2022
        expand(os.path.join(RESULTS_DIR,"files","gsea","Urbanski2022-{ontology_oi}.tsv.gz"), ontology_oi=ONTOLOGIES),

        # protein activity Urbanski2022
        os.path.join(RESULTS_DIR,"files","protein_activity","Urbanski2022-EX.tsv.gz"),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","gsea_carcinogenesis")
        

rule run_gsea:
    input:
        signature = os.path.join(PREP_DIR,"signatures","Urbanski2022-genexpr_tpm.tsv.gz"),
        msigdb_dir = os.path.join(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs"),
        gene_info = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","gsea","Urbanski2022-{ontology_oi}.tsv.gz")
    params:
        ontology_oi = "{ontology_oi}",
        random_seed = 1234,
        script_dir = SRC_DIR
    threads: 10
    shell:
        """
        Rscript {params.script_dir}/gsea_on_matrix.R \
                    --msigdb_dir={input.msigdb_dir} \
                    --signature_file={input.signature} \
                    --gene_info_file={input.gene_info} \
                    --ontology_oi={params.ontology_oi} \
                    --n_jobs={threads} \
                    --output_file={output}
        """
        
rule compute_protein_activity:
    input:
        signature = os.path.join(PREP_DIR,"signatures","Urbanski2022-EX.tsv.gz"),
        regulons_path = os.path.join(VIPER_SPLICING_DIR,"data","empirical_sf_networks-EX")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","Urbanski2022-EX.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
        
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
        

    
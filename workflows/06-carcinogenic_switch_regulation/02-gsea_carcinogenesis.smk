import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
CARCINOGENESIS_BULK_DIR = os.path.join(ROOT,"results","new_empirical_network")
CARCINOGENESIS_SC_DIR = os.path.join(ROOT,"results","activity_estimation_w_genexpr")
RESULTS_DIR = os.path.join(ROOT,"results","carcinogenic_switch_regulation")
VIPER_SPLICING_DIR = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","datasets")

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}
ONTOLOGIES = ["hallmarks","hallmarks_nomyc"]

PERT_GENEXPR_FILES = {
    "ReplogleWeissman2022_rpe1": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    "Urbanski2022": os.path.join(PREP_DIR,"signatures","Urbanski2022-genexpr_tpm.tsv.gz"),
    "Hodis2022-invitro_eng_melanoc": os.path.join(CARCINOGENESIS_SC_DIR,"files","signatures","Hodis2022-invitro_eng_melanoc-genexpr.tsv.gz"),
    "tumorigenesis-genexpr": os.path.join(CARCINOGENESIS_BULK_DIR,"files","signatures","tumorigenesis-genexpr.tsv.gz"),
    "CardosoMoreira2020": os.path.join(ROOT,"results","sf_programs_in_differentiation","files","signatures","CardosoMoreira2020-genexpr.tsv.gz")
}

REGULON_DIR = os.path.join(ROOT,"results","activity_estimation_w_genexpr")
REGULON_DIRS = {
    "EX": os.path.join(ROOT,"results","new_empirical_network","files","experimentally_derived_regulons_pruned_w_viper_networks-EX")
}

##### RULES #####
rule all:
    input:
        # GSEA
        expand(os.path.join(RESULTS_DIR,"files","gsea","{dataset}-{ontology_oi}.tsv.gz"), dataset=PERT_GENEXPR_FILES.keys(), ontology_oi=ONTOLOGIES),

        # protein activity Urbanski2022
        os.path.join(RESULTS_DIR,"files","protein_activity","Urbanski2022-EX.tsv.gz"),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","gsea_carcinogenesis")
        
        
rule run_gsea:
    input:
        signature = lambda wildcards: PERT_GENEXPR_FILES[wildcards.dataset],
        msigdb_dir = os.path.join(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs"),
        gene_info = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","gsea","{dataset}-{ontology_oi}.tsv.gz")
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
        regulons_path = REGULON_DIRS["EX"]
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
        carcinogenesis_bulk_genexpr = os.path.join(VIPER_SPLICING_DIR,'genexpr_tpm','tumorigenesis.tsv.gz'),
        carcinogenesis_bulk_activity = os.path.join(CARCINOGENESIS_BULK_DIR,"figures","carcinogenesis","figdata","carcinogenesis","protein_activity.tsv.gz"),
        carcinogenesis_bulk_hallmarks = os.path.join(RESULTS_DIR,"files","gsea","tumorigenesis-genexpr-hallmarks.tsv.gz"),
        carcinogenesis_bulk_metadata = os.path.join(VIPER_SPLICING_DIR,"metadata","tumorigenesis.tsv.gz"),
        
        carcinogenesis_singlecell_genexpr = os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-pseudobulk.tsv.gz"),
        carcinogenesis_singlecell_activity = os.path.join(CARCINOGENESIS_SC_DIR,"figures","eval_carcinogenesis","figdata","eval_carcinogenesis","protein_activity_singlecell.tsv.gz"),
        carcinogenesis_singlecell_hallmarks = os.path.join(RESULTS_DIR,"files","gsea","Hodis2022-invitro_eng_melanoc-hallmarks.tsv.gz"),
        carcinogenesis_singlecell_metadata = os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-conditions.tsv.gz"),
        
        pertseq_activity = os.path.join(RESULTS_DIR,"figures","upstream_regulators","figdata","upstream_regulators","cancer_program_activity.tsv.gz"),
        pertseq_hallmarks = os.path.join(RESULTS_DIR,"files","gsea","ReplogleWeissman2022_rpe1-hallmarks.tsv.gz"),
        pertseq_genexpr = os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
        
        urbanski_metadata = os.path.join(PREP_DIR,"metadata","Urbanski2022.tsv.gz"),
        urbanski_genexpr = os.path.join(PREP_DIR,'genexpr_tpm',"Urbanski2022.tsv.gz"),
        urbanski_ex = os.path.join(PREP_DIR,'event_psi',"Urbanski2022-EX.tsv.gz"),
        urbanski_activity = os.path.join(RESULTS_DIR,"files","protein_activity","Urbanski2022-EX.tsv.gz"),
        urbanski_hallmarks = os.path.join(RESULTS_DIR,"files","gsea","Urbanski2022-hallmarks.tsv.gz"),
        
        msigdb_dir = os.path.join(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs"),
        chea = os.path.join(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz"),
        splicing_factors = os.path.join(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv"),
        cancer_program = os.path.join(CARCINOGENESIS_BULK_DIR,'files','PANCAN','cancer_program.tsv.gz'),
        gene_annot = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","gsea_carcinogenesis"))
    shell:
        """
        Rscript scripts/figures_gsea_carcinogenesis.R \
                    --carcinogenesis_bulk_genexpr_file={input.carcinogenesis_bulk_genexpr} \
                    --carcinogenesis_bulk_activity_file={input.carcinogenesis_bulk_activity} \
                    --carcinogenesis_bulk_hallmarks_file={input.carcinogenesis_bulk_hallmarks} \
                    --carcinogenesis_bulk_metadata_file={input.carcinogenesis_bulk_metadata} \
                    --carcinogenesis_singlecell_genexpr_file={input.carcinogenesis_singlecell_genexpr} \
                    --carcinogenesis_singlecell_activity_file={input.carcinogenesis_singlecell_activity} \
                    --carcinogenesis_singlecell_hallmarks_file={input.carcinogenesis_singlecell_hallmarks} \
                    --carcinogenesis_singlecell_metadata_file={input.carcinogenesis_singlecell_metadata} \
                    --pertseq_activity_file={input.pertseq_activity} \
                    --pertseq_hallmarks_file={input.pertseq_hallmarks} \
                    --pertseq_genexpr_file={input.pertseq_genexpr} \
                    --urbanski_metadata_file={input.urbanski_metadata} \
                    --urbanski_genexpr_file={input.urbanski_genexpr} \
                    --urbanski_ex_file={input.urbanski_ex} \
                    --urbanski_activity_file={input.urbanski_activity} \
                    --urbanski_hallmarks_file={input.urbanski_hallmarks} \
                    --msigdb_dir={input.msigdb_dir} \
                    --chea_file={input.chea} \
                    --splicing_factors_file={input.splicing_factors} \
                    --cancer_program_file={input.cancer_program} \
                    --gene_annot_file={input.gene_annot} \
                    --figs_dir={output}
        """
    
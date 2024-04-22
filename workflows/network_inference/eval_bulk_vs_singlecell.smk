import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","network_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

PERT_GENEXPR_FILES = {
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ReplogleWeissman2022_K562_essential": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-log2_fold_change_cpm.tsv.gz"),
    "ReplogleWeissman2022_K562_gwps": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_gwps-log2_fold_change_cpm.tsv.gz")
}


##### RULES #####
rule all:
    input:
        # estimate splicing factor activity for K562 in bulk and single cell
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-genexpr.tsv.gz"), dataset=PERT_GENEXPR_FILES.keys()),
        
        # make figures
        os.path.join(RESULTS_DIR,"figures","eval_bulk_vs_singlecell"),
        
        
rule compute_protein_activity:
    input:
        signature = lambda wildcards: PERT_GENEXPR_FILES[wildcards.dataset],
        regulons_path = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-genexpr")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-genexpr.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
        
rule figures_regulon_evaluation:
    input:
        protein_activity_bulk = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-genexpr.tsv.gz"),
        protein_activity_singlecell = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_K562_essential-genexpr.tsv.gz"),
        cancer_program = os.path.join(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","eval_bulk_vs_singlecell"))
    shell:
        """
        Rscript scripts/figures_eval_bulk_vs_singlecell.R \
                    --protein_activity_bulk_file={input.protein_activity_bulk} \
                    --protein_activity_singlecell_file={input.protein_activity_singlecell} \
                    --cancer_program_file={input.cancer_program} \
                    --figs_dir={output}
        """

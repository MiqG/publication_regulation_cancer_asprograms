import os
import copy

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","activity_estimation_w_genexpr")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

PERT_SPLICING_FILES = {
    "ENASFS": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_psi_ena.tsv.gz'),
    "ENCOREKD_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_psi_encorekd_hepg2.tsv.gz'),
    "ENCOREKD_K562": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_psi_encorekd_k562.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_psi_encoreko_hepg2.tsv.gz'),
    "ENCOREKO_K562": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_psi_encoreko_k562.tsv.gz'),
    "Rogalska2024": os.path.join(PREP_DIR,'delta_psi','Rogalska2024-EX.tsv.gz')
}

PERT_GENEXPR_FILES = {
    "ENASFS": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_tpm_ena.tsv.gz'),
    "ENCOREKD_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_tpm_encorekd_hepg2.tsv.gz'),
    "ENCOREKD_K562": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_tpm_encorekd_k562.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_tpm_encoreko_hepg2.tsv.gz'),
    "ENCOREKO_K562": os.path.join(RAW_DIR,'viper_splicing_intermediate_files','benchmark','signatures_tpm_encoreko_k562.tsv.gz'),
    "Rogalska2024": os.path.join(PREP_DIR,'log2_fold_change_tpm','Rogalska2024-genexpr_tpm.tsv.gz'),
    "ReplogleWeissman2022_K562_essential-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    "ReplogleWeissman2022_rpe1-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    "ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz")
}

PERT_EVAL_FILES = {
    "EX": PERT_SPLICING_FILES,
    "genexpr": PERT_GENEXPR_FILES
}

OMIC_PERT_DICT = {
    "EX": "EX",
    "genexpr": "genexpr",
    "scgenexpr": "genexpr"
}

REGULON_SETS = {
    "EX": [
        "experimentally_derived_regulons_pruned_w_viper_networks-EX",
        "viper_networks-EX"
    ],
    "genexpr": [
        "experimentally_derived_regulons_pruned-bulkgenexpr",
        "experimentally_derived_regulons_pruned-scgenexpr",
        "experimentally_derived_regulons_pruned-bulkscgenexpr"
    ]
}

METHODS_ACTIVITY = ["viper"]

EVENT_TYPES = ["EX"]
OMIC_TYPES = ["genexpr","scgenexpr"] + EVENT_TYPES

OMIC_TYPE_LIST = []
METHOD_ACTIVITY_LIST = []
REGULON_SETS_LIST = []
DATASETS_LIST = []
for o in OMIC_TYPES:
    for d in PERT_EVAL_FILES[OMIC_PERT_DICT[o]]:
        for r in REGULON_SETS[OMIC_PERT_DICT[o]]:
            for m in METHODS_ACTIVITY:
                OMIC_TYPE_LIST.append(o)
                METHOD_ACTIVITY_LIST.append(m)
                REGULON_SETS_LIST.append(r)
                DATASETS_LIST.append(d)
                
##### RULES #####
rule all:
    input:
        # evaluate regulons
        ## run
        expand(
            os.path.join(RESULTS_DIR,"files","network_evaluation_scores","{method_activity}","{regulon_set}__{dataset}__{omic_type}.tsv.gz"), 
            zip, 
            method_activity=METHOD_ACTIVITY_LIST, regulon_set=REGULON_SETS_LIST, dataset=DATASETS_LIST, omic_type=OMIC_TYPE_LIST
        ),
        ## merge
        os.path.join(RESULTS_DIR,"files","network_evaluation_scores","merged.tsv.gz"),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","network_evaluation")
        

rule evaluate_regulons:
    input:
        signature = lambda wildcards: PERT_EVAL_FILES[OMIC_PERT_DICT[wildcards.omic_type]][wildcards.dataset],
        regulons = os.path.join(RESULTS_DIR,"files","{regulon_set}"),
        eval_labels = os.path.join(ROOT,"results","new_empirical_network","files","network_evaluation_labels")
    output:
        os.path.join(RESULTS_DIR,"files","network_evaluation_scores","{method_activity}","{regulon_set}__{dataset}__{omic_type}.tsv.gz")
    params:
        eval_labels = os.path.join(ROOT,"results","new_empirical_network","files","network_evaluation_labels","{dataset}.tsv.gz"),
        script_dir = SRC_DIR,
        shadow = "no",
        n_tails = "two",
        method_activity = "{method_activity}"
    threads: 1
    resources:
        # runtime = 3600*6, # h in seconds
        runtime = 60*24, # h in minutes 
        memory = 300, # GB
    shell:
        """
        nice Rscript {params.script_dir}/evaluate_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons} \
                    --eval_labels_file={params.eval_labels} \
                    --output_file={output} \
                    --method_activity={params.method_activity} \
                    --shadow_correction={params.shadow} \
                    --n_tails={params.n_tails}
        """
        
rule combine_evaluations:
    input:
        evaluations = [
            os.path.join(RESULTS_DIR,"files","network_evaluation_scores","{method_activity}","{regulon_set}__{dataset}__{omic_type}.tsv.gz").format(method_activity=m, regulon_set=r, dataset=d, omic_type=o) 
            for (m, r, d, o) in zip(METHOD_ACTIVITY_LIST, REGULON_SETS_LIST, DATASETS_LIST, OMIC_TYPE_LIST)
        ]
    output:
        os.path.join(RESULTS_DIR,"files","network_evaluation_scores","merged.tsv.gz")
    run:
        import os
        import pandas as pd
        
        print("Reading...")
        cols_oi = [
            "PERT_ID","auc_roc","auc_pr","avg_pr","mean_rank_percentile","median_rank_percentile",
            "grouping_var","n_pos_class","n_total","curves_type","eval_direction","eval_type","regulator",
            "n_networks_per_regulator","regulon_set_id","signature_id","shadow_correction","n_tails",
            "method_activity"
        ]
        evaluation = pd.concat([
            pd.read_table(f, usecols=cols_oi, low_memory=False).drop_duplicates().assign(
                omic_type = os.path.basename(f).replace(".tsv.gz","").split("__")[-1]
            ) 
            for f in input.evaluations
        ])
        
        print("Saving...")
        evaluation.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
      
    
rule figures_network_evaluation:
    input:
        evaluation = os.path.join(RESULTS_DIR,"files","network_evaluation_scores","merged.tsv.gz"),
        splicing_factors = os.path.join(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","network_evaluation"))
    shell:
        """
        Rscript scripts/figures_network_evaluation.R \
                    --evaluation_file={input.evaluation} \
                    --splicing_factors_file={input.splicing_factors} \
                    --figs_dir={output}
        """
       
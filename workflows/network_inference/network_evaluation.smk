import os
import copy

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","network_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

VIPER_SPLICING_DIR = os.path.join(ROOT,"../../repositories/viper_splicing")

PERT_SPLICING_FILES = {
    "ENCOREKD_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_psi_encorekd_hepg2.tsv.gz'),
    "ENCOREKD_K562": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_psi_encorekd_k562.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_psi_encoreko_hepg2.tsv.gz'),
    "ENCOREKO_K562": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_psi_encoreko_k562.tsv.gz'),
    "ENASFS": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_tpm_ena.tsv.gz'),
    "Rogalska2024": os.path.join(PREP_DIR,'delta_psi','Rogalska2024-EX.tsv.gz')
}

PERT_GENEXPR_FILES = {
    "ENCOREKD_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_tpm_encorekd_hepg2.tsv.gz'),
    "ENCOREKD_K562": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_tpm_encorekd_k562.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_tpm_encoreko_hepg2.tsv.gz'),
    "ENCOREKO_K562": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_tpm_encoreko_k562.tsv.gz'),
    "ENASFS": os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark','signatures_tpm_ena.tsv.gz'),
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

METADATA_FILES = [
    os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark',"metadata_encoreko.tsv.gz"),
    os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark',"metadata_encorekd.tsv.gz"),
    os.path.join(RAW_DIR,'viper_splicing_intermediary_files','benchmark',"metadata_ena.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","Rogalska2024.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-conditions.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-conditions.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches-conditions.tsv.gz")
]

REGULON_SETS = {
    "EX": ["experimentally_derived_regulons_pruned_w_viper_networks-EX"],
    "genexpr": [
        "experimentally_derived_regulons_pruned-genexpr",
        "experimentally_derived_regulons_pruned-scgenexpr",
    ]
}

METHODS_ACTIVITY = ["viper","correlation_pearson","correlation_spearman"]#,"gsea"]

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
                
# OMIC_TYPE_LIST = OMIC_TYPE_LIST[:4]
# METHOD_ACTIVITY_LIST = METHOD_ACTIVITY_LIST[:4]
# REGULON_SETS_LIST = REGULON_SETS_LIST[:4]
# DATASETS_LIST = DATASETS_LIST[:4]

##### RULES #####
rule all:
    input:
        # prepare evaluation labels
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"),
        
        # evaluate regulons
        ## run
        expand(
            os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{method_activity}","{regulon_set}__{dataset}__{omic_type}.tsv.gz"), 
            zip, 
            method_activity=METHOD_ACTIVITY_LIST, regulon_set=REGULON_SETS_LIST, dataset=DATASETS_LIST, omic_type=OMIC_TYPE_LIST
        ),
        ## merge
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged.tsv.gz"),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","network_evaluation")
        
        
rule make_evaluation_labels:
    input:
        metadatas = METADATA_FILES
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"))
    run:
        import pandas as pd
        
        perts_oi = ["KNOCKDOWN","KNOCKOUT","OVEREXPRESSION"]
        
        for f in input.metadatas:
            print(f)
            
            # load
            metadata = pd.read_table(f)
            
            os.makedirs(output.output_dir, exist_ok=True)
            if "encore" in f:
                for cell_line in metadata["cell_line"].unique():
                    # make labels
                    labels = metadata.loc[
                        metadata["cell_line"]==cell_line, 
                        ["cell_line","PERT_ENSEMBL","PERT_TYPE"]
                    ].drop_duplicates().copy()
                    
                    dataset = "ENCOREKD" if "kd" in os.path.basename(f) else "ENCOREKO"
                    labels["study_accession"] = dataset
                    
                    labels["PERT_ID"] = labels[
                        ["study_accession","cell_line","PERT_ENSEMBL","PERT_TYPE"]
                    ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                    
                    # only simple perturbations
                    labels = labels.loc[labels["PERT_TYPE"].isin(perts_oi)]
                    
                    # save
                    labels.dropna().to_csv(os.path.join(output.output_dir,"%s_%s.tsv.gz") % (dataset, cell_line), **SAVE_PARAMS)
            
            elif "Rogalska2024" in f:
                    cell_line = "HELA_CERVIX"
                    dataset = "Rogalska2024"
                    
                    # make labels
                    labels = metadata[["cell_line","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates().copy()
                    
                    labels["study_accession"] = dataset
                    
                    labels["PERT_ID"] = labels[
                        ["study_accession","cell_line","PERT_ENSEMBL","PERT_TYPE"]
                    ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                    
                    # only simple perturbations
                    labels = labels.loc[labels["PERT_TYPE"].isin(perts_oi)]
                    
                    # save
                    labels.dropna().to_csv(os.path.join(output.output_dir,"%s_%s.tsv.gz") % (dataset, cell_line), **SAVE_PARAMS)
            
            elif "singlecell" in f:
                metadata["study_accession"] = "ReplogleWeissman2022_"+dataset
                metadata["cell_line"] = "K562" if "K562" in dataset else "RPE1"
                metadata["PERT_TYPE"] = "KNOCKDOWN"
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line","PERT_ENSEMBL","PERT_TYPE"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)

                # "PERT_GENE" columns would be nice to have
                labels = metadata[["PERT_ID","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates()
                
                # only simple perturbations
                labels = labels.loc[labels["PERT_TYPE"].isin(perts_oi)]

                # save
                dataset_file = os.path.basename(f).replace("-conditions","")
                labels.dropna().to_csv(os.path.join(output.output_dir,dataset_file), **SAVE_PARAMS)
                
            elif "ENASFS" in f:
                # prepare labels
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line_name","PERT_ENSEMBL","PERT_TYPE"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                labels = metadata[["PERT_ID","PERT_GENE","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates()

                # only simple perturbations
                labels = labels.loc[labels["PERT_TYPE"].isin(perts_oi)]

                # save
                labels.dropna().to_csv(os.path.join(output.output_dir,"ENASFS.tsv.gz"), **SAVE_PARAMS)
                
                
        print("Done!")
        

rule evaluate_regulons:
    input:
        signature = lambda wildcards: PERT_EVAL_FILES[OMIC_PERT_DICT[wildcards.omic_type]][wildcards.dataset],
        regulons = os.path.join(RESULTS_DIR,"files","{regulon_set}"),
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels")
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{method_activity}","{regulon_set}__{dataset}__{omic_type}.tsv.gz")
    params:
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels","{dataset}.tsv.gz"),
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
            os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{method_activity}","{regulon_set}__{dataset}__{omic_type}.tsv.gz").format(method_activity=m, regulon_set=r, dataset=d, omic_type=o) 
            for (m, r, d, o) in zip(METHOD_ACTIVITY_LIST, REGULON_SETS_LIST, DATASETS_LIST, OMIC_TYPE_LIST)
        ]
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged.tsv.gz")
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
        evaluation_ex = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz"),
        evaluation_genexpr = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-genexpr.tsv.gz"),
        evaluation_scgenexpr = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-scgenexpr.tsv.gz"),
        protein_activity_ex = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_K562-EX.tsv.gz"),
        protein_activity_genexpr = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_K562-genexpr.tsv.gz"),
        protein_activity_scgenexpr = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-scgenexpr.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","eval_genexpr_vs_splicing"))
    shell:
        """
        Rscript scripts/figures_eval_genexpr_vs_splicing.R \
                    --evaluation_ex_file={input.evaluation_ex} \
                    --evaluation_genexpr_file={input.evaluation_genexpr} \
                    --protein_activity_ex_file={input.protein_activity_ex} \
                    --protein_activity_genexpr_file={input.protein_activity_genexpr} \
                    --figs_dir={output}
        """
       
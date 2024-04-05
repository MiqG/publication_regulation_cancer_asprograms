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

VIPER_SPLICING_DIR = os.path.join(ROOT,"../../repositories/viper_splicing")

EVENT_TYPES = ["EX"]
OMIC_TYPES = ["genexpr"] + EVENT_TYPES

PERT_SPLICING_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{omic_type}.tsv.gz')
}

PERT_GENEXPR_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz')
}

EVAL_DATASETS = list(PERT_GENEXPR_FILES.keys())

PERT_FILES = {
    "EX": PERT_SPLICING_FILES,
    "genexpr": PERT_GENEXPR_FILES
}

METADATA_FILES = [
    os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz")
]

REGULON_DIRS = {
    "genexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-genexpr"),
    "EX": os.path.join(VIPER_SPLICING_DIR,"data","empirical_sf_networks-EX")
}

SHADOWS = ["no"] # bug in viper does not allow shadow correction
N_TAILS = ["two"]

DATASETS = {
    "CCLE": {
        "genexpr": os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz"), 
        "EX": os.path.join(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
    }
}

##### RULES #####
rule all:
    input:
        # prepare evaluation labels
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"),
        
        # evaluate regulons
        ## run
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz"), dataset=EVAL_DATASETS, omic_type=OMIC_TYPES, shadow=SHADOWS, n_tails=N_TAILS),
        ## merge
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # estimate splicing factor activities in CCLE
        ## compute signatures within
        # expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"), zip, dataset=DATASETS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"), dataset=EVAL_DATASETS, omic_type=OMIC_TYPES),
        
        # make figures
        os.path.join(RESULTS_DIR,"figures","eval_genexpr_vs_splicing")
        
        
rule make_evaluation_labels:
    input:
        metadatas = METADATA_FILES
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"))
    run:
        import pandas as pd
        
        for f in input.metadatas:
            # load
            metadata = pd.read_table(f)
            
            os.makedirs(output.output_dir, exist_ok=True)
            if "ENCORE" in f:
                for cell_line in metadata["cell_line"].unique():
                    # make labels
                    labels = metadata.loc[
                        metadata["cell_line"]==cell_line, ["PERT_GENE","PERT_ENSEMBL"]
                    ].drop_duplicates().copy()
                    labels["PERT_ID"] = metadata["PERT_ENSEMBL"]
                    labels["PERT_TYPE"] = "KNOCKDOWN" if "KD" in os.path.basename(f) else "KNOCKOUT"
                    
                    dataset = "ENCOREKD" if "KD" in os.path.basename(f) else "ENCOREKO"
                    
                    # save
                    labels.dropna().to_csv(os.path.join(output.output_dir,"%s_%s.tsv.gz") % (dataset, cell_line), **SAVE_PARAMS)
                
            elif "ENASFS" in f:
                # prepare labels
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line_name","PERT_ENSEMBL"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                labels = metadata[["PERT_ID","PERT_GENE","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates()

                # save
                labels.dropna().to_csv(os.path.join(output.output_dir,"ENASFS.tsv.gz"), **SAVE_PARAMS)
                
                
        print("Done!")
        

rule evaluate_regulons:
    input:
        signature = lambda wildcards: PERT_FILES[wildcards.omic_type][wildcards.dataset],
        regulons = lambda wildcards: REGULON_DIRS[wildcards.omic_type],
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels","{dataset}.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz")
    params:
        script_dir = SRC_DIR,
        shadow = "{shadow}",
        n_tails = "{n_tails}"
    shell:
        """
        nice Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons} \
                    --eval_labels_file={input.eval_labels} \
                    --output_file={output} \
                    --shadow_correction={params.shadow} \
                    --n_tails={params.n_tails}
        """
        
        
rule combine_evaluations:
    input:
        evaluations = [os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz").format(dataset=d, omic_type="{omic_type}", shadow=s, n_tails=n) for d in EVAL_DATASETS for s in SHADOWS for n in N_TAILS]
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz")
    params:
        omic_type = "{omic_type}"
    run:
        import pandas as pd
    
        evaluation = pd.concat([pd.read_table(f) for f in input.evaluations])
        evaluation["omic_type"] = params.omic_type
        
        evaluation.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
        
        
# rule compute_signature_within:
#     input:
#         data = lambda wildcards: DATASETS[wildcards.dataset][wildcards.omic_type]
#     output:
#         signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz")
#     run:
#         import pandas as pd
        
#         data = pd.read_table(input.data, index_col=0)
        
#         # subtract median
#         signature = data
#         signature = signature - signature["ACH-001086"].values.reshape(-1,1)
#         #signature = signature - signature.median(axis=1).values.reshape(-1,1)
        
#         # save
#         signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
#         print("Done!")
        
        
rule compute_protein_activity:
    input:
        #signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"),
        signature = lambda wildcards: PERT_FILES[wildcards.omic_type][wildcards.dataset],
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_type]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
        
rule figures_network_evaluation:
    input:
        evaluation_ex = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz"),
        evaluation_genexpr = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-genexpr.tsv.gz"),
        protein_activity_ex = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_K562-EX.tsv.gz"),
        protein_activity_genexpr = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_K562-genexpr.tsv.gz")
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
       
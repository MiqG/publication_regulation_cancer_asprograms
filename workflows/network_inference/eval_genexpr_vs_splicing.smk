import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

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

REGULON_SETS = [
    "experimentally_derived_regulons_pruned",
    "aracne_regulons_development",
    "mlr_regulons_development",
    "aracne_and_experimental_regulons",
    "mlr_and_experimental_regulons",
    "aracne_and_mlr_regulons"
]

TOP_N = [100, 90, 80, 70, 60, 50, 40]
ROBUSTNESS_EVAL_SETS = ["top{N}_experimentally_derived_regulons_pruned".format(N=n) for n in TOP_N]
REGULON_SETS = REGULON_SETS + ROBUSTNESS_EVAL_SETS

SHADOWS = ["no"] # bug in viper does not allow shadow correction
N_TAILS = ["one","two"]

##### RULES #####
rule all:
    input:
        # prepare evaluation labels
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"),
        
        # evaluate regulons
        ## run
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{regulon_set}-{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz"), regulon_set=REGULON_SETS, dataset=EVAL_DATASETS, omic_type=OMIC_TYPES, shadow=SHADOWS, n_tails=N_TAILS),
        ## merge
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # make figures
        os.path.join(RESULTS_DIR,"figures","regulon_evaluation"),
        os.path.join(RESULTS_DIR,"figures","regulon_inference")
        
        
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
        regulons = os.path.join(RESULTS_DIR,"files","{regulon_set}-{omic_type}"),
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels","{dataset}.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{regulon_set}-{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz")
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
        evaluations = [os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{regulon_set}-{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz").format(regulon_set=r, dataset=d, omic_type="{omic_type}", shadow=s, n_tails=n) for r in REGULON_SETS for d in EVAL_DATASETS for s in SHADOWS for n in N_TAILS]
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
        
    
rule figures_regulon_evaluation:
    input:
        evaluation_ex = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz"),
        evaluation_genexpr = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-genexpr.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","regulon_evaluation"))
    shell:
        """
        Rscript scripts/figures_regulon_evaluation.R \
                    --evaluation_ex_file={input.evaluation_ex} \
                    --evaluation_genexpr_file={input.evaluation_genexpr} \
                    --figs_dir={output}
        """
        
        
rule figures_inference_troubleshooting:
    input:
        experimental_pruned_path = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX"),
        aracne_and_experimental_path = os.path.join(RESULTS_DIR,"files","aracne_and_experimental_regulons-EX"),
        mlr_and_experimental_path = os.path.join(RESULTS_DIR,"files","mlr_and_experimental_regulons-EX")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","regulon_inference"))
    shell:
        """
        Rscript scripts/figures_inference_troubleshooting.R \
                    --experimental_pruned_path={input.experimental_pruned_path} \
                    --aracne_and_experimental_path={input.aracne_and_experimental_path} \
                    --mlr_and_experimental_path={input.mlr_and_experimental_path} \
                    --figs_dir={output}
        """
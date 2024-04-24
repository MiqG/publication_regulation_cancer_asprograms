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

PREP_VIPER_DIR = os.path.join(os.path.dirname(ROOT),"publication_viper_splicing","data","prep")

OMIC_TYPES = ["genexpr","scgenexpr"]

REGULON_DIRS = {
    "genexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-genexpr"),
    "scgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr")
}

##### RULES #####
rule all:
    input:
        # detection of tumorigenesis with single cell networks
        ## signature
        os.path.join(RESULTS_DIR,"files","signatures","tumorigenesis-genexpr.tsv.gz"),
        ## compute protein activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),

        # estimate splicing factor activity for K562 in bulk and single cell
        # expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-genexpr.tsv.gz"), dataset=PERT_GENEXPR_FILES.keys()),
        
        
        # make figures
        # os.path.join(RESULTS_DIR,"figures","eval_bulk_vs_singlecell"),
        
        
rule compute_signatures:
    input:
        metadata = os.path.join(PREP_VIPER_DIR,'metadata','tumorigenesis.tsv.gz'),
        splicing = os.path.join(PREP_VIPER_DIR,'genexpr_tpm','tumorigenesis.tsv.gz'),
    output:
        signatures = os.path.join(RESULTS_DIR,"files","signatures","tumorigenesis-genexpr.tsv.gz")
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # load
        metadata = pd.read_table(input.metadata)
        splicing = pd.read_table(input.splicing, index_col=0)
        
        # subset
        common_samples = set(metadata["sampleID"]).intersection(splicing.columns)
        metadata = metadata.loc[metadata["sampleID"].isin(common_samples)].copy()
        splicing = splicing[common_samples].copy()
        
        # delta PSI as the difference between conditions and the mean of the conditions
        signatures = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                
                # there may be empty controls
                if any(splicing.columns.isin(ctls)):
                    psi_ctls = splicing.loc[:,splicing.columns.isin(ctls)].mean(axis=1)

                    # compute delta PSI
                    dpsi = splicing[sample_oi] - psi_ctls

                    signatures[sample_oi] = dpsi

                    del dpsi, psi_ctls, ctls
                else:
                    print(splicing.columns, ctls)
                    print(splicing.columns.isin(ctls))
                    continue

        signatures = pd.DataFrame(signatures)
        
        # save
        signatures.reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","tumorigenesis-genexpr.tsv.gz"),
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_type]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-{omic_type}.tsv.gz")
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
        protein_activity_bulk = os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-genexpr.tsv.gz"),
        protein_activity_singlecell = os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-scgenexpr.tsv.gz"),
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

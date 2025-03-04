import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","sf_programs_in_differentiation")
DATASET_DIR = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","datasets")
REGULONS_PATH = os.path.join(ROOT,"results","new_empirical_network","files","experimentally_derived_regulons_pruned_w_viper_networks-EX")
DRIVER_TYPES_PATH = os.path.join(ROOT,"results","new_empirical_network",'files','PANCAN','cancer_program.tsv.gz')

METADATA_FILES = {
    "CardosoMoreira2020": os.path.join(DATASET_DIR,'metadata','CardosoMoreira2020.tsv.gz')
}

SPLICING_FILES = {
    "CardosoMoreira2020": os.path.join(DATASET_DIR,'event_psi','CardosoMoreira2020-EX.tsv.gz')
}

GENEXPR_FILES = {
    "CardosoMoreira2020": os.path.join(DATASET_DIR,'genexpr_tpm','CardosoMoreira2020.tsv.gz')
}

DATA_FILES = {
    "EX": SPLICING_FILES,
    "genexpr": GENEXPR_FILES
}

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}
DATASETS = ["CardosoMoreira2020"]
OMIC_TYPES = ["EX","genexpr"]

##### RULES #####
rule all:
    input:
        # signatures
        expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),

        # protein activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX.tsv.gz"), dataset=DATASETS),

        # figures
        os.path.join(RESULTS_DIR,"figures","CardosoMoreira2020")
    
    
rule compute_signatures:
    input:
        metadata = lambda wildcards: METADATA_FILES[wildcards.dataset],
        splicing = lambda wildcards: DATA_FILES[wildcards.omic_type][wildcards.dataset],
    output:
        signatures = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz")
    params:
        dataset = "{dataset}"
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # load
        metadata = pd.read_table(input.metadata)
        splicing = pd.read_table(input.splicing, index_col=0)
        dataset = params.dataset
        
        # subset
        common_samples = list(set(metadata["sampleID"]).intersection(splicing.columns))
        metadata = metadata.loc[metadata["sampleID"].isin(common_samples)].copy()
        splicing = splicing[common_samples].copy()
        
        signatures = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            if dataset=="CardosoMoreira2020":
                tissue = metadata.loc[metadata["sampleID"]==sample_oi,"sample_title"]
                tissue = tissue.values[0].split(".")[2]
                ctls = metadata.loc[metadata["sample_title"].str.contains(tissue),"sampleID"]
            
            # controls will be np.nan
            if any(splicing.columns.isin(ctls)):
                psi_ctls = splicing.loc[:,splicing.columns.isin(ctls)].mean(axis=1)

                # compute delta PSI
                dpsi = splicing[sample_oi] - psi_ctls

                signatures[sample_oi] = dpsi

                del dpsi, psi_ctls, ctls
            else:
                print("Skipping %s" % sample_oi)

        signatures = pd.DataFrame(signatures)
        
        # save
        signatures.reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-EX.tsv.gz"),
        regulons_path = REGULONS_PATH
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
       """
        
rule figures_CardosoMoreira2020:
    input:
        metadata = os.path.join(DATASET_DIR,"metadata","CardosoMoreira2020.tsv.gz"),
        genexpr = os.path.join(DATASET_DIR,"genexpr_tpm","CardosoMoreira2020.tsv.gz"),
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","CardosoMoreira2020-EX.tsv.gz"),
        driver_types = DRIVER_TYPES_PATH
    output:
        directory(os.path.join(RESULTS_DIR,"figures","CardosoMoreira2020"))
    shell:
        """
        Rscript scripts/figures_CardosoMoreira2020.R \
                    --genexpr_file={input.genexpr} \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --driver_types_file={input.driver_types} \
                    --figs_dir={output}
        """    
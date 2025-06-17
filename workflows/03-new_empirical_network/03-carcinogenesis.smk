import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
SRC_DIR = os.path.join(ROOT,"src")
RESULTS_DIR = os.path.join(ROOT,"results","new_empirical_network")
DATASET_DIR = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","datasets")
REGULONS_PATH = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned_w_viper_networks-EX")


# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}
OMIC_TYPES = ["EX"]
DATASETS = ["carcinogenesis"]

###### RULES ######
rule all:
    input:
        # compute signatures
        expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # compute protein activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","carcinogenesis")
        

rule compute_signatures:
    input:
        metadata = os.path.join(DATASET_DIR,'metadata','tumorigenesis.tsv.gz'),
        splicing = os.path.join(DATASET_DIR,'event_psi','tumorigenesis-{omic_type}.tsv.gz'),
    output:
        signatures = os.path.join(RESULTS_DIR,"files","signatures","carcinogenesis-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # load
        metadata = pd.read_table(input.metadata)
        splicing = pd.read_table(input.splicing, index_col=0)
        
        # metadata - only data from Darmanis
        metadata = metadata.loc[metadata["study_accession"]=="PRJNA193487"]
        
        # subset
        common_samples = list(set(metadata["sampleID"]).intersection(splicing.columns))
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
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"),
        regulons_path = REGULONS_PATH
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

        
rule figures_tumorigenesis:
    input:
        metadata = os.path.join(DATASET_DIR,"metadata","tumorigenesis.tsv.gz"),
        genexpr = os.path.join(DATASET_DIR,"genexpr_tpm","tumorigenesis.tsv.gz"),
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","carcinogenesis-EX.tsv.gz"),
        driver_types = os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
    output:
        directory(os.path.join(RESULTS_DIR,"figures","carcinogenesis"))
    shell:
        """
        Rscript scripts/figures_carcinogenesis.R \
                    --genexpr_file={input.genexpr} \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --driver_types_file={input.driver_types} \
                    --figs_dir={output}
        """
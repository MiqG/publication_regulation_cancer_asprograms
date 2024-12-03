import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","activity_estimation_w_genexpr")
DATASET_DIR = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","datasets")

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}
REGULON_DIRS = {
    "bulkgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-bulkgenexpr"),
    "scgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr"),
    "bulkscgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-bulkscgenexpr")
}

OMIC_GENEXPR_REGULONS = ["bulkgenexpr","scgenexpr","bulkscgenexpr"]
MODEL_ARCHS = ["fclayer","ewlayer"]
K_CROSS_VALIDATION = 5
DATASETS = ["carcinogenesis","Hodis2022-invitro_eng_melanoc"]
##### RULES #####
rule all:
    input:
        # comppute genexpr signatures carcinogenesis
        os.path.join(RESULTS_DIR,"files","signatures","carcinogenesis-genexpr.tsv.gz"),
        os.path.join(RESULTS_DIR,"files","signatures","Hodis2022-invitro_eng_melanoc-genexpr.tsv.gz"),
        
        # compute activity carcinogenesis with genexpr networks
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, dataset=DATASETS),
        
        # adjust activities from genexpr networks
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}-adjusted_{model_type}.tsv.gz"), model_type=MODEL_ARCHS, omic_regulon=OMIC_GENEXPR_REGULONS, dataset=DATASETS),
        
        # combine activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-merged.tsv.gz"), dataset=DATASETS),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","eval_carcinogenesis")
        
        
rule compute_signatures_bulk:
    input:
        metadata = os.path.join(DATASET_DIR,'metadata','tumorigenesis.tsv.gz'),
        genexpr = os.path.join(DATASET_DIR,'genexpr_tpm','tumorigenesis.tsv.gz'),
    output:
        signatures = os.path.join(RESULTS_DIR,"files","signatures","carcinogenesis-genexpr.tsv.gz")
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # load
        metadata = pd.read_table(input.metadata)
        genexpr = pd.read_table(input.genexpr, index_col=0)
        
        # metadata - only data from Darmanis
        metadata = metadata.loc[metadata["study_accession"]=="PRJNA193487"]
        
        # subset
        common_samples = list(set(metadata["sampleID"]).intersection(genexpr.columns))
        metadata = metadata.loc[metadata["sampleID"].isin(common_samples)].copy()
        genexpr = genexpr[common_samples].copy()
        
        # log2 fold change as the difference between conditions and the mean of the conditions
        signatures = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                
                # there may be empty controls
                if any(genexpr.columns.isin(ctls)):
                    tpm_ctls = genexpr.loc[:,genexpr.columns.isin(ctls)].mean(axis=1)

                    # compute delta PSI
                    fctpm = genexpr[sample_oi] - tpm_ctls

                    signatures[sample_oi] = fctpm

                    del fctpm, tpm_ctls, ctls
                else:
                    print(genexpr.columns, ctls)
                    print(genexpr.columns.isin(ctls))
                    continue

        signatures = pd.DataFrame(signatures)
        
        # save
        signatures.reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_signatures_singlecell:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-pseudobulk.h5ad")
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","Hodis2022-invitro_eng_melanoc-genexpr.tsv.gz")
    params:
        dataset = "Hodis2022-invitro_eng_melanoc"
    resources:
        memory = 20, # GB
        runtime = 3600*2 # h
    run:
        import pandas as pd
        import scanpy as sc
        import numpy as np
        from tqdm import tqdm
        
        # load
        adata = sc.read_h5ad(input.adata)
        dataset = params.dataset
        
        # compute fold changes by cell type
        if dataset=="Boiarsky2022-myeloma":
            # in this case all CTL cells are of CD138+ and cell type "normal"
            ctl_cells = adata.obs["is_ctl"]
            print("We have %s CTL cells." % ctl_cells.sum() )

            genexpr_ctl = adata[ctl_cells,:].to_df().mean(axis=0).values.reshape(1,-1)
            genexpr = adata.to_df()
            genexpr = genexpr - genexpr_ctl
            
        else:
            genexpr = []
            cell_types = adata.obs["cell_type"].unique()
            for cell_type_oi in tqdm(cell_types):

                ctl_cells = (adata.obs["cell_type"]==cell_type_oi) & adata.obs["is_ctl"]
                print("We have %s CTL cells for cell type %s." % (ctl_cells.sum(), cell_type_oi))

                genexpr_ctl = adata[ctl_cells,:].to_df().mean(axis=0).values.reshape(1,-1)
                genexpr_batch = adata[adata.obs["cell_type"]==cell_type_oi].to_df()
                genexpr_batch = genexpr_batch - genexpr_ctl
                genexpr.append(genexpr_batch)

            genexpr = pd.concat(genexpr)

        # save
        genexpr.T.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")    
        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-genexpr.tsv.gz"),
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_regulon]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
       """

rule adjust_activity:
    input:
        activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"),
        weights = [os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth").format(k=k, omic_regulon="{omic_regulon}", model_type="{model_type}") for k in range(K_CROSS_VALIDATION)],
        input_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","input_regulators.tsv.gz"),
        output_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","output_regulators.tsv.gz")
    params:
        weights = ",".join([os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth").format(k=k, omic_regulon="{omic_regulon}", model_type="{model_type}") for k in range(K_CROSS_VALIDATION)]),
        model_type = "{model_type}",
        script_dir = SRC_DIR
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}-adjusted_{model_type}.tsv.gz")
    shell:
        """
        python {params.script_dir}/vipersp/scripts/adjust_genexpr_sf_activity.py \
                    --activity_file={input.activity} \
                    --weights_files="{params.weights}" \
                    --input_regulators_file={input.input_regulators} \
                    --output_regulators_file={input.output_regulators} \
                    --model_type={params.model_type} \
                    --output_file={output}
        """
        
rule combine_protein_activity:
    input:
        protein_activity_ex = lambda wildcards: [os.path.join(ROOT,"results","new_empirical_network","files","protein_activity","carcinogenesis-EX.tsv.gz")] if wildcards.dataset=="carcinogenesis" else [],
        protein_activity_genexpr = [os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz").format(omic_regulon=o, dataset="{dataset}") for o in OMIC_GENEXPR_REGULONS],
        protein_activity_genexpr_adjusted = [os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}-adjusted_{model_type}.tsv.gz").format(omic_regulon=o, model_type=m, dataset="{dataset}") for o in OMIC_GENEXPR_REGULONS for m in MODEL_ARCHS]
    output:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-merged.tsv.gz")
    run:
        import os
        import pandas as pd
        
        files = input.protein_activity_ex + input.protein_activity_genexpr + input.protein_activity_genexpr_adjusted
        
        protein_activity = []
        for f in files:
            df = pd.read_table(f).melt(id_vars=["regulator"])
            df["dataset_id"] = os.path.basename(f).replace(".tsv.gz","")
            df = df.rename(columns={"variable":"sampleID", "value":"activity"})
            protein_activity.append(df)
            
        protein_activity = pd.concat(protein_activity)
        
        protein_activity.to_csv(output.protein_activity, **SAVE_PARAMS)
        
        print("Done!")
    
rule figures_eval_carcinogenesis:
    input:
        metadata_bulk = os.path.join(DATASET_DIR,"metadata","tumorigenesis.tsv.gz"),
        metadata_singlecell = os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-conditions.tsv.gz"),
        protein_activity_bulk = os.path.join(RESULTS_DIR,"files","protein_activity","carcinogenesis-merged.tsv.gz"),
        protein_activity_singlecell = os.path.join(RESULTS_DIR,"files","protein_activity","Hodis2022-invitro_eng_melanoc-merged.tsv.gz"),
        driver_types = os.path.join(ROOT,"results","new_empirical_network",'files','PANCAN','cancer_program.tsv.gz')
    output:
        directory(os.path.join(RESULTS_DIR,"figures","eval_carcinogenesis"))
    shell:
        """
        Rscript scripts/figures_eval_carcinogenesis.R \
                    --protein_activity_bulk_file={input.protein_activity_bulk} \
                    --protein_activity_singlecell_file={input.protein_activity_singlecell} \
                    --metadata_bulk_file={input.metadata_bulk} \
                    --metadata_singlecell_file={input.metadata_singlecell} \
                    --driver_types_file={input.driver_types} \
                    --figs_dir={output}
        """
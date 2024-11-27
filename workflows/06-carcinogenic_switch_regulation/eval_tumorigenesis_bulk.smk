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
PREP_VIPER_DIR = os.path.join(os.path.dirname(ROOT),"publication_viper_splicing","data","prep")

OMIC_TYPES = ["genexpr","scgenexpr","EX"]
OMIC_GENEXPR_TYPES = ["genexpr","scgenexpr"]

REGULON_DIRS = {
    "genexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-genexpr"),
    "EX": os.path.join(VIPER_SPLICING_DIR,"data","empirical_sf_networks-EX"),
    "scgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr")
}
MODEL_TYPES = ["fclayer","ewlayer"]

DATASETS = {
    "EX": os.path.join(PREP_VIPER_DIR,'event_psi','tumorigenesis-EX.tsv.gz'),
    "genexpr": os.path.join(PREP_VIPER_DIR,'genexpr_tpm','tumorigenesis.tsv.gz')
}

OMIC_PERT_DICT = {
    "EX": "EX",
    "genexpr": "genexpr",
    "scgenexpr": "genexpr"
}

ONTOLOGIES = ["hallmarks","reactome"]

##### RULES #####
rule all:
    input:
        # detection of tumorigenesis with single cell networks
        ## signature
        expand(os.path.join(RESULTS_DIR,"files","signatures","tumorigenesis-{omic_signature}.tsv.gz"), omic_signature=["EX","genexpr"]),
        ## GSEA
        expand(os.path.join(RESULTS_DIR,"files","gsea","tumorigenesis-{omic_signature}-{ontology_oi}.tsv.gz"), ontology_oi=ONTOLOGIES, omic_signature=["genexpr"]),
        ## compute protein activities
        ### classic
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-{omic_regulon}.tsv.gz"), omic_regulon=REGULON_DIRS.keys()),
        ### models
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz"), model_type=MODEL_TYPES, omic_regulon=OMIC_GENEXPR_TYPES)
        # make figures
        # os.path.join(RESULTS_DIR,"figures","eval_tumorigenesis"),
        
        
rule compute_signatures:
    input:
        metadata = os.path.join(PREP_VIPER_DIR,'metadata','tumorigenesis.tsv.gz'),
        splicing = lambda wildcards: DATASETS[wildcards.omic_signature],
    output:
        signatures = os.path.join(RESULTS_DIR,"files","signatures","tumorigenesis-{omic_signature}.tsv.gz")
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
        
rule run_gsea:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","tumorigenesis-{omic_signature}.tsv.gz"),
        msigdb_dir = os.path.join(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs"),
        gene_info = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","gsea","tumorigenesis-{omic_signature}-{ontology_oi}.tsv.gz")
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
        signature = lambda wildcards: os.path.join(RESULTS_DIR,"files","signatures","tumorigenesis-{omic_signature}.tsv.gz").format(omic_signature=OMIC_PERT_DICT[wildcards.omic_regulon]),
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_regulon]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-{omic_regulon}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
        
rule predict_sf_activity_from_model:
    input:
        activity = os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-{omic_regulon}.tsv.gz"),
        weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"),
        genes_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","genes_train.tsv.gz"),
        regulators_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","regulators_train.tsv.gz")
    output:
        activity_pred = os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz")
    params:
        model_type = "{model_type}"        
    run:
        import torch
        import pandas as pd
        from vipersp.model import EWlayer, FClayer
        
        # load
        activity = pd.read_table(input.activity, index_col=0)
        weights = torch.load(input.weights)
        genes_train = list(pd.read_table(input.genes_train, header=None)[0])
        regulators_train = list(pd.read_table(input.regulators_train, header=None)[0])
        model_type = params.model_type
        
        # prep data
        activity = pd.merge(
            pd.DataFrame(index=genes_train),
            activity,
            how="left", left_index=True, right_index=True
        )
        X = torch.tensor(activity.fillna(0).T.values, dtype=torch.float32)
        
        if model_type=="fclayer":
            model = FClayer(input_size=len(genes_train), output_size=len(regulators_train))
        
        elif model_type=="ewlayer":
            model = EWlayer(input_size=len(genes_train))        
            
        model.load_state_dict(weights)        
        
        # make predictions
        model.eval()
        with torch.no_grad():
            Y_hat = model(X)
        
        # prep outputs
        activity_pred = pd.DataFrame(Y_hat.detach().numpy().T, index=regulators_train, columns=activity.columns)
        activity_pred.index.name = "regulator"
        
        # save
        activity_pred.reset_index().to_csv(output.activity_pred, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule figures_regulon_evaluation:
    input:
        protein_activity_bulk = os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-genexpr.tsv.gz"),
        protein_activity_singlecell = os.path.join(RESULTS_DIR,"files","protein_activity","tumorigenesis-scgenexpr.tsv.gz"),
        cancer_program = os.path.join(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","eval_tumorigenesis"))
    shell:
        """
        Rscript scripts/figures_eval_tumorigenesis.R \
                    --protein_activity_bulk_file={input.protein_activity_bulk} \
                    --protein_activity_singlecell_file={input.protein_activity_singlecell} \
                    --cancer_program_file={input.cancer_program} \
                    --figs_dir={output}
        """

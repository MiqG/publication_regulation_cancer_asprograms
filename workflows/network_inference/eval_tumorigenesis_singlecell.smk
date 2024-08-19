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

OMIC_TYPES = ["genexpr","scgenexpr"]
OMIC_GENEXPR_TYPES = ["genexpr","scgenexpr"]

REGULON_DIRS = {
    "genexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-genexpr"),
    "scgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr")
}
MODEL_TYPES = ["fclayer","ewlayer"]

DATASETS = ["Hodis2022-invitro_eng_melanoc","Becker2021-adenoma","Boiarsky2022-myeloma"]

ONTOLOGIES = ["hallmarks","reactome"]

##### RULES #####
rule all:
    input:
        # detection of tumorigenesis with single cell networks
        ## signature
        expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-genexpr.tsv.gz"), dataset=DATASETS),
        ## GSEA
        expand(os.path.join(RESULTS_DIR,"files","gsea","{dataset}-{ontology_oi}.tsv.gz"), dataset=DATASETS[:1], ontology_oi=ONTOLOGIES),
        ## compute protein activities
        ### classic
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"), dataset=DATASETS, omic_regulon=REGULON_DIRS.keys()),
        ### models
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz"), dataset=DATASETS, model_type=MODEL_TYPES, omic_regulon=OMIC_GENEXPR_TYPES),
        # make figures
        expand(os.path.join(RESULTS_DIR,"figures","eval_tumorigenesis_singlecell-{dataset}"), dataset=DATASETS),
        
        
rule compute_signatures:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","{dataset}-pseudobulk.h5ad")
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-genexpr.tsv.gz")
    params:
        dataset = "{dataset}"
    resources:
        memory = 20, # GB
        runtime = 3600*2 # h
    run:
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
        
rule run_gsea:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-genexpr.tsv.gz"),
        msigdb_dir = os.path.join(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs"),
        gene_info = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","gsea","{dataset}-{ontology_oi}.tsv.gz")
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
        
rule predict_sf_activity_from_model:
    input:
        activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"),
        weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"),
        genes_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","genes_train.tsv.gz"),
        regulators_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","regulators_train.tsv.gz")
    output:
        activity_pred = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz")
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
        
        
rule figures_eval_tumorigenesis_singlecell:
    input:
        protein_activity_genexpr = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-genexpr.tsv.gz"),
        protein_activity_scgenexpr = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-scgenexpr.tsv.gz"),
        protein_activity_ew_model_genexpr = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_ewlayer_and_genexpr.tsv.gz"),
        protein_activity_ew_model_scgenexpr = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_ewlayer_and_scgenexpr.tsv.gz"),
        protein_activity_fc_model_genexpr = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_fclayer_and_genexpr.tsv.gz"),
        protein_activity_fc_model_scgenexpr = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_fclayer_and_scgenexpr.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"singlecell","{dataset}-conditions.tsv.gz"),
        cancer_program = os.path.join(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
    params:
        dataset = "{dataset}"
    output:
        directory(os.path.join(RESULTS_DIR,"figures","eval_tumorigenesis_singlecell-{dataset}"))
    shell:
        """
        Rscript scripts/figures_eval_tumorigenesis_singlecell.R \
                    --protein_activity_genexpr_file={input.protein_activity_genexpr} \
                    --protein_activity_scgenexpr_file={input.protein_activity_scgenexpr} \
                    --protein_activity_ew_model_genexpr_file={input.protein_activity_ew_model_genexpr} \
                    --protein_activity_ew_model_scgenexpr_file={input.protein_activity_ew_model_scgenexpr} \
                    --protein_activity_fc_model_genexpr_file={input.protein_activity_fc_model_genexpr} \
                    --protein_activity_fc_model_scgenexpr_file={input.protein_activity_fc_model_scgenexpr} \
                    --cancer_program_file={input.cancer_program} \
                    --metadata_file={input.metadata} \
                    --dataset={params.dataset} \
                    --figs_dir={output}
        """

import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
NETWORKS_DIR = os.path.join(ROOT,"results","network_inference")
RESULTS_DIR = os.path.join(ROOT,"results","program_regulation")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

PERT_GENEXPR_FILES = {
    "ReplogleWeissman2022_rpe1": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-log2_fold_change_cpm.tsv.gz"),
    #"ReplogleWeissman2022_K562_gwps": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_gwps-log2_fold_change_cpm.tsv.gz")
    #"ReplogleWeissman2022_K562_essential": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-log2_fold_change_cpm.tsv.gz")
}
MODEL_TYPES = ["fclayer"]
OMIC_REGULONS = ["scgenexpr"]

##### RULES #####
rule all:
    input:
        # estimate splicing factor activity
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-scgenexpr.tsv.gz"), dataset=PERT_GENEXPR_FILES.keys()),
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz"), dataset=PERT_GENEXPR_FILES.keys(), model_type=MODEL_TYPES, omic_regulon=OMIC_REGULONS)
        
        # make figures
        #os.path.join(RESULTS_DIR,"figures","upstream_regulators"),
        
        
        
rule compute_protein_activity:
    input:
        signature = lambda wildcards: PERT_GENEXPR_FILES[wildcards.dataset],
        regulons_path = os.path.join(NETWORKS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-scgenexpr.tsv.gz")
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
        weights = os.path.join(NETWORKS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"),
        genes_train = os.path.join(NETWORKS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","genes_train.tsv.gz"),
        regulators_train = os.path.join(NETWORKS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","regulators_train.tsv.gz")
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
        
        
        
rule figures_upstream_regulators:
    input:
        protein_activity_rpe1 = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_rpe1-genexpr.tsv.gz"),
        protein_activity_k562 = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_K562_gwps-genexpr.tsv.gz"),
        gene_info = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz"),
        cancer_program = os.path.join(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz"),
        msigdb_dir = os.path.join(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs"),
        cosmic_genes = os.path.join(RAW_DIR,"COSMIC","cancer_gene_census.tsv")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","upstream_regulators"))
    shell:
        """
        Rscript scripts/figures_upstream_regulators.R \
                    --protein_activity_rpe1_file={input.protein_activity_rpe1} \
                    --protein_activity_k562_file={input.protein_activity_k562} \
                    --gene_info_file={input.gene_info} \
                    --cancer_program_file={input.cancer_program} \
                    --msigdb_dir={input.msigdb_dir} \
                    --cosmic_genes_file={input.cosmic_genes} \
                    --figs_dir={output}
        """

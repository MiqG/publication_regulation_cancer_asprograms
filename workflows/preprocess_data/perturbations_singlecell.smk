"""
Author: Miquel Anglada Girotto
Contact: miquelangladagirotto [at] gmail [dot] com
"""

import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SRC_DIR = os.path.join(ROOT,"src")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","preprocess_data")

SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

DATASETS = ["K562_essential","K562_gwps","rpe1"]

##### RULES #####
rule all:
    input:
        # ReplogleWeissman2022
        ## preprocess
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad"), dataset=DATASETS),
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-cell_summary.tsv.gz"), dataset=DATASETS),
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-gene_summary.tsv.gz"), dataset=DATASETS),
        ## summarize
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-pseudobulk.h5ad"), dataset=DATASETS),
        ## fold change
        expand(os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_{dataset}-log2_fold_change_cpm.tsv.gz"), dataset=DATASETS)
        
        
rule preprocess_scperturb:
    input:
        adata = os.path.join(RAW_DIR,"scPerturb","ReplogleWeissman2022_{dataset}.h5ad"),
    output:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad"),
        cell_summary = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-cell_summary.tsv.gz"),
        gene_summary = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-gene_summary.tsv.gz")
    resources:
        memory = 300, # GB
        runtime = 3600*2 # h
    run:
        import scanpy as sc
        import gc
        
        # load
        adata = sc.read_h5ad(input.adata)
        gc.collect()
        
        # prep
        ## metadata
        ### column with gene name and ensembl of perturbed gene
        adata.obs["PERT_ENSEMBL"] = adata.obs["gene_id"]
        adata.obs["PERT_GENE"] = adata.obs["gene"]
        
        ## gene expression
        ### use ensembl as gene identifier
        adata.var["gene_name"] = adata.var.index
        adata.var = adata.var.set_index("ensembl_id")
        ### filter out cells
        adata = adata[adata.obs["ngenes"]>=200,:]
        ### filter out genes
        adata = adata[:,adata.var["ncells"]>=3]
        
        # save
        adata.write(output.adata)
        adata.var.reset_index().to_csv(output.gene_summary, **SAVE_PARAMS)
        adata.obs.reset_index().to_csv(output.cell_summary, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule summarize_genexpr_by_perturbation:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad"),
    output:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-pseudobulk.h5ad")
    resources:
        memory = 300, # GB
        runtime = 3600*2 # h
    run:
        import scanpy as sc
        import numpy as np
        from tqdm import tqdm
        
        # load
        adata = sc.read_h5ad(input.adata)
        
        # prep
        ## normalize counts to CPMs
        adata.X = 1e6 * (adata.X / np.sum(adata.X, axis=0).reshape(1,-1))
        ## log scale
        adata.X = np.log2(adata.X + 1)
        
        # summarize gene expression by perturbation
        def summarize_genexpr_perts(adata):
            conditions = adata.obs["PERT_ENSEMBL"].value_counts()
            replicated_conditions = conditions.loc[conditions>1]
            unique_conditions = conditions.loc[conditions==1]

            # add non unique conditions
            idx = adata.obs["PERT_ENSEMBL"].isin(unique_conditions.index)
            genexpr = adata[idx,:].to_df()
            genexpr = genexpr.rename(index=adata[idx,:].obs["PERT_ENSEMBL"].to_dict())
            
            # add replicated conditions
            summarized_genexpr = {}
            for condition_oi in tqdm(replicated_conditions.index):
                cells_oi = adata.obs[adata.obs["PERT_ENSEMBL"]==condition_oi].index
                summarized_genexpr[condition_oi] = adata[cells_oi,:].to_df().mean(axis=0)
            summarized_genexpr = pd.DataFrame(summarized_genexpr).T

            # prepare outputs
            summarized_genexpr = pd.concat([summarized_genexpr, genexpr], axis=0)

            return summarized_genexpr
        
        genexpr = summarize_genexpr_perts(adata)

        # make adata
        adata = sc.AnnData(genexpr)
        
        # save
        adata.write(output.adata)
        
        print("Done!")
        
        
rule compute_signatures:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-pseudobulk.h5ad")
    output:
        signature = os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_{dataset}-log2_fold_change_cpm.tsv.gz")
    resources:
        memory = 20, # GB
        runtime = 3600*2 # h
    run:
        import scanpy as sc
        
        # load
        adata = sc.read_h5ad(input.adata)
        
        # compute fold changes
        adata.X = adata.X - adata["non-targeting"].X.reshape(1,-1)
        
        # save
        adata.to_df().T.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
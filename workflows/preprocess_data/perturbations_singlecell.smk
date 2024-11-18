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
#PSEUDOBULK_TYPES = ["pseudobulk_by_batch","pseudobulk_across_batches"]
PSEUDOBULK_TYPES = ["pseudobulk_across_batches"]

##### RULES #####
rule all:
    input:
        # ReplogleWeissman2022
        ## preprocess
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad"), dataset=DATASETS),
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-cell_summary.tsv.gz"), dataset=DATASETS),
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-gene_summary.tsv.gz"), dataset=DATASETS),
        ## summarize
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-{pseudobulk_type}.h5ad"), dataset=DATASETS, pseudobulk_type=PSEUDOBULK_TYPES),
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-{pseudobulk_type}-conditions.tsv.gz"), dataset=DATASETS, pseudobulk_type=PSEUDOBULK_TYPES),
        # fold change
        expand(os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_{dataset}-{pseudobulk_type}-log2_fold_change_cpm.tsv.gz"), dataset=DATASETS, pseudobulk_type=PSEUDOBULK_TYPES)
        
        
rule preprocess_scperturb:
    input:
        adata = os.path.join(RAW_DIR,"scPerturb","ReplogleWeissman2022_{dataset}.h5ad")
    output:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad"),
        cell_summary = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-cell_summary.tsv.gz"),
        gene_summary = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-gene_summary.tsv.gz")
    params:
        thresh_ngenes = 200,
        thresh_ncells = 3,
        thresh_mito = 0.15
    resources:
        memory = 150, # GB
        runtime = 3600*2 # h
    run:
        import scanpy as sc
        import gc
        import numpy as np
        
        # load
        adata = sc.read_h5ad(input.adata)
        thresh_ngenes = params.thresh_ngenes
        thresh_ncells = params.thresh_ncells
        thresh_mito = params.thresh_mito
        gc.collect()
        
        # QC metrics
        ## total read counts per cell
        adata.obs["total_read_count"] = np.sum(adata.X, axis=1)
        ## what percentage of the total reads per cell belongs to mitochondrial genes?
        adata.var["is_mitochondrial"] = adata.var_names.str.startswith("MT-")
        adata.obs["total_mitochondrial_reads"] = np.sum(adata[:,adata.var["is_mitochondrial"]].X, axis=1)
        adata.obs["perc_mitochondrial_reads"] = adata.obs["total_mitochondrial_reads"] / adata.obs["total_read_count"]
        ## in how many cells is the gene expressed?
        adata.var["ncells_detected"] = np.sum(adata.X > 0, axis=0)
        ## how many genes are expressed per cell?
        adata.obs["ngenes_detected"] = np.sum(adata.X > 0, axis=1)
        
        # Perturb-seq specific
        ## add column with gene name and ensembl of perturbed gene
        adata.obs["PERT_ENSEMBL"] = adata.obs["gene_id"]
        adata.obs["PERT_GENE"] = adata.obs["gene"]
        adata.obs["condition"] = adata.obs["PERT_ENSEMBL"].astype(str)+"-"+adata.obs["batch"].astype(str)
        ## use ensembl as gene identifier
        adata.var["gene_name"] = adata.var.index
        adata.var = adata.var.set_index("ensembl_id")

        # filter out cells based on detected genes and mitochondrial reads
        to_keep = (adata.obs["ngenes_detected"]>=thresh_ngenes) & (adata.obs["perc_mitochondrial_reads"]<=thresh_mito)
        adata = adata[to_keep,:]
        
        # filter out genes based on detected cells
        to_keep = (adata.var["ncells_detected"]>=thresh_ncells)
        adata = adata[:,to_keep]
        
        # save
        adata.write(output.adata)
        adata.var.reset_index().to_csv(output.gene_summary, **SAVE_PARAMS)
        adata.obs.reset_index().to_csv(output.cell_summary, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule summarize_genexpr_by_perturbation:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad"),
    output:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-{pseudobulk_type}.h5ad"),
        conditions = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-{pseudobulk_type}-conditions.tsv.gz")
    resources:
        memory = 150, # GB
        runtime = 3600*2 # h
    params:
        pseudobulk_type = "{pseudobulk_type}"
    run:
        import scanpy as sc
        import numpy as np
        from tqdm import tqdm
        import gc
        
        # load
        adata = sc.read_h5ad(input.adata)
        pseudobulk_type = params.pseudobulk_type

        # # normalize counts to CPMs
        # ## Normalizing to median total counts
        # adata.X = 1e6 * (adata.X / adata.obs["total_read_count"].values.reshape(-1,1))
        # ## log scale
        # adata.X = np.log2(adata.X + 1)
        
        # summarize gene expression by some grouping variable(s)
        def summarize_genexpr_perts(adata, obs_oi, chunk_size, sep="___"):
            adata.obs["condition"] = adata.obs[obs_oi].apply(lambda x: sep.join(x.astype("str")), axis=1)
            conditions = adata.obs.value_counts("condition")
            replicated_conditions = conditions.loc[conditions>1]
            unique_conditions = conditions.loc[conditions==1]

            # add non unique conditions
            cells_oi = adata.obs.index[adata.obs["condition"].isin(unique_conditions.index)]
            genexpr = adata[cells_oi,:].to_df()
            genexpr = genexpr.rename(index=adata[cells_oi,:].obs["condition"].to_dict())
            
            # add replicated conditions in chunks
            chunk_size = chunk_size
            n_samples = len(replicated_conditions)
            n_chunks = int(np.ceil(n_samples/chunk_size))
            chunk_size = int(np.ceil(n_samples/n_chunks))
            
            print("Summing replicates...")
            summarized_genexpr = []
            for chunk_i in tqdm(range(n_chunks), total=n_chunks):
                conds_chunk = replicated_conditions[chunk_i*chunk_size:chunk_i*chunk_size+chunk_size]
                cells_oi = adata.obs[adata.obs["condition"].isin(conds_chunk.index)].index
                adata_chunk = adata[cells_oi].to_df()
                adata_chunk["condition"] = adata[cells_oi].obs["condition"]
                adata_chunk = adata_chunk.groupby(["condition"]).sum()
                
                # normalize pseudobulk gene counts
                adata_chunk = 1e6 * (adata_chunk / adata_chunk.sum(axis=1).values.reshape(-1,1))
                adata_chunk = np.log2(adata_chunk + 1)
                
                summarized_genexpr.append(adata_chunk)
                
                del adata_chunk
                gc.collect()
            summarized_genexpr = pd.concat(summarized_genexpr)

            # prepare outputs
            summarized_genexpr = pd.concat([summarized_genexpr, genexpr], axis=0)
            conditions.name = "n_cells"
            conditions = pd.DataFrame(conditions)
            
            # add condition information
            ## split obs_oi into new columns
            conditions[obs_oi] = conditions.reset_index()["condition"].str.split(sep, expand=True).values
            ## average expression of target gene
            ### we can only get the gene expression for the targeted gene if it was detected
            conditions["detected_pert"] = conditions["PERT_ENSEMBL"].isin(summarized_genexpr.columns)
            genes_oi = summarized_genexpr.columns.isin(conditions["PERT_ENSEMBL"])
            cells_oi = conditions.loc[conditions["detected_pert"]].index
            perts = summarized_genexpr.loc[cells_oi, genes_oi]
            ### prepare chunks
            chunk_size = 10_000
            n_samples = perts.shape[0]
            n_chunks = int(np.ceil(n_samples/chunk_size))
            chunk_size = int(np.ceil(n_samples/n_chunks))
            
            print("Getting expression of target gene...")
            genexpr_perts = []
            for chunk_i in tqdm(range(n_chunks), total=n_chunks):
                perts_chunk = perts.iloc[chunk_i*chunk_size:chunk_i*chunk_size+chunk_size,:]
                genes_pert = conditions.loc[perts_chunk.index,"PERT_ENSEMBL"].tolist()
                perts_chunk = pd.Series(
                    np.diagonal(perts_chunk.loc[:,genes_pert]),
                    index = perts_chunk.index
                )
                genexpr_perts.append(perts_chunk)

                del perts_chunk
                gc.collect()

            genexpr_perts = pd.concat(genexpr_perts)
            genexpr_perts.name = "genexpr_target"
            conditions = pd.concat([conditions, genexpr_perts], axis=1)
            
            ## average expression of target gene in corresponding control
            genes_oi = summarized_genexpr.columns.isin(conditions["PERT_ENSEMBL"])
            cells_oi = conditions.loc[conditions["PERT_ENSEMBL"]=="non-targeting"].index
            ctls = summarized_genexpr.loc[cells_oi, genes_oi].melt(ignore_index=False, value_name="genexpr_nontargeting")
            ctls[obs_oi] = ctls.reset_index()["condition"].str.split(sep, expand=True).values
            ctls = ctls.drop(columns="PERT_ENSEMBL").rename(columns={"ensembl_id":"PERT_ENSEMBL"})
            conditions = pd.merge(conditions.reset_index(), ctls, on=obs_oi, how="left")
            conditions = conditions.set_index("condition")
            
            ## fold change - knockdown efficiency
            conditions["pert_efficiency_fc"] = conditions["genexpr_target"] - conditions["genexpr_nontargeting"]
            
            # prepare outputs
            summarized_genexpr = summarized_genexpr.loc[conditions.index]
            
            return summarized_genexpr, conditions
        
        if pseudobulk_type=="pseudobulk_by_batch":
            obs_oi = ["PERT_ENSEMBL","batch"]
            chunk_size = 5_000
            
        elif pseudobulk_type=="pseudobulk_across_batches":
            obs_oi = ["PERT_ENSEMBL"]
            chunk_size = 50
        
        genexpr, conditions = summarize_genexpr_perts(adata, obs_oi, chunk_size)
        
        # make adata
        adata = sc.AnnData(X=genexpr, obs=conditions)
        
        # save
        adata.write(output.adata)
        conditions.reset_index().to_csv(output.conditions, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_signatures:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}-{pseudobulk_type}.h5ad")
    output:
        signature = os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_{dataset}-{pseudobulk_type}-log2_fold_change_cpm.tsv.gz")
    params:
        pseudobulk_type = "{pseudobulk_type}",
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
        pseudobulk_type = params.pseudobulk_type
        dataset = params.dataset
        
        # add info
        adata.obs["study_accession"] = "ReplogleWeissman2022_"+dataset
        adata.obs["cell_line"] = "K562" if "K562" in dataset else "RPE1"
        adata.obs["PERT_TYPE"] = "KNOCKDOWN"
        
        # compute fold changes by batch
        if pseudobulk_type=="pseudobulk_by_batch":
            all_batches = adata.obs["batch"].unique()
            ctl_batches = adata.obs.loc[adata.obs["PERT_ENSEMBL"]=="non-targeting","batch"].unique()
            print("We have CTL cells for %s out of %s batches." % (len(all_batches), len(ctl_batches)))
            
            genexpr = []
            for batch_i in tqdm(ctl_batches):
                ctl_cells = (adata.obs["PERT_ENSEMBL"]=="non-targeting") & (adata.obs["batch"]==batch_i)
                genexpr_ctl = adata[ctl_cells,:].X
                genexpr_batch = adata[adata.obs["batch"]==batch_i].to_df()
                genexpr_batch = genexpr_batch - genexpr_ctl
                genexpr.append(genexpr_batch)
            genexpr = pd.concat(genexpr)

        elif pseudobulk_type=="pseudobulk_across_batches":            
            ctl_cells = (adata.obs["PERT_ENSEMBL"]=="non-targeting")
            adata.obs["PERT_ID"] = adata.obs[
                ["study_accession","cell_line","PERT_ENSEMBL","PERT_TYPE"]
            ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
            adata.obs = adata.obs.set_index("PERT_ID")
            genexpr_ctl = adata[ctl_cells,:].X
            genexpr = adata[~ctl_cells,:].to_df()
            genexpr = genexpr - genexpr_ctl
            
        # save
        genexpr.T.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
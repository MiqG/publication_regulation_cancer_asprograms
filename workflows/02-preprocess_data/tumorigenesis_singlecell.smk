import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SRC_DIR = os.path.join(ROOT,"src")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","preprocess_data")

SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

DATASETS = ["Hodis2022-invitro_eng_melanoc"]

##### RULES #####
rule all:
    input:
        # preprocess
        ## Hodis2022
        os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc.h5ad"),
        os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-cell_summary.tsv.gz"),
        
        # summarize
        expand(os.path.join(PREP_DIR,"singlecell","{dataset}-pseudobulk.h5ad"), dataset=DATASETS),
        expand(os.path.join(PREP_DIR,"singlecell","{dataset}-conditions.tsv.gz"), dataset=DATASETS),
        
        # convert to tsv
        expand(os.path.join(PREP_DIR,"singlecell","{dataset}-pseudobulk.tsv.gz"), dataset=DATASETS)
        
        
rule preprocess_Hodis2022:
    input:
        genexpr = os.path.join(RAW_DIR,"articles","Hodis2022","invitro_eng_melanoc_logTP10K.txt.gz"),
        metadata = os.path.join(RAW_DIR,"articles","Hodis2022","invitro_invivo_all_metadatafile_mod_withCelltypes.csv"),
        gene_annot = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        adata = os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc.h5ad"),
        metadata = os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-cell_summary.tsv.gz")
    params:
        thresh_ngenes = 200,
        thresh_ncells = 3,
        thresh_mito = 0.15
    resources:
        memory = 150, # GB
        runtime = 3600*2 # h
    run:
        import pandas as pd
        import numpy as np
        import scanpy as sc
        import gc
        
        # load data
        genexpr = pd.read_table(input.genexpr)
        metadata = pd.read_csv(input.metadata)
        gene_annot = pd.read_table(input.gene_annot)
        gc.collect()
        
        # prep
        ## gene annotation
        gene_annot = gene_annot[["Approved symbol","Ensembl gene ID"]]
        gene_annot.columns = ["GENE","ENSEMBL"]
        
        ## gene names as ENSEMBL identifiers
        genexpr = pd.merge(genexpr, gene_annot, on="GENE", how="left")
        genexpr = genexpr.loc[~genexpr["ENSEMBL"].isnull()].copy()
        
        ## Transform log10 genexpr to log2
        genexpr = genexpr.drop(columns=["GENE"]).set_index("ENSEMBL").copy()
        genexpr = np.log2(10**genexpr)
        
        ## remove weird row in metadata
        metadata = metadata.iloc[1:].copy()
        
        ## add condition
        metadata["is_ctl"] = metadata["donor_id"] == "WT"
        obs_oi = ["donor_id","biosample_id","annotation","is_ctl"] # "treatment","replicate","cell_type","is_ctl"
        sep = "___"
        metadata["condition"] = metadata[obs_oi].apply(lambda x: sep.join(x.astype("str")), axis=1)
        
        # create adata
        adata = sc.AnnData(genexpr.T)
        adata.obs = pd.concat([adata.obs, metadata.set_index("NAME")], axis=1, join="inner")
        
        # save data
        adata.write(output.adata)
        metadata.to_csv(output.metadata, **SAVE_PARAMS)
        
        print("Done!")
        

rule summarize_genexpr:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","{dataset}.h5ad"),
    output:
        adata = os.path.join(PREP_DIR,"singlecell","{dataset}-pseudobulk.h5ad"),
        conditions = os.path.join(PREP_DIR,"singlecell","{dataset}-conditions.tsv.gz")
    resources:
        memory = 150, # GB
        runtime = 3600*2 # h
    params:
        chunk_size = 5000,
    run:
        import scanpy as sc
        import pandas as pd
        import numpy as np
        from tqdm import tqdm
        import gc
        
        # load
        adata = sc.read_h5ad(input.adata)
        chunk_size = params.chunk_size
        
        # summarize gene expression by some grouping variable(s)
        def summarize_genexpr_perts(adata, chunk_size, sep="___"):
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
                adata_chunk = adata_chunk.groupby(["condition"]).mean() # average genexpr by condition
                
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
            obs_oi = ["treatment","replicate","cell_type","is_ctl"]
            conditions[obs_oi] = conditions.reset_index()["condition"].astype(str).str.split(sep, expand=True).values
            conditions["is_ctl"] = conditions["is_ctl"].map({'True': True, 'False': False})
            
            # prepare outputs
            summarized_genexpr = summarized_genexpr.loc[conditions.index]
            
            return summarized_genexpr, conditions
        
        genexpr, conditions = summarize_genexpr_perts(adata, chunk_size)
        
        # make adata
        adata = sc.AnnData(X=genexpr, obs=conditions)
        
        # save
        adata.write(output.adata)
        conditions.reset_index().to_csv(output.conditions, **SAVE_PARAMS)
        
        print("Done!")
        
rule h5ad_to_tsv:
    input:
        adata = os.path.join(PREP_DIR,"singlecell","{dataset}-pseudobulk.h5ad")
    output:
        adata = os.path.join(PREP_DIR,"singlecell","{dataset}-pseudobulk.tsv.gz")
    run:
        import scanpy as sc
        
        adata = sc.read_h5ad(input.adata)
        
        adata = adata.to_df().T
        
        adata.reset_index().to_csv(output.adata, **SAVE_PARAMS)
        
        print("Done!")
    
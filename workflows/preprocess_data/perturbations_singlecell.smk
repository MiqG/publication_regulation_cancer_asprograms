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
        # preprocess ReplogleWeissman2022
        expand(os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad"), dataset=DATASETS)
        
        
rule preprocess_scperturb:
    input:
        adata = os.path.join(DATA_DIR,"scPerturb","ReplogleWeissman2022_{dataset}.h5ad"),
    output:
        adata = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_{dataset}.h5ad")
    run:
        import scanpy as sc
        
        # load
        adata = sc.read_h5ad(input.adata)
        
        # prep
        
        ## metadata
        ### column with gene name and ensembl of perturbed gene
        
        ## gene expression
        ### filter out cells
        sc.pp.filter_cells(adata, min_genes=200)
        ### filter out genes
        sc.pp.filter_genes(adata, min_cells=3)
        ### make QC metrics
        adata.var["mt"] = adata.var["gene_name"].str.startswith("MT-").fillna(False)
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        ### normalize counts
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        ### mark highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
        ## save as raw attribute up to here
        adata.raw = adata
        
        # save
        adata.write(output.adata)
        
        print("Done!")
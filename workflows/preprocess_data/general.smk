import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

OMIC_TYPES = ["EX","genexpr_tpm"]

##### RULES #####
rule all:
    input:
        # preprocess STRINGDB
        os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz'),
        
        # prep Urbanski2022
        os.path.join(PREP_DIR,"metadata","Urbanski2022.tsv.gz"),
        os.path.join(PREP_DIR,'event_psi','Urbanski2022-EX.tsv.gz'),
        os.path.join(PREP_DIR,'event_psi','Urbanski2022-ALTA.tsv.gz'),
        os.path.join(PREP_DIR,'event_psi','Urbanski2022-ALTD.tsv.gz'),
        os.path.join(PREP_DIR,'event_psi','Urbanski2022-INT.tsv.gz'),
        os.path.join(PREP_DIR,'genexpr_tpm','Urbanski2022.tsv.gz'),
        
        # signatures Urbanski2022
        os.path.join(PREP_DIR,"signatures","Urbanski2022-EX.tsv.gz"),
        os.path.join(PREP_DIR,"signatures","Urbanski2022-genexpr_tpm.tsv.gz"),
        
        # prep Rogalska2024
        os.path.join(PREP_DIR,"metadata","Rogalska2024.tsv.gz"),
        os.path.join(PREP_DIR,'event_psi','Rogalska2024-EX.tsv.gz'),
        os.path.join(PREP_DIR,'event_psi','Rogalska2024-ALTA.tsv.gz'),
        os.path.join(PREP_DIR,'event_psi','Rogalska2024-ALTD.tsv.gz'),
        os.path.join(PREP_DIR,'event_psi','Rogalska2024-INT.tsv.gz'),
        os.path.join(PREP_DIR,'genexpr_tpm','Rogalska2024.tsv.gz'),
        
        # signatures Rogalska2024
        os.path.join(PREP_DIR,"signatures","Rogalska2024-EX.tsv.gz"),
        os.path.join(PREP_DIR,"signatures","Rogalska2024-genexpr_tpm.tsv.gz"), 
        
        
rule preprocess_stringdb:
    input:
        ppi = os.path.join(RAW_DIR,'STRINGDB','9606.protein.links.full.v11.5.txt.gz'),
        aliases = os.path.join(RAW_DIR,'STRINGDB','9606.protein.aliases.v11.5.txt.gz')
    output:
        os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz')
    shell:
        """
        python scripts/preprocess_stringdb.py \
                    --raw_ppi_file={input.ppi} \
                    --raw_aliases_file={input.aliases} \
                    --prep_ppi_file={output}
        """
        
        
rule preprocess_Urbanski2022:
    input:
        metadata = os.path.join(SUPPORT_DIR,"ENA_filereport-PRJNA754112-Urbanski2022.tsv"),
        psi = os.path.join(RAW_DIR,"articles","Urbanski2022",'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        genexpr = os.path.join(RAW_DIR,"articles","Urbanski2022",'vast_out','TPM-hg38-24.tab.gz')
    output:
        metadata = os.path.join(PREP_DIR,"metadata","Urbanski2022.tsv.gz"),
        psi_EX = os.path.join(PREP_DIR,'event_psi','Urbanski2022-EX.tsv.gz'),
        psi_ALTA = os.path.join(PREP_DIR,'event_psi','Urbanski2022-ALTA.tsv.gz'),
        psi_ALTD = os.path.join(PREP_DIR,'event_psi','Urbanski2022-ALTD.tsv.gz'),
        psi_INT = os.path.join(PREP_DIR,'event_psi','Urbanski2022-INT.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','Urbanski2022.tsv.gz')      
    run:
        import gc
        import pandas as pd
        import numpy as np
        
        # load
        print("Loading data...")
        metadata = pd.read_table(input.metadata)
        psi = pd.read_table(input.psi, index_col=0)
        genexpr = pd.read_table(input.genexpr, index_col=[0,1])
        
        gc.collect()
        
        # metadata
        metadata["sampleID"] = metadata["run_accession"]
        
        # PSI
        print("Processing PSI matrix...")
        ## drop empty rows
        is_na = psi.isnull()
        non_missing = is_na.shape[1] - is_na.sum(1)
        to_keep = non_missing >= 1
        psi = psi.loc[to_keep]
        
        ## remove vast-tools' suffix
        psi.columns = [c.replace('_1','') for c in psi.columns]
        
        ## split by event type
        event_types = ["EX","ALTA","ALTD","INT"]
        psis = {e: psi.loc[psi.index.str.contains(e)] for e in event_types}
        
        # TPM
        print("Processing TPM matrix...")
        ## remove vast-tools' suffix
        genexpr.columns = [c.replace('_1','') for c in genexpr.columns]
        
        ## log-transform
        genexpr = np.log2(genexpr + 1)
        
        # subset
        ## find common samples
        common_samples = list(set(metadata["run_accession"]).intersection(
            psis["EX"].columns
        ).intersection(
            genexpr.columns
        ))
        psis = {e: psis[e][common_samples].copy() for e in event_types}
        genexpr = genexpr[common_samples].copy()
        metadata = metadata.loc[metadata["run_accession"].isin(common_samples)]
        
        # save
        print("Saving...")
        ## metadata
        metadata.to_csv(output.metadata, **SAVE_PARAMS)

        ## PSIs
        psis["EX"].reset_index().to_csv(output.psi_EX, **SAVE_PARAMS)
        psis["ALTD"].reset_index().to_csv(output.psi_ALTD, **SAVE_PARAMS)
        psis["ALTA"].reset_index().to_csv(output.psi_ALTA, **SAVE_PARAMS)
        psis["INT"].reset_index().to_csv(output.psi_INT, **SAVE_PARAMS)
        
        ## TPMs
        genexpr.reset_index().drop(columns='NAME').to_csv(output.genexpr, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule preprocess_Rogalska2024:
    input:
        metadata = os.path.join(SUPPORT_DIR,'ENA_filereport-PRJEB49033-Rogalska2024.tsv'),
        psi = os.path.join(RAW_DIR,"articles","Rogalska2024",'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        genexpr = os.path.join(RAW_DIR,"articles","Rogalska2024",'vast_out','TPM-hg38-319.tab.gz'),
        annot = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        metadata = os.path.join(PREP_DIR,"metadata","Rogalska2024.tsv.gz"),
        psi_EX = os.path.join(PREP_DIR,'event_psi','Rogalska2024-EX.tsv.gz'),
        psi_ALTA = os.path.join(PREP_DIR,'event_psi','Rogalska2024-ALTA.tsv.gz'),
        psi_ALTD = os.path.join(PREP_DIR,'event_psi','Rogalska2024-ALTD.tsv.gz'),
        psi_INT = os.path.join(PREP_DIR,'event_psi','Rogalska2024-INT.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','Rogalska2024.tsv.gz')      
    run:
        import gc
        import pandas as pd
        import numpy as np
        
        # load
        print("Loading data...")
        metadata = pd.read_table(input.metadata)
        annot = pd.read_table(input.annot)
        psi = pd.read_table(input.psi, index_col=0)
        genexpr = pd.read_table(input.genexpr, index_col=[0,1])
        
        gc.collect()
        
        # metadata
        metadata["sampleID"] = metadata["run_accession"]
        idx_ctls = metadata["sample_title"].str.startswith("AA")
        metadata["control_samples"] = "||".join(metadata.loc[idx_ctls,"sampleID"])
        metadata["cell_line"] = "HELA_CERVIX"
        metadata["PERT_GENE"] = metadata["sample_title"].str.split("_", expand=True)[0]
        metadata["PERT_ENSEMBL"] = pd.merge(metadata[["PERT_GENE"]], annot[["Approved symbol","Ensembl gene ID"]], left_on="PERT_GENE", right_on="Approved symbol", how="left")["Ensembl gene ID"]
        metadata["PERT_TYPE"] = "KNOCKDOWN"
        metadata.loc[idx_ctls,"PERT_TYPE"] = np.nan
        
        # PSI
        print("Processing PSI matrix...")
        ## drop empty rows
        is_na = psi.isnull()
        non_missing = is_na.shape[1] - is_na.sum(1)
        to_keep = non_missing >= 1
        psi = psi.loc[to_keep]
        
        ## remove vast-tools' suffix
        psi.columns = [c.replace('_1','') for c in psi.columns]
        
        ## split by event type
        event_types = ["EX","ALTA","ALTD","INT"]
        psis = {e: psi.loc[psi.index.str.contains(e)] for e in event_types}
        
        # TPM
        print("Processing TPM matrix...")
        ## remove vast-tools' suffix
        genexpr.columns = [c.replace('_1','') for c in genexpr.columns]
        
        ## log-transform
        genexpr = np.log2(genexpr + 1)
        
        # subset
        ## find common samples
        common_samples = list(set(metadata["run_accession"]).intersection(
            psis["EX"].columns
        ).intersection(
            genexpr.columns
        ))
        psis = {e: psis[e][common_samples].copy() for e in event_types}
        genexpr = genexpr[common_samples].copy()
        metadata = metadata.loc[metadata["run_accession"].isin(common_samples)]
        
        # save
        print("Saving...")
        ## metadata
        metadata.to_csv(output.metadata, **SAVE_PARAMS)

        ## PSIs
        psis["EX"].reset_index().to_csv(output.psi_EX, **SAVE_PARAMS)
        psis["ALTD"].reset_index().to_csv(output.psi_ALTD, **SAVE_PARAMS)
        psis["ALTA"].reset_index().to_csv(output.psi_ALTA, **SAVE_PARAMS)
        psis["INT"].reset_index().to_csv(output.psi_INT, **SAVE_PARAMS)
        
        ## TPMs
        genexpr.reset_index().drop(columns='NAME').to_csv(output.genexpr, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_signatures_splicing:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','{dataset}.tsv.gz'),
        splicing = os.path.join(PREP_DIR,'event_psi','{dataset}-EX.tsv.gz'),
    output:
        signatures = os.path.join(PREP_DIR,"signatures","{dataset}-EX.tsv.gz")
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # load
        metadata = pd.read_table(input.metadata)
        splicing = pd.read_table(input.splicing, index_col=0)
        
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
        
        
rule compute_signatures_genexpr:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','{dataset}.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','{dataset}.tsv.gz'),
    output:
        signatures = os.path.join(PREP_DIR,"signatures","{dataset}-genexpr_tpm.tsv.gz")
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # load
        metadata = pd.read_table(input.metadata)
        genexpr = pd.read_table(input.genexpr, index_col=0)
        
        # subset
        common_samples = list(set(metadata["sampleID"]).intersection(genexpr.columns))
        metadata = metadata.loc[metadata["sampleID"].isin(common_samples)].copy()
        genexpr = genexpr[common_samples].copy()
        
        # delta PSI as the difference between conditions and the mean of the conditions
        signatures = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                
                # there may be empty controls
                if any(genexpr.columns.isin(ctls)):
                    genexpr_ctls = genexpr.loc[:,genexpr.columns.isin(ctls)].mean(axis=1)

                    # compute log fold change
                    fold_change = genexpr[sample_oi] - genexpr_ctls

                    signatures[sample_oi] = fold_change

                    del fold_change, genexpr_ctls, ctls
                else:
                    print(genexpr.columns, ctls)
                    print(genexpr.columns.isin(ctls))
                    continue

        signatures = pd.DataFrame(signatures)
        
        # save
        signatures.reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
        print("Done!")
        
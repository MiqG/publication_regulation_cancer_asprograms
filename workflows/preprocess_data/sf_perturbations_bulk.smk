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

EVENT_TYPES = ["EX"]
CELL_LINES = ["K562","HepG2"]
ENCORE_DATASETS = ["ENCOREKD","ENCOREKO"]

##### RULES #####
rule all:
    input:
        # preprocess ENCORE
        os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKD-EX.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKD-ALTA.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKD-ALTD.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKD-INT.tsv.gz"),
        os.path.join(PREP_DIR,"genexpr_tpm","ENCOREKD.tsv.gz"),    
        
        # preprocess ENCOREKO
        os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKO-EX.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKO-ALTA.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKO-ALTD.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENCOREKO-INT.tsv.gz"),
        os.path.join(PREP_DIR,"genexpr_tpm","ENCOREKO.tsv.gz"),    

        # preprocess ENASFS
        os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENASFS-EX.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENASFS-ALTA.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENASFS-ALTD.tsv.gz"),
        os.path.join(PREP_DIR,"event_psi","ENASFS-INT.tsv.gz"),
        os.path.join(PREP_DIR,"genexpr_tpm","ENASFS.tsv.gz"),    
        
        # PERTURBATION SCREENS
        ### ENCORE
        ### compute log FC TPM
        expand(os.path.join(PREP_DIR,'pert_transcriptomes','{dataset}',"{cell_line}",'log2_fold_change_tpm.tsv.gz'), cell_line=CELL_LINES, dataset=ENCORE_DATASETS),
        ### compute delta PSI
        expand(os.path.join(PREP_DIR,'pert_transcriptomes','{dataset}',"{cell_line}",'delta_psi-{event_type}.tsv.gz'), event_type=EVENT_TYPES, cell_line=CELL_LINES, dataset=ENCORE_DATASETS),
        ### simplify delta PSI and log FC TPM
        expand(os.path.join(PREP_DIR,'ground_truth_pert','{dataset}',"{cell_line}",'delta_psi-{event_type}.tsv.gz'), event_type=EVENT_TYPES, cell_line=CELL_LINES, dataset=ENCORE_DATASETS),
        expand(os.path.join(PREP_DIR,'ground_truth_pert','{dataset}',"{cell_line}",'log2_fold_change_tpm.tsv.gz'), cell_line=CELL_LINES, dataset=ENCORE_DATASETS),
        
        ## ENASFS
        ### delta PSI
        expand(os.path.join(PREP_DIR,'pert_transcriptomes','ENASFS','delta_psi-{event_type}.tsv.gz'), event_type=EVENT_TYPES),
        ### genexpr fold change
        os.path.join(PREP_DIR,'pert_transcriptomes','ENASFS','log2_fold_change_tpm.tsv.gz'),
        ### simplify delta PSI
        expand(os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{event_type}.tsv.gz'), event_type=EVENT_TYPES),
        ### simplify genexpr fold change
        os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz'),
        
        # EDA data figures
        # os.path.join(RESULTS_DIR,'figures','eda')
        

rule preprocess_encorekd:
    input:
        metadata = os.path.join(RAW_DIR,'ENCODE','ENCORE','metadata','ENCORE.tsv'),
        psi = os.path.join(RAW_DIR,'ENCODE','ENCORE','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        genexpr = os.path.join(RAW_DIR,'ENCODE','ENCORE','vast_out','TPM-hg38-1097.tab.gz')
    output:
        metadata = os.path.join(PREP_DIR,'metadata','ENCOREKD.tsv.gz'),
        psi_EX = os.path.join(PREP_DIR,'event_psi','ENCOREKD-EX.tsv.gz'),
        psi_ALTA = os.path.join(PREP_DIR,'event_psi','ENCOREKD-ALTA.tsv.gz'),
        psi_ALTD = os.path.join(PREP_DIR,'event_psi','ENCOREKD-ALTD.tsv.gz'),
        psi_INT = os.path.join(PREP_DIR,'event_psi','ENCOREKD-INT.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','ENCOREKD.tsv.gz')
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
        print("Processing metadata...")
        ## sample names
        metadata["sampleID"] = metadata["dbxrefs"].str.replace("SRA:","")
        ## cell lines info
        metadata["cell_line"] = metadata["Biosample term name"]
        depmapids = {"K562":"ACH-000551", "HepG2":"ACH-000739"}
        metadata["DepMap_ID"] = [depmapids[c] for c in metadata["cell_line"]]
        ## perturbation info
        metadata["PERT_GENE"] = metadata["Experiment target"].str.replace("-human","")
        gene_annot = genexpr.index.to_frame().rename(columns={"ID":"PERT_ENSEMBL", "NAME":"PERT_GENE"})
        metadata = pd.merge(metadata, gene_annot, how="left", on="PERT_GENE")
        metadata["PERT_TYPE"] = "SHRNAKD"
        ## experiment
        metadata["experiment"] = metadata["Experiment accession"]
        ## replicate
        metadata["replicate"] = metadata["Biological replicate(s)"]
        
        ## controls
        ctls_exps = []
        ctls_samps = []
        for idx, row in metadata.iterrows():
            if isinstance(row["Controlled by"], str):
                # get file accession controls
                accs = row["Controlled by"]\
                        .replace("files","")\
                        .replace("/","")\
                        .replace(" ","")\
                        .split(",")
                idx = metadata["File accession"].isin(accs)

                # get experiment accession
                exps = metadata.loc[idx, "experiment"].unique()
                
                # get sample accession
                samps = metadata.loc[idx, "sampleID"].unique()
                
                # save
                exps = ','.join(np.sort(exps))
                samps = ','.join(np.sort(samps))
                ctls_exps.append(exps)
                ctls_samps.append(samps)
            else:
                ctls_exps.append(np.nan)
                ctls_samps.append(np.nan)
        metadata["control_experiment"] = ctls_exps
        metadata["control_samples"] = ctls_samps
        
        cols_oi = ['sampleID','cell_line', 'DepMap_ID', 'PERT_GENE', 'PERT_ENSEMBL', 'experiment', 
                   'control_experiment', 'control_samples','replicate']
        metadata = metadata[cols_oi].drop_duplicates()
        
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
        
        
rule preprocess_encoreko:
    input:
        metadata = os.path.join(RAW_DIR,'ENCODE','ENCORE','CRISPRKO','metadata','CRISPRKO.tsv'),
        psi = os.path.join(RAW_DIR,'ENCODE','ENCORE','CRISPRKO','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        genexpr = os.path.join(RAW_DIR,'ENCODE','ENCORE','CRISPRKO','vast_out','TPM-hg38-924.tab.gz')
    output:
        metadata = os.path.join(PREP_DIR,'metadata','ENCOREKO.tsv.gz'),
        psi_EX = os.path.join(PREP_DIR,'event_psi','ENCOREKO-EX.tsv.gz'),
        psi_ALTA = os.path.join(PREP_DIR,'event_psi','ENCOREKO-ALTA.tsv.gz'),
        psi_ALTD = os.path.join(PREP_DIR,'event_psi','ENCOREKO-ALTD.tsv.gz'),
        psi_INT = os.path.join(PREP_DIR,'event_psi','ENCOREKO-INT.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','ENCOREKO.tsv.gz')
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
        print("Processing metadata...")
        ## sample names
        metadata["sampleID"] = metadata["dbxrefs"].str.replace("SRA:","")
        ## cell lines info
        metadata["cell_line"] = metadata["Biosample term name"]
        depmapids = {"K562":"ACH-000551", "HepG2":"ACH-000739"}
        metadata["DepMap_ID"] = [depmapids[c] for c in metadata["cell_line"]]
        ## perturbation info
        metadata["PERT_GENE"] = metadata["Experiment target"].str.replace("-human","")
        gene_annot = genexpr.index.to_frame().rename(columns={"ID":"PERT_ENSEMBL", "NAME":"PERT_GENE"})
        metadata = pd.merge(metadata, gene_annot, how="left", on="PERT_GENE")
        metadata["PERT_TYPE"] = "CRISPRKO"
        ## experiment
        metadata["experiment"] = metadata["Experiment accession"]
        ## replicate
        metadata["replicate"] = metadata["Biological replicate(s)"]
        
        ## controls
        ctls_exps = []
        ctls_samps = []
        for idx, row in metadata.iterrows():
            if isinstance(row["Controlled by"], str):
                # get file accession controls
                accs = row["Controlled by"]\
                        .replace("files","")\
                        .replace("/","")\
                        .replace(" ","")\
                        .split(",")
                idx = metadata["File accession"].isin(accs)

                # get experiment accession
                exps = metadata.loc[idx, "experiment"].unique()
                
                # get sample accession
                samps = metadata.loc[idx, "sampleID"].unique()
                
                # save
                exps = ','.join(np.sort(exps))
                samps = ','.join(np.sort(samps))
                ctls_exps.append(exps)
                ctls_samps.append(samps)
            else:
                ctls_exps.append(np.nan)
                ctls_samps.append(np.nan)
        metadata["control_experiment"] = ctls_exps
        metadata["control_samples"] = ctls_samps
        
        cols_oi = ['sampleID','cell_line', 'DepMap_ID', 'PERT_GENE', 'PERT_ENSEMBL', 'experiment', 
                   'control_experiment', 'control_samples','replicate']
        metadata = metadata[cols_oi].drop_duplicates()
        
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
        

rule diff_tpm_encore:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','{dataset}.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','{dataset}.tsv.gz')
    output:
        diff_tpm = os.path.join(PREP_DIR,'pert_transcriptomes','{dataset}',"{cell_line}",'log2_fold_change_tpm.tsv.gz')
    params:
        cell_line = "{cell_line}"
    run:
        import pandas as pd
        import numpy as np
        
        metadata = pd.read_table(input.metadata)
        genexpr = pd.read_table(input.genexpr, index_col=0)
        cell_line = params.cell_line
        
        # subset by cell line
        metadata = metadata.loc[metadata["cell_line"]==cell_line].copy()
        
        # log transform (already done)
        #genexpr = np.log2(genexpr + 1)
        #genexpr.columns = [c.replace("_1","") for c in genexpr.columns]
        
        # as the difference between conditions and the mean of the conditions
        diff_tpm = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split(",")
                tpm_ctls = genexpr[ctls].mean(axis=1)
                
                # compute log2 fold-change
                dtpm = genexpr[sample_oi] - tpm_ctls
                
                diff_tpm[sample_oi] = dtpm
                
                del dtpm, tpm_ctls, ctls
                
        diff_tpm = pd.DataFrame(diff_tpm)
        
        # save
        diff_tpm.reset_index().to_csv(output.diff_tpm, sep="\t", index=False, compression="gzip")
        
        print("Done!")
        
        
rule delta_psi_encore:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','{dataset}.tsv.gz'),
        psi = os.path.join(PREP_DIR,'event_psi','{dataset}-{event_type}.tsv.gz')
    output:
        delta_psi = os.path.join(PREP_DIR,'pert_transcriptomes','{dataset}','{cell_line}','delta_psi-{event_type}.tsv.gz'),
        delta_psi_rel = os.path.join(PREP_DIR,'pert_transcriptomes','{dataset}','{cell_line}','delta_psi_rel-{event_type}.tsv.gz')
    params:
        cell_line = "{cell_line}"
    run:
        import pandas as pd
        import numpy as np
        
        metadata = pd.read_table(input.metadata)
        psi = pd.read_table(input.psi, index_col=0)
        cell_line = params.cell_line
        
        # subset by cell line
        metadata = metadata.loc[metadata["cell_line"]==cell_line].copy()
        
        # delta PSI as the difference between conditions and the mean of the conditions
        delta_psi = {}
        delta_psi_rel = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split(",")
                psi_ctls = psi[ctls].mean(axis=1)
                
                # compute delta PSI
                dpsi = psi[sample_oi] - psi_ctls
                
                # compute relative delta PSI (% changed of possible change w.r.t control)
                extreme_psi = dpsi.copy()
                extreme_psi[dpsi < 0] = 0 - psi_ctls # if it has decreased, the minimum inclusion
                extreme_psi[dpsi > 0] = 100 - psi_ctls # if it has increased, the maximum inclusion
                extreme_psi[dpsi == 0] = 1 # no change (for completeness)
                
                ## what percentage of the maximum possible change in inclusion has occured?
                rel_dpsi = np.sign(dpsi) * (dpsi / extreme_psi) * 100
                
                delta_psi[sample_oi] = dpsi
                delta_psi_rel[sample_oi] = rel_dpsi
                
                del dpsi, rel_dpsi, psi_ctls, ctls
                
        
        delta_psi = pd.DataFrame(delta_psi)
        delta_psi_rel = pd.DataFrame(delta_psi_rel)
        
        # save
        delta_psi.reset_index().to_csv(output.delta_psi, sep="\t", index=False, compression="gzip")
        delta_psi_rel.reset_index().to_csv(output.delta_psi_rel, sep="\t", index=False, compression="gzip")
        
        print("Done!")
        
        
rule prepare_ground_truth_pert_dpsi:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','{dataset}.tsv.gz'),
        dpsi = os.path.join(PREP_DIR,'pert_transcriptomes','{dataset}',"{cell_line}",'delta_psi-{event_type}.tsv.gz'),
    params:
        cell_line = "{cell_line}"
    output:
        dpsi = os.path.join(PREP_DIR,'ground_truth_pert','{dataset}',"{cell_line}",'delta_psi-{event_type}.tsv.gz'),
    run:
        import pandas as pd
        
        metadata = pd.read_table(input.metadata)
        dpsi = pd.read_table(input.dpsi, index_col=0)
        cell_line = params.cell_line
        
        # drop control samples
        metadata = metadata.loc[~metadata["PERT_ENSEMBL"].isnull()].copy()
        
        # subset by cell line
        metadata = metadata.loc[metadata["cell_line"]==cell_line].copy()
        
        dpsis = []
        for ensembl in metadata["PERT_ENSEMBL"].unique():
            samples_oi = metadata.loc[metadata["PERT_ENSEMBL"]==ensembl, "sampleID"]
            
            dpsi_avg = dpsi[samples_oi].mean(axis=1)
            dpsi_avg.name = ensembl
            dpsis.append(dpsi_avg)
            
        dpsis = pd.concat(dpsis, axis=1)

        dpsis.reset_index().to_csv(output.dpsi, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule prepare_ground_truth_pert_logFCtpm:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','{dataset}.tsv.gz'),
        diff_tpm = os.path.join(PREP_DIR,'pert_transcriptomes','{dataset}',"{cell_line}",'log2_fold_change_tpm.tsv.gz')
    params:
        cell_line = "{cell_line}"
    output:
        diff_tpm = os.path.join(PREP_DIR,'ground_truth_pert','{dataset}',"{cell_line}",'log2_fold_change_tpm.tsv.gz')
    run:
        import pandas as pd
        
        metadata = pd.read_table(input.metadata)
        diff_tpm = pd.read_table(input.diff_tpm, index_col=0)
        cell_line = params.cell_line
        
        # drop control samples
        metadata = metadata.loc[~metadata["PERT_ENSEMBL"].isnull()].copy()
        
        # subset by cell line
        metadata = metadata.loc[metadata["cell_line"]==cell_line].copy()
        
        diff_tpms = []
        for ensembl in metadata["PERT_ENSEMBL"].unique():
            samples_oi = metadata.loc[metadata["PERT_ENSEMBL"]==ensembl, "sampleID"]
            
            diff_tpm_avg = diff_tpm[samples_oi].mean(axis=1)
            diff_tpm_avg.name = ensembl
            diff_tpms.append(diff_tpm_avg)
            
        diff_tpms = pd.concat(diff_tpms, axis=1)

        diff_tpms.reset_index().to_csv(output.diff_tpm, **SAVE_PARAMS)
        
        print("Done!")
    
    
rule preprocess_ena_splicing_factors:
    input:
        metadata = os.path.join(SUPPORT_DIR,'ENA_filereport-selected_sf_experiments_handcurated_w_control_samples.tsv'),
        metadata_groups = os.path.join(SUPPORT_DIR,"ENA_filereport-selected_sf_experiments_handcurated-low_read_count_groups.tsv"),
        psi = os.path.join(RAW_DIR,"ENA","splicing_factors",'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        genexpr = os.path.join(RAW_DIR,"ENA","splicing_factors",'vast_out','TPM-hg38-1597.tab.gz')
    output:
        metadata = os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"),
        psi_EX = os.path.join(PREP_DIR,'event_psi','ENASFS-EX.tsv.gz'),
        psi_ALTA = os.path.join(PREP_DIR,'event_psi','ENASFS-ALTA.tsv.gz'),
        psi_ALTD = os.path.join(PREP_DIR,'event_psi','ENASFS-ALTD.tsv.gz'),
        psi_INT = os.path.join(PREP_DIR,'event_psi','ENASFS-INT.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','ENASFS.tsv.gz')      
    run:
        import gc
        import pandas as pd
        import numpy as np
        
        # load
        print("Loading data...")
        metadata = pd.read_table(input.metadata)
        metadata_groups = pd.read_table(input.metadata_groups)
        psi = pd.read_table(input.psi, index_col=0)
        genexpr = pd.read_table(input.genexpr, index_col=[0,1])
        
        gc.collect()
        
        # metadata
        ## not grouped metadata
        metadata_not_merged = metadata.loc[~metadata["run_accession"].isin(metadata_groups["run_accession"])].copy()
        metadata_not_merged["sampleID"] = metadata_not_merged["run_accession"]
        metadata_not_merged["control_samples"] = metadata_not_merged["control_samples"].str.replace(",","||")
        
        ## grouped metadata
        ### subset
        metadata_merged = metadata.loc[metadata["run_accession"].isin(metadata_groups["run_accession"])].copy()
        metadata_merged = pd.merge(metadata_merged, metadata_groups[["group_label","run_accession"]], how="left", on="run_accession")
        ### control_samples at the group level
        control_samples_group = []
        for ctls in metadata_merged["control_samples"].dropna().unique():
            # get run accessions
            run_accessions = ctls.split(",")
            
            # get corresponding group labels
            group_labs = metadata_merged.loc[metadata_merged["run_accession"].isin(run_accessions),"group_label"].unique()
            
            # join group labels
            ctls_group = "||".join(group_labs)
            
            # save
            control_samples_group.append({
                "control_samples": ctls,
                "control_samples_group": ctls_group
            })
        control_samples_group = pd.DataFrame(control_samples_group)
        metadata_merged = pd.merge(metadata_merged, control_samples_group, on="control_samples", how="left")
        metadata_merged["control_samples"] = metadata_merged["control_samples_group"]
        
        ### columns to keep
        cols_oi = [
            "study_accession","cell_line_name","condition","pert_time",
            "pert_time_units","pert_concentration","pert_concentration_units",
            'found_sfs_in_run_alias', 'found_sfs_in_sample_alias', 'found_sfs_in_sample_title', 
            'found_sfs_in_experiment_title', 'found_sfs_in_study_title', 'is_match', 'group_label',
            "DepMap_ID","CCLE_Name","PERT_GENE","PERT_ENSEMBL","IS_USEFUL","comments","PERT_TYPE",
            "control_samples"
        ]
        metadata_merged = metadata_merged[cols_oi].drop_duplicates().copy()
        metadata_merged["sampleID"] = metadata_merged["group_label"]
        
        ## combine metadatas
        metadata = pd.concat([metadata_not_merged, metadata_merged])
        
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
        common_samples = set(metadata["sampleID"]).intersection(
            psis["EX"].columns
        ).intersection(
            genexpr.columns
        )
        psis = {e: psis[e][common_samples].copy() for e in event_types}
        genexpr = genexpr[common_samples].copy()
        metadata = metadata.loc[metadata["sampleID"].isin(common_samples)]
        
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
        
        
rule delta_psi_ena_splicing_factors:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        psi = os.path.join(PREP_DIR,'event_psi','ENASFS-{event_type}.tsv.gz')
    output:
        delta_psi = os.path.join(PREP_DIR,'pert_transcriptomes','ENASFS','delta_psi-{event_type}.tsv.gz')
    run:
        import pandas as pd
        
        # load
        metadata = pd.read_table(input.metadata)
        psi = pd.read_table(input.psi, index_col=0)
        
        # delta PSI as the difference between conditions and the mean of the conditions
        delta_psi = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                psi_ctls = psi[ctls].mean(axis=1)
                
                # compute delta PSI
                dpsi = psi[sample_oi] - psi_ctls
                
                delta_psi[sample_oi] = dpsi
                
                del dpsi, psi_ctls, ctls
                
        
        delta_psi = pd.DataFrame(delta_psi)
        
        # save
        delta_psi.reset_index().to_csv(output.delta_psi, **SAVE_PARAMS)
        
        print("Done!")

        
rule diff_tpm_ena_splicing_factors:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','ENASFS.tsv.gz')
    output:
        diff_tpm = os.path.join(PREP_DIR,'pert_transcriptomes','ENASFS','log2_fold_change_tpm.tsv.gz')
    run:
        import pandas as pd
        import numpy as np
        
        metadata = pd.read_table(input.metadata)
        genexpr = pd.read_table(input.genexpr, index_col=0)
        
        # as the difference between conditions and the mean of the conditions
        diff_tpm = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                tpm_ctls = genexpr[ctls].mean(axis=1)
                
                # compute log2 fold-change
                dtpm = genexpr[sample_oi] - tpm_ctls
                
                diff_tpm[sample_oi] = dtpm
                
                del dtpm, tpm_ctls, ctls
                
        diff_tpm = pd.DataFrame(diff_tpm)
        
        # save
        diff_tpm.reset_index().to_csv(output.diff_tpm, sep="\t", index=False, compression="gzip")
        
        print("Done!")
        
        
rule prepare_ground_truth_pert_dpsi_ena_splicing_factors:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        dpsi = os.path.join(PREP_DIR,'pert_transcriptomes','ENASFS','delta_psi-{event_type}.tsv.gz')
    output:
        dpsi = os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{event_type}.tsv.gz')
    run:
        import pandas as pd
        
        metadata = pd.read_table(input.metadata)
        dpsi = pd.read_table(input.dpsi, index_col=0)
        
        # drop control samples
        metadata = metadata.loc[~metadata["PERT_ENSEMBL"].isnull()].copy()
        metadata["PERT_ID"] = metadata[
            ["study_accession","cell_line_name","PERT_ENSEMBL"]
        ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
        
        dpsis = []
        for pert_id in metadata["PERT_ID"].unique():
            samples_oi = metadata.loc[metadata["PERT_ID"]==pert_id, "sampleID"]

            dpsi_avg = dpsi[samples_oi].mean(axis=1)
            dpsi_avg.name = pert_id
            dpsis.append(dpsi_avg)

        dpsis = pd.concat(dpsis, axis=1)

        dpsis.reset_index().to_csv(output.dpsi, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule prepare_ground_truth_pert_logFCtpm_ena_splicing_factors:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        diff_tpm = os.path.join(PREP_DIR,'pert_transcriptomes','ENASFS','log2_fold_change_tpm.tsv.gz')
    output:
        diff_tpm = os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz')
    run:
        import pandas as pd
        
        metadata = pd.read_table(input.metadata)
        diff_tpm = pd.read_table(input.diff_tpm, index_col=0)
        
        # drop control samples
        metadata = metadata.loc[~metadata["PERT_ENSEMBL"].isnull()].copy()
        metadata["PERT_ID"] = metadata[
            ["study_accession","cell_line_name","PERT_ENSEMBL"]
        ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
        
        diff_tpms = []
        for pert_id in metadata["PERT_ID"].unique():
            samples_oi = metadata.loc[metadata["PERT_ID"]==pert_id, "sampleID"]
            
            diff_tpm_avg = diff_tpm[samples_oi].mean(axis=1)
            diff_tpm_avg.name = pert_id
            diff_tpms.append(diff_tpm_avg)
            
        diff_tpms = pd.concat(diff_tpms, axis=1)

        diff_tpms.reset_index().to_csv(output.diff_tpm, **SAVE_PARAMS)
        
        print("Done!")
            
            
rule figures_eda:
    input:
        splicing_factors = os.path.join(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv"),
        metadata_encore_kd = os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
        metadata_encore_ko = os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
        kd_screen = os.path.join(SUPPORT_DIR,"kd_screen-symbol.txt"),
        ena_sfs = os.path.join(SUPPORT_DIR,"ENA_filereport-selected_sf_experiments_handcurated.tsv")
    output:
        directory(os.path.join(RESULTS_DIR,'figures','eda'))
    shell:
        """
        Rscript scripts/figures_eda.R \
                    --splicing_factors_file={input.splicing_factors} \
                    --metadata_encore_kd_file={input.metadata_encore_kd} \
                    --metadata_encore_ko_file={input.metadata_encore_ko} \
                    --kd_screen_file={input.kd_screen} \
                    --ena_sfs_file={input.ena_sfs} \
                    --figs_dir={output}
        """
        
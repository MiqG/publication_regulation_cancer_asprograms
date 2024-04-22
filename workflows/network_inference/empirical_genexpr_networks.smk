import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","network_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

PERT_SPLICING_FILES = [
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-{omic_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-{omic_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{omic_type}.tsv.gz')
]

PERT_GENEXPR_FILES = [
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'log2_fold_change_tpm.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'log2_fold_change_tpm.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz')
]

PERT_SCGENEXPR_FILES = [
    os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    #os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-pseudobulk_by_batch-log2_fold_change_cpm.tsv.gz"),
    #os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_by_batch-log2_fold_change_cpm.tsv.gz"),
    #os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
]

PERT_FILES = {
    "EX": PERT_SPLICING_FILES,
    "genexpr": PERT_GENEXPR_FILES,
    "scgenexpr": PERT_SCGENEXPR_FILES
}


OMIC_TYPES = ["genexpr","scgenexpr"]

##### RULES #####
rule all:
    input:
        # make regulons
        expand(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{omic_type}"), omic_type=OMIC_TYPES),
        
        # prune regulons
        expand(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}"), omic_type=OMIC_TYPES)
        
        
rule make_regulons:
    input:
        perts = lambda wildcards: PERT_FILES[wildcards.omic_type],
        metadata = os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"), # only for ENASFS
        regulators = os.path.join(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{omic_type}"))
    params:
        omic_type = "{omic_type}"
    run:
        import os
        import pandas as pd
        import numpy as np
        
        regulators = pd.read_table(input.regulators)
        omic_type = params.omic_type
        feature_name = "EVENT" if omic_type not in ["genexpr","scgenexpr"] else "ENSEMBL_ID"
        value_name = "delta_psi" if omic_type not in ["genexpr","scgenexpr"] else "log2fc_genexpr"
        
        # prep regulators
        regulators = regulators[["GENE","ENSEMBL"]]

        for f in input.perts:
            print("Loading %s..." % f)
            # load
            perts = pd.read_table(f, index_col=0)
            
            # get info
            if "ENCORE" in f:
                dataset = os.path.basename(os.path.dirname(os.path.dirname(f)))
                cell_line = os.path.basename(os.path.dirname(f))
            
            elif "Replogle" in f:
                dataset = os.path.basename(f).split("-")[0].split("_")[0]
                cell_line = os.path.basename(f).split("-")[0].replace("ReplogleWeissman2022_","")
            
            elif "ENASFS" in f:
                dataset = "ENASFS"
                metadata = pd.read_table(input.metadata)
                
            # prep perturbations
            perts.index.name = feature_name
            perts = perts.melt(
                ignore_index=False, var_name="PERT_ID", value_name=value_name
            ).dropna().reset_index().copy()

            # format
            perts["regulator"] = perts["PERT_ID"]
            perts["target"] = perts[feature_name]
            perts["likelihood"] = np.abs(perts[value_name])
            perts["tfmode"] = (-1)*np.sign(perts[value_name]) # they come from KD or KO, decrease activity
            
            os.makedirs(output.output_dir, exist_ok=True)
            if ("ENCORE" in dataset) | ("Replogle" in dataset):
                # subset
                perts = perts.loc[perts["PERT_ID"].isin(regulators["ENSEMBL"])].copy()

                # add gene symbols
                perts = pd.merge(perts, regulators, left_on="PERT_ID", right_on="ENSEMBL", how="left")

                # save
                output_file = os.path.join(output.output_dir,"%s-%s-%s.tsv.gz") % (dataset, cell_line, value_name)
                print("Saving %s..." % output_file)
                perts.to_csv(output_file, **SAVE_PARAMS)
            
            elif dataset=="ENASFS":
                # correct PERT_ID column
                X = perts["PERT_ID"].str.split("___", expand=True)
                X.columns = ["study_accession","cell_line_name","PERT_ENSEMBL"]
                perts[["study_accession","cell_line_name","PERT_ENSEMBL"]] = X
                perts["regulator"] = perts["PERT_ENSEMBL"]
                
                # subset
                perts = perts.loc[perts["PERT_ENSEMBL"].isin(regulators["ENSEMBL"])].copy()
                
                # add gene symbols
                perts = pd.merge(perts, regulators, left_on="PERT_ENSEMBL", right_on="ENSEMBL", how="left")
                
                # add pert type
                metadata = metadata.loc[~metadata["PERT_ENSEMBL"].isnull()].copy()
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line_name","PERT_ENSEMBL"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                perts = pd.merge(perts, metadata[["PERT_ID","PERT_TYPE"]].drop_duplicates(), on="PERT_ID", how="left")
                
                # subset pert types
                pert_types_oi = ["KNOCKDOWN","KNOCKOUT","OVEREXPRESSION"]
                perts = perts.loc[perts["PERT_TYPE"].isin(pert_types_oi)].copy()
                
                # correct overexpression sign as increase in activity
                idx = perts["PERT_TYPE"]=="OVEREXPRESSION"
                perts.loc[idx,"tfmode"] = -perts.loc[idx,"tfmode"]
                
                # create metaexperiments with a perturbation in each splicing factor
                # try to put perturbations from the same project together
                cols_oi = ["PERT_ENSEMBL","cell_line_name","study_accession"]
                gene_study = perts[cols_oi].drop_duplicates().groupby(
                    cols_oi
                ).size().reset_index().rename(columns={0:"n"}).sort_values(["n","cell_line_name"])
                gene_study["index"] = np.arange(0,len(gene_study))
                
                metaexperiments = {}
                it = 0
                while len(gene_study)>0:
                    to_keep = gene_study["index"].isin(gene_study.groupby('PERT_ENSEMBL')["index"].min())
                    metaexperiments["metaexperiment%s" % it] = gene_study.loc[to_keep].sort_values("PERT_ENSEMBL")
                    gene_study = gene_study.loc[~to_keep].sort_values("PERT_ENSEMBL").copy()
                    it = it + 1
                
                # save
                for metaexperiment_oi in metaexperiments.keys():
                    metaexperiment = metaexperiments[metaexperiment_oi][cols_oi]
                    perts_oi = pd.merge(metaexperiment, perts, on=cols_oi, how="left")

                    # save
                    if len(metaexperiment)>1:
                        output_file = os.path.join(output.output_dir,"%s-%s-%s.tsv.gz") % (dataset, metaexperiment_oi, value_name)
                        print("Saving %s..." % output_file)
                        perts_oi.to_csv(output_file, **SAVE_PARAMS)

            
        print("Done!")

        
rule prune_regulons:
    input:
        regulons_dir = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}"))
    params:
        thresh = lambda wildcards: 15 if wildcards.omic_type not in ["genexpr","scgenexpr"] else 1
    run:
        import os
        import pandas as pd
        
        regulon_files = [os.path.join(input.regulons_dir,f) for f in os.listdir(input.regulons_dir) if f.endswith(".tsv.gz")]
        thresh = params.thresh
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in regulon_files:
            print(regulon_file)
            
            regulon = pd.read_table(regulon_file)
            
            regulon = regulon.loc[regulon["likelihood"]>=thresh].copy()
            
            output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
            regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")
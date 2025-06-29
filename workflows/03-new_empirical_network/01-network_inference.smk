import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","new_empirical_network")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

OMIC_TYPES = ["EX"]

##### RULES #####
rule all:
    input:
        # make regulons
        os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-EX"),
        
        # prune regulons
        os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX"),
        
        # copy viper networks
        os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned_w_viper_networks-EX"),
        os.path.join(RESULTS_DIR,"files","viper_networks-EX"),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","empirical_psi_networks")
        
rule make_regulons:
    input:
        perts = [os.path.join(PREP_DIR,'delta_psi','Rogalska2024-EX.tsv.gz')],
        regulators = os.path.join(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-EX"))
    params:
        omic_type = "EX"
    run:
        import os
        import pandas as pd
        import numpy as np
        
        regulators = pd.read_table(input.regulators)
        omic_type = params.omic_type
        feature_name = "EVENT"
        value_name = "delta_psi"
        
        # prep regulators
        regulators = regulators[["GENE","ENSEMBL"]]

        for f in input.perts:
            print("Loading %s..." % f)
            # load
            perts = pd.read_table(f, index_col=0)
            
            # get info
            dataset = "Rogalska2024"
            cell_line = "HELA_CERVIX"
                
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

            # correct PERT_ID column
            X = perts["PERT_ID"].str.split("___", expand=True)
            X.columns = ["study_accession","cell_line_name","PERT_ENSEMBL","PERT_TYPE"]
            perts[["study_accession","cell_line_name","PERT_ENSEMBL","PERT_TYPE"]] = X
            perts["regulator"] = perts["PERT_ENSEMBL"]

            # subset
            perts = perts.loc[perts["PERT_ENSEMBL"].isin(regulators["ENSEMBL"])].copy()

            # add gene symbols
            perts = pd.merge(perts, regulators, left_on="PERT_ENSEMBL", right_on="ENSEMBL", how="left")

            # save
            output_file = os.path.join(output.output_dir,"%s-%s-%s.tsv.gz") % (dataset, cell_line, value_name)
            print("Saving %s..." % output_file)
            perts.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")

        
rule prune_regulons:
    input:
        regulons_dir = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-EX")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX"))
    params:
        thresh = 15
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
        
        
rule copy_networks_from_viper_splicing_publication:
    input:
        regulons_dir = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX"),
        viper_splicing_networks = [
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment0-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment1-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment2-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment3-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-HepG2-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-K562-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-HepG2-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-K562-delta_psi.tsv.gz")
        ]
    output:
        outdir = directory(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned_w_viper_networks-EX"))
    run:
        import os
        import subprocess
        
        regulons_dir = input.regulons_dir
        
        os.makedirs(output.outdir, exist_ok=True)
        
        # copy regulons
        for f in os.listdir(regulons_dir):
            if f.endswith(".tsv.gz"):
                filename = os.path.join(regulons_dir,f)
                outfile = os.path.join(output.outdir,f)
                cmd = ["cp",filename,outfile]
                print(cmd)
                subprocess.call(cmd)            
            
        # copy viper networks
        for f in input.viper_splicing_networks:
            outfile = os.path.join(output.outdir,os.path.basename(f))
            cmd = ["cp",f,outfile]
            print(cmd)
            subprocess.call(cmd)
        
        print("Done!")

        
rule copy_networks_from_viper_alone:
    input:
        viper_splicing_networks = [
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment0-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment1-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment2-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment3-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-HepG2-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-K562-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-HepG2-delta_psi.tsv.gz"),
            os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-K562-delta_psi.tsv.gz")
        ]
    output:
        outdir = directory(os.path.join(RESULTS_DIR,"files","viper_networks-EX"))
    run:
        import os
        import subprocess
        
        os.makedirs(output.outdir, exist_ok=True)
        
        # copy viper networks
        for f in input.viper_splicing_networks:
            outfile = os.path.join(output.outdir,os.path.basename(f))
            cmd = ["cp",f,outfile]
            print(cmd)
            subprocess.call(cmd)
        
        print("Done!")

        
rule figures_empirical_psi_networks:
    input:
        splicing_factors = os.path.join(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv"),
        networks_viper_alone_dir = os.path.join(RESULTS_DIR,"files","viper_networks-EX"),
        networks_new_w_viper_dir = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned_w_viper_networks-EX"),
        encore_metadata = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","benchmark","metadata_encorekd.tsv.gz"),
        encore_ex = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_ex_encorekd.tsv.gz"),
        encore_alta = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_alta_encorekd.tsv.gz"),
        encore_altd = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_altd_encorekd.tsv.gz"),
        encore_int = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_int_encorekd.tsv.gz"),
    output:
        directory(os.path.join(RESULTS_DIR,"figures","empirical_psi_networks"))
    shell:
        """
        Rscript scripts/figures_empirical_psi_networks.R \
                        --splicing_factors_file={input.splicing_factors} \
                        --networks_viper_alone_dir={input.networks_viper_alone_dir} \
                        --encore_metadata_file={input.encore_metadata} \
                        --encore_ex_file={input.encore_ex} \
                        --encore_alta_file={input.encore_alta} \
                        --encore_altd_file={input.encore_altd} \
                        --encore_int_file={input.encore_int} \
                        --figs_dir={output}
        """
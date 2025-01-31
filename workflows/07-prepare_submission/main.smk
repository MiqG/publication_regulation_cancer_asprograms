import os
import pandas as pd

SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### VARIABLES #####
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
SUPPORT_DIR = os.path.join(ROOT,'support')
RESULTS_DIR = os.path.join(ROOT,'results')
DATA_DIR = os.path.join(ROOT,'data')
RAW_DIR = os.path.join(DATA_DIR,'raw')
PREP_DIR = os.path.join(DATA_DIR,'prep')

PROGRAM_DIR = os.path.join(RESULTS_DIR,"new_empirical_network")

OMIC_GENEXPR_REGULONS = ["bulkgenexpr","scgenexpr","bulkscgenexpr"]
MODEL_ARCHS = ["fclayer","ewlayer"]
K_CROSS_VALIDATION = 5
K_CROSS_VALIDATION = list(range(K_CROSS_VALIDATION))

##### RULES #####
rule all:
    input:
        # supplementary tables
        os.path.join(RESULTS_DIR,'prepare_submission','files','supplementary_tables'),
        # intermediate files
        ## networks: exon, bulk, single-cell, combined
        os.path.join(RESULTS_DIR,'prepare_submission','files','intermediate_files','networks'),
        ## models 
        ### weights
        expand(os.path.join(RESULTS_DIR,'prepare_submission','files','intermediate_files','models',"from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS, k=K_CROSS_VALIDATION),
        ### regulators
        expand(os.path.join(RESULTS_DIR,'prepare_submission',"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","input_regulators.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),
        expand(os.path.join(RESULTS_DIR,'prepare_submission',"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","output_regulators.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),
        ## datasets
        ### Rogalska2024
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-EX.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-ALTA.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-ALTD.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-INT.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","genexpr_tpm","Rogalska2024.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","Rogalska2024.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Rogalska2024-EX.tsv.gz"),
        ### Urbanski2022
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-EX.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-ALTA.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-ALTD.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-INT.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","genexpr_tpm","Urbanski2022.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","Urbanski2022.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Urbanski2022-EX.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Urbanski2022-genexpr_tpm.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Urbanski2022-EX.tsv.gz"),      
        ### Danielsson2013
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Danielsson2013-EX.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Danielsson2013-genexpr_tpm.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Danielsson2013-EX.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Danielsson2013-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        ### Hodis2022
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Hodis2022-invitro_eng_melanoc-genexpr_cpm.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Hodis2022-invitro_eng_melanoc-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","Hodis2022-invitro_eng_melanoc.tsv.gz"),
        ### Replogle2022
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","ReplogleWeissman2022_rpe1-genexpr_cpm.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","ReplogleWeissman2022_rpe1-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","ReplogleWeissman2022_rpe1.tsv.gz"),
        ### Cardoso-Moreira2019
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","CardosoMoreira2019-EX.tsv.gz"),
        os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","CardosoMoreira2019-EX.tsv.gz")        

        
rule supplementary_tables:
    input:
        # Identified cancer splicing programs
        suptab01_cancer_splicing_programs = os.path.join(PROGRAM_DIR,'files','PANCAN','cancer_program.tsv.gz'),
    output:
        directory(os.path.join(RESULTS_DIR,'prepare_submission','files','supplementary_tables'))
    run:
        import os
        import subprocess
        
        outdir = output[0]
        os.makedirs(outdir, exist_ok=True)
        
        for key, f in input.items():
            filename = os.path.basename(f)
            extension = ".".join(filename.split(".")[1:])
            outfile = os.path.join(outdir,key+"."+extension)
            cmd = ["cp", f, outfile]
            print(cmd)
            subprocess.call(cmd)
            
            if not filename.endswith(".gz"):
                cmd = ["gzip", outfile]
                print(cmd)
                subprocess.call(cmd)
            
        print("Done!")
        
        
rule intermediate_files_networks:
    input:
        exon = os.path.join(RESULTS_DIR,"new_empirical_network","files","experimentally_derived_regulons_pruned_w_viper_networks-EX"),
        bulk = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","experimentally_derived_regulons_pruned-bulkgenexpr"),
        singlecell = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","experimentally_derived_regulons_pruned-scgenexpr"),
        combined = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","experimentally_derived_regulons_pruned-bulkscgenexpr")
    output:
        directory(os.path.join(RESULTS_DIR,'prepare_submission','files','intermediate_files','networks')) 
    run:
        import os
        import subprocess
        
        output_dir = output[0]
        os.makedirs(output_dir, exist_ok=True)
        
        for key, regulon_dir in input.items():
            
            regulon_outdir = os.path.join(output_dir,os.path.basename(regulon_dir))
            os.makedirs(regulon_outdir, exist_ok=True)
            
            regulon_files = [os.path.join(regulon_dir,f) for f in os.listdir(regulon_dir) if f.endswith(".tsv.gz")]
            
            for regulon_file in regulon_files:
                filename = os.path.basename(regulon_file)
                outfile = os.path.join(regulon_outdir, filename)
                cmd = ["cp", regulon_file, outfile]
                print(cmd)
                subprocess.call(cmd)
            
        print("Done!")
        
        
rule intermediate_files_models_weights:
    input:
        weights = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth")
    output:
        weights = os.path.join(RESULTS_DIR,'prepare_submission','files','intermediate_files','models',"from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth")
    run:
        import os
        import subprocess
        
        cmd = ["cp", input.weights, output.weights]
        print(cmd)
        subprocess.call(cmd)
        
        print("Done!")
        
rule intermediate_files_models_regulators:
    input:
        input_regulators = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","input_regulators.tsv.gz"),
        output_regulators = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","output_regulators.tsv.gz")
    output:
        input_regulators = os.path.join(RESULTS_DIR,'prepare_submission',"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","input_regulators.tsv.gz"),
        output_regulators = os.path.join(RESULTS_DIR,'prepare_submission',"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","output_regulators.tsv.gz")
    run:
        import os
        import subprocess
        
        cmd = ["cp", input.input_regulators, output.input_regulators]
        print(cmd)
        subprocess.call(cmd)
        
        cmd = ["cp", input.output_regulators, output.output_regulators]
        print(cmd)
        subprocess.call(cmd)
        
        print("Done!")
        

rule intermediate_files_datasets_rogalska:
    input:
        psi_ex = os.path.join(PREP_DIR,"event_psi","Rogalska2024-EX.tsv.gz"),
        psi_alta = os.path.join(PREP_DIR,"event_psi","Rogalska2024-ALTA.tsv.gz"),
        psi_altd = os.path.join(PREP_DIR,"event_psi","Rogalska2024-ALTD.tsv.gz"),
        psi_int = os.path.join(PREP_DIR,"event_psi","Rogalska2024-INT.tsv.gz"),
        tpm = os.path.join(PREP_DIR,"genexpr_tpm","Rogalska2024.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Rogalska2024.tsv.gz"),
        signature_ex = os.path.join(PREP_DIR,'delta_psi','Rogalska2024-EX.tsv.gz')
    output:
        psi_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-EX.tsv.gz"),
        psi_alta = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-ALTA.tsv.gz"),
        psi_altd = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-ALTD.tsv.gz"),
        psi_int = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Rogalska2024-INT.tsv.gz"),
        tpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","genexpr_tpm","Rogalska2024.tsv.gz"),
        metadata = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","Rogalska2024.tsv.gz"),
        signature_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Rogalska2024-EX.tsv.gz"),
    run:
        import os
        import subprocess

        for key, f in input.items():
            outfile = output[key]
            cmd = ["cp", f, outfile]
            print(cmd)
            subprocess.call(cmd)
            
        print("Done!")
        
        
rule intermediate_files_datasets_urbanski:
    input:
        psi_ex = os.path.join(PREP_DIR,"event_psi","Urbanski2022-EX.tsv.gz"),
        psi_alta = os.path.join(PREP_DIR,"event_psi","Urbanski2022-ALTA.tsv.gz"),
        psi_altd = os.path.join(PREP_DIR,"event_psi","Urbanski2022-ALTD.tsv.gz"),
        psi_int = os.path.join(PREP_DIR,"event_psi","Urbanski2022-INT.tsv.gz"),
        tpm = os.path.join(PREP_DIR,"genexpr_tpm","Urbanski2022.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Urbanski2022.tsv.gz"),
        signature_ex = os.path.join(PREP_DIR,'signatures','Urbanski2022-EX.tsv.gz'),
        signature_tpm = os.path.join(PREP_DIR,'signatures','Urbanski2022-genexpr_tpm.tsv.gz'),
        activity_ex = os.path.join(RESULTS_DIR,"carcinogenic_switch_regulation","files","protein_activity","Urbanski2022-EX.tsv.gz")
    output:
        psi_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-EX.tsv.gz"),
        psi_alta = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-ALTA.tsv.gz"),
        psi_altd = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-ALTD.tsv.gz"),
        psi_int = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","event_psi","Urbanski2022-INT.tsv.gz"),
        tpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","genexpr_tpm","Urbanski2022.tsv.gz"),
        metadata = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","Urbanski2022.tsv.gz"),
        signature_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Urbanski2022-EX.tsv.gz"),
        signature_tpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Urbanski2022-genexpr_tpm.tsv.gz"),
        activity_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Urbanski2022-EX.tsv.gz")
    run:
        import os
        import subprocess

        for key, f in input.items():
            outfile = output[key]
            cmd = ["cp", f, outfile]
            print(cmd)
            subprocess.call(cmd)
            
        print("Done!")
        

rule intermediate_files_datasets_danielsson:
    input:
        signature_ex = os.path.join(RESULTS_DIR,"new_empirical_network","files","signatures","carcinogenesis-EX.tsv.gz"),
        signature_tpm = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","signatures","carcinogenesis-genexpr.tsv.gz"),
        activity_ex = os.path.join(RESULTS_DIR,"new_empirical_network","files","protein_activity","carcinogenesis-EX.tsv.gz"),
        activity_tpm = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","protein_activity","carcinogenesis-bulkgenexpr-adjusted_fclayer.tsv.gz")
    output:
        signature_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Danielsson2013-EX.tsv.gz"),
        signature_tpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Danielsson2013-genexpr_tpm.tsv.gz"),
        activity_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Danielsson2013-EX.tsv.gz"),
        activity_tpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Danielsson2013-bulkgenexpr-adjusted_fclayer.tsv.gz")
    run:
        import os
        import subprocess

        for key, f in input.items():
            outfile = output[key]
            cmd = ["cp", f, outfile]
            print(cmd)
            subprocess.call(cmd)
            
        print("Done!")

        
rule intermediate_files_datasets_hodis:
    input:
        signature_cpm = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","signatures","Hodis2022-invitro_eng_melanoc-genexpr.tsv.gz"),
        activity_cpm = os.path.join(RESULTS_DIR,"activity_estimation_w_genexpr","files","protein_activity","Hodis2022-invitro_eng_melanoc-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-conditions.tsv.gz")
    output:
        signature_cpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","Hodis2022-invitro_eng_melanoc-genexpr_cpm.tsv.gz"),
        activity_cpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","Hodis2022-invitro_eng_melanoc-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        metadata = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","Hodis2022-invitro_eng_melanoc.tsv.gz")
    run:
        import os
        import subprocess

        for key, f in input.items():
            outfile = output[key]
            cmd = ["cp", f, outfile]
            print(cmd)
            subprocess.call(cmd)
            
        print("Done!")


rule intermediate_files_datasets_replogle:
    input:
        signature_cpm = os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
        activity_cpm = os.path.join(RESULTS_DIR,"carcinogenic_switch_regulation","files","protein_activity","ReplogleWeissman2022_rpe1-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        metadata =  os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-conditions.tsv.gz")
    output:
        signature_cpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","ReplogleWeissman2022_rpe1-genexpr_cpm.tsv.gz"),
        activity_cpm = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","ReplogleWeissman2022_rpe1-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        metadata = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","metadata","ReplogleWeissman2022_rpe1.tsv.gz")
    run:
        import os
        import subprocess

        for key, f in input.items():
            outfile = output[key]
            cmd = ["cp", f, outfile]
            print(cmd)
            subprocess.call(cmd)
            
        print("Done!")
        
        
rule intermediate_files_datasets_cardosomoreira:
    input:
        signature_ex = os.path.join(RESULTS_DIR,"sf_programs_in_differentiation","files","signatures","CardosoMoreira2020-EX.tsv.gz"),
        activity_ex = os.path.join(RESULTS_DIR,"sf_programs_in_differentiation","files","protein_activity","CardosoMoreira2020-EX.tsv.gz"),
    output:
        signature_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","signatures","CardosoMoreira2019-EX.tsv.gz"),
        activity_ex = os.path.join(RESULTS_DIR,"prepare_submission","files","intermediate_files","datasets","activity","CardosoMoreira2019-EX.tsv.gz")
    run:
        import os
        import subprocess

        for key, f in input.items():
            outfile = output[key]
            cmd = ["cp", f, outfile]
            print(cmd)
            subprocess.call(cmd)
            
        print("Done!")

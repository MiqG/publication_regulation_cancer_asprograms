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

##### RULES #####
rule all:
    input:
        # supplementary tables
        os.path.join(RESULTS_DIR,'prepare_submission','files','supplementary_tables'),

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
        
        

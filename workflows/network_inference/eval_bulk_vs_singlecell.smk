import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
BIN_DIR = os.path.join(ROOT,"bin")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","network_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### RULES #####
rule all:
    input:
        # estimate splicing factor activity using 
        ## bulk samples
        ## single cell samples
        
        # make figures
        os.path.join(RESULTS_DIR,"figures","eval_genexpr_vs_splicing"),
        
    
rule figures_regulon_evaluation:
    input:
        evaluation_ex = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz"),
        evaluation_genexpr = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-genexpr.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","regulon_evaluation"))
    shell:
        """
        Rscript scripts/figures_regulon_evaluation.R \
                    --evaluation_ex_file={input.evaluation_ex} \
                    --evaluation_genexpr_file={input.evaluation_genexpr} \
                    --figs_dir={output}
        """

import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")

##### RULES #####
rule all:
    input:
        # preprocess STRINGDB
        os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz')
        
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
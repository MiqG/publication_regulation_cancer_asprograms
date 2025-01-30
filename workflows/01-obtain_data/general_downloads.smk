import os
import pandas as pd
import numpy as np

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
ARTICLES_DIR = os.path.join(RAW_DIR,'articles')
ARTICLE_DIR = os.path.join(ARTICLES_DIR,'Rogalska2024')
VASTDB_DIR = os.path.join(RAW_DIR,'VastDB')

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### RULES #####
rule all:
    input:
        # download metadata
        os.path.join(ARTICLE_DIR,"metadata.tsv"),
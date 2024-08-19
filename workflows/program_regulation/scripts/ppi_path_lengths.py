#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# For each drug, measure the shortest path length from each of its targets
# to each of the genes significantly associated to the drug sensitivity profile

import argparse
import pandas as pd
import numpy as np
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
FDR_THRESH = 0.1
N_JOBS = 1
RANDOM_SEED = 1234

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_viper_splicing'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
SUPPORT_DIR = os.path.join(ROOT,'support')
ppi_file = os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz')
targets_file = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-symbol.txt")
sources_file = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-symbol.txt")
n_jobs=10
n_random_sources = 100
"""

##### FUNCTIONS #####
def load_data(ppi_file, targets_file, sources_file):
    # load
    ppi = pd.read_table(ppi_file)
    targets = pd.read_table(targets_file, header=None)[0].to_list()
    sources = pd.read_table(sources_file, header=None)[0].to_list()

    return ppi, targets, sources


def prepare_data(ppi, targets, sources):
    # make network
    ppi = nx.from_pandas_edgelist(ppi)
    avail_genes = list(ppi.nodes())
    print("Total nodes in PPI network: %s" % len(avail_genes))
    
    # subset
    targets = list(set(targets).intersection(avail_genes))
    sources = list(set(sources).intersection(avail_genes))

    return ppi, targets, sources


def single_shortest_path(G, source, target, weight=None):
    try:
        l = nx.shortest_path_length(G, source=source, target=target, weight=weight)
    except nx.NetworkXNoPath:
        print("Unreachable from %s to %s." % (source, target))
        l = np.nan

    return {"source": source, "target": target, "shortest_path_length": l}


def compute_shortest_paths_real(ppi, targets, sources, n_jobs=None):

    pairs = [(source, target) for target in targets for source in sources]

    result = Parallel(n_jobs=n_jobs)(
        delayed(single_shortest_path)(ppi, s, t, weight="distance")
        for s, t in tqdm(pairs)
    )
    result = pd.DataFrame(result)
    result["type"] = "real"

    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ppi_file", type=str)
    parser.add_argument("--targets_file", type=str)
    parser.add_argument("--sources_file", type=str, default=None)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    ppi_file = args.ppi_file
    targets_file = args.targets_file
    sources_file = args.sources_file
    n_jobs = args.n_jobs
    output_file = args.output_file

    print(args)

    # load
    print("Loading data...")
    ppi, targets, sources = load_data(
        ppi_file, targets_file, sources_file
    )

    # prepare
    ppi, targets, sources = prepare_data(
        ppi, targets, sources
    )

    # compute shortest paths
    result = compute_shortest_paths_real(
        ppi, targets, sources, n_jobs
    )

    # save
    print("Saving data...")
    result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")


import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### RULES #####
rule all:
    input:
        # make regulons
        os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment0-delta_psi.tsv.gz"),
        os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment1-delta_psi.tsv.gz"),
        os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment2-delta_psi.tsv.gz"),
        os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment3-delta_psi.tsv.gz"),
        os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-HepG2-delta_psi.tsv.gz"),
        os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-K562-delta_psi.tsv.gz"),
        os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-HepG2-delta_psi.tsv.gz"),
        os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-K562-delta_psi.tsv.gz")

rule download_networks:
    params:
        ena_0 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENASFS-metaexperiment0-delta_psi.tsv.gz",
        ena_1 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENASFS-metaexperiment1-delta_psi.tsv.gz",
        ena_2 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENASFS-metaexperiment2-delta_psi.tsv.gz",
        ena_3 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENASFS-metaexperiment3-delta_psi.tsv.gz",
        encorekd_hepg2 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENCOREKD-HepG2-delta_psi.tsv.gz",
        encorekd_k562 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENCOREKD-K562-delta_psi.tsv.gz",
        encoreko_hepg2 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENCOREKO-HepG2-delta_psi.tsv.gz",
        encoreko_k562 = "https://github.com/MiqG/viper_splicing/raw/refs/heads/master/data/empirical_sf_networks-EX/ENCOREKO-K562-delta_psi.tsv.gz"
    output:
        ena_0 = os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment0-delta_psi.tsv.gz"),
        ena_1 = os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment1-delta_psi.tsv.gz"),
        ena_2 = os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment2-delta_psi.tsv.gz"),
        ena_3 = os.path.join(RAW_DIR,"viper_splicing_networks","ENASFS-metaexperiment3-delta_psi.tsv.gz"),
        encorekd_hepg2 = os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-HepG2-delta_psi.tsv.gz"),
        encorekd_k562 = os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKD-K562-delta_psi.tsv.gz"),
        encoreko_hepg2 = os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-HepG2-delta_psi.tsv.gz"),
        encoreko_k562 = os.path.join(RAW_DIR,"viper_splicing_networks","ENCOREKO-K562-delta_psi.tsv.gz")
    shell:
        """
        set -eo pipefail
        
        wget {params.ena_0} -O {output.ena_0}
        wget {params.ena_1} -O {output.ena_1}
        wget {params.ena_2} -O {output.ena_2}
        wget {params.ena_3} -O {output.ena_3}
        wget {params.encorekd_hepg2} -O {output.encorekd_hepg2}
        wget {params.encorekd_k562} -O {output.encorekd_k562}
        wget {params.encoreko_hepg2} -O {output.encoreko_hepg2}
        wget {params.encoreko_k562} -O {output.encoreko_k562}
        
        echo "Done!"
        """
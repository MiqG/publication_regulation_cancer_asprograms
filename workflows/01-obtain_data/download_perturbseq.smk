import os

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
DATA_DIR = os.path.join(ROOT,'data','raw')
SUPPORT_DIR = os.path.join(ROOT,'support')

##### RULES #####
rule all:
    input:
        ## ReplogleWeissman2022
        os.path.join(DATA_DIR,"scPerturb","ReplogleWeissman2022_K562_essential.h5ad"),
        os.path.join(DATA_DIR,"scPerturb","ReplogleWeissman2022_K562_gwps.h5ad"),
        os.path.join(DATA_DIR,"scPerturb","ReplogleWeissman2022_rpe1.h5ad")
        
rule download_scperturb:
    params:
        replogle_k562_essential = "https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad?download=1",
        replogle_k562_gwps = "https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_gwps.h5ad?download=1",
        replogle_rpe1 = "https://zenodo.org/records/7041849/files/ReplogleWeissman2022_rpe1.h5ad?download=1",
    output:
        replogle_k562_essential = os.path.join(DATA_DIR,"scPerturb","ReplogleWeissman2022_K562_essential.h5ad"),
        replogle_k562_gwps = os.path.join(DATA_DIR,"scPerturb","ReplogleWeissman2022_K562_gwps.h5ad"),
        replogle_rpe1 = os.path.join(DATA_DIR,"scPerturb","ReplogleWeissman2022_rpe1.h5ad")
    shell:
        """
        set -eo pipefail
        
        wget {params.replogle_k562_essential} -O {output.replogle_k562_essential}
        wget {params.replogle_k562_gwps} -O {output.replogle_k562_gwps}
        wget {params.replogle_rpe1} -O {output.replogle_rpe1}
        
        echo "Done!"
        """
        
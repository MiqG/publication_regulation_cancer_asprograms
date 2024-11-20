#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------

require(optparse)
require(tidyverse)
require(Seurat)
require(SeuratDisk)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# merged_seurat_file = file.path(RAW_DIR,"articles","Becker2021","Final_scHTAN_colon_all_epithelial_220213.rds")

##### FUNCTIONS #####
parseargs = function(){
    
    option_list = list( 
        make_option("--merged_seurat_file", type="character"),
        make_option("--adata_file", type="character"),
        make_option("--h5seurat_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    merged_seurat_file = args[["merged_seurat_file"]]
    adata_file = args[["adata_file"]]
    h5seurat_file = args[["h5seurat_file"]]
    
    # load
    seurat = readRDS(merged_seurat_file)
    seurat = UpdateSeuratObject(seurat)
    
    # normalize
    seurat = NormalizeData(seurat)

    # save
    SaveH5Seurat(seurat, filename = h5seurat_file)
    Convert(h5seurat_file, dest = adata_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
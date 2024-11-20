#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# make figures of evaluate target inference algorithms.
# median accuracy of 0.973 and 0.975 in HepG2 and K562 with threshold of 0.2

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(ggvenn)
require(extrafont)

# variables

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_SINGLE_DARK = "darkgreen"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","network_inference")
# SUPPORT_DIR = file.path(ROOT,"support")

# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
# networks_viper_alone_dir = file.path(RESULTS_DIR,"files","viper_networks-EX")
# networks_new_w_viper_dir = file.path(RESULTS_DIR,"files","experimentally_derived_regulons_pruned_w_viper_networks-EX")
# figs_dir = file.path(RESULTS_DIR,'figures','empirical_psi_networks')

##### FUNCTIONS #####
load_network_tsv = function(regulons_dir){
    net = lapply(list.files(regulons_dir, full.names=TRUE), function(regulons_file){
            regulon_id = basename(regulons_file) %>% gsub("-delta_psi.tsv.gz","",.)
            regulons = read_tsv(regulons_file) %>%
                mutate(regulon_id = regulon_id)
            return(regulons)
        }) %>% bind_rows()
    return(net)
}

plot_empirical_networks = function(splicing_factors){
    plts = list()
    
    X = splicing_factors
    
    sfs_oi = list(
        "Version 1" = X %>% filter(in_viper_alone) %>% pull(ENSEMBL),
        "Version 2" = X %>% filter(in_new_w_viper) %>% pull(ENSEMBL)
    )
    
    plts[["empirical_networks-versions-venn"]] = sfs_oi %>%
        ggvenn(
            fill_color = get_palette("rickandmorty",length(sfs_oi)),
            stroke_color = NA,
            set_name_size = FONT_SIZE+0.5,
            text_size = FONT_SIZE
        )
    
    return(plts)
}

make_plots = function(splicing_factors){
    plts = list(
        plot_empirical_networks(splicing_factors)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(splicing_factors){
    
    figdata = list(
        "eda" = list(
            "splicing_factors" = splicing_factors
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)   
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    save_plt(plts, "empirical_networks-versions-venn", '.pdf', figs_dir, width=5, height=5)
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        d = file.path(dir,'figdata',x)
        dir.create(d, recursive=TRUE)
        lapply(names(figdata[[x]]), function(nm){
            df = figdata[[x]][[nm]]
            filename = file.path(d, paste0(nm,'.tsv.gz'))
            write_tsv(df, filename)
            print(filename)
        })
    })
}


parseargs = function(){
    
    option_list = list( 
        make_option("--splicing_factors_file", type="character"),
        make_option("--networks_viper_alone_dir", type="character"),
        make_option("--networks_new_w_viper_dir", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    splicing_factors_file = args[["splicing_factors_file"]]
    networks_viper_alone_dir = args[["networks_viper_alone_dir"]]
    networks_new_w_viper_dir = args[["networks_new_w_viper_dir"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    splicing_factors = read_tsv(splicing_factors_file)
    networks_viper_alone = load_network_tsv(networks_viper_alone_dir)  
    networks_new_w_viper = load_network_tsv(networks_new_w_viper_dir)
    
    # prep
    splicing_factors = splicing_factors %>%
        mutate(
            in_viper_alone = ENSEMBL %in% networks_viper_alone[["regulator"]],
            in_new_w_viper = ENSEMBL %in% networks_new_w_viper[["regulator"]]
        )
    
    # stats
    ## overview
    n_total_sfs = splicing_factors %>% nrow()
    n_sfs_screened = splicing_factors %>% filter(in_viper_alone) %>% nrow()
    print(sprintf("VIPER splicing alone | Screened %s SFs out of %s", n_sfs_screened, n_total_sfs))
    n_sfs_screened = splicing_factors %>% filter(in_new_w_viper) %>% nrow()
    print(sprintf("New + VIPER splicing | Screened %s SFs out of %s", n_sfs_screened, n_total_sfs))
    
    # plot
    plts = make_plots(splicing_factors)
    
    # make figdata
    figdata = make_figdata(splicing_factors)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
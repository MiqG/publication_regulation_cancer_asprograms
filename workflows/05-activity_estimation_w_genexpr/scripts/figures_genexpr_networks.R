require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(ggvenn)
require(extrafont)
require(scattermore)

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
# RESULTS_DIR = file.path(ROOT,"results","activity_estimation_w_genexpr")
# SUPPORT_DIR = file.path(ROOT,"support")
# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
# networks_bulk_dir = file.path(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-bulkgenexpr")
# networks_singlecell_dir = file.path(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr")
# figs_dir = file.path(RESULTS_DIR,"figures","genexpr_networks")

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


make_plots = function(splicing_factors, signatures, metadata_encore){
    plts = list(
        plot_empirical_networks(splicing_factors),
        plot_encore(signatures, metadata_encore)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(splicing_factors, signatures, metadata_encore){
    
    figdata = list(
        "eda" = list(
            "splicing_factors" = splicing_factors
        ),
        "encore" = list(
            "signatures" = signatures,
            "metadata_encore" = metadata_encore
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
        make_option("--networks_bulk_dir", type="character"),
        make_option("--networks_singlecell_dir", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    splicing_factors_file = args[["splicing_factors_file"]]
    networks_bulk_dir = args[["networks_bulk_dir"]]
    networks_singlecell_dir = args[["networks_singlecell_dir"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    splicing_factors = read_tsv(splicing_factors_file)
    networks_bulk = load_network_tsv(networks_bulk_dir)  
    networks_singlecell = load_network_tsv(networks_singlecell_dir)
    
    # prep
    ## single cell from 2 datasets of Replogle study
    ## perturbs 473 different splicing factors and 7,840 genes
    ## 113,495 interactions
    networks_singlecell = networks_singlecell %>%
        filter(cell_line_name=="K562") %>%
        distinct(regulator, target, likelihood, tfmode)
    
    ## bulk from three studies ENCORE KD, ENCORE KO, PRJNA540831
    ## perturbs 157 different splicing factors and 13,889 genes
    ## 104,274 interactions
    networks_bulk = networks_bulk %>%
        filter(cell_line_name=="K562" | cell_line_name=="K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE") %>%
        distinct(regulator, target, likelihood, tfmode)
    
    networks = merge(networks_singlecell, networks_bulk, by=c("regulator","target"), suffixes=c("_singlecell","_bulk"), all=TRUE)
    
    ## Only 2,611 interactions in common
    ## 1,390 target genes in common
    ## 134 splicing factors in common
    networks %>% drop_na() %>% nrow()
    networks %>% drop_na() %>% count(target) %>% nrow()
    networks %>% drop_na() %>% count(regulator) %>% nrow()
    
    ## how many interactions do common splicing factors have?
    common_sfs = networks %>% drop_na() %>% pull(regulator)
    networks_singlecell %>% filter(regulator%in%common_sfs) %>% nrow()
    networks_bulk %>% filter(regulator%in%common_sfs) %>% nrow()
    
    ## Correlations
    networks %>% drop_na() %>% summarize(
        correl_pearson = cor(likelihood_singlecell, likelihood_bulk, method="pearson"), # 0.03153418
        correl_spearman = cor(likelihood_singlecell, likelihood_bulk, method="spearman"), # 0.06684701
        correl_pearson = cor(likelihood_singlecell*tfmode_singlecell, likelihood_bulk*tfmode_bulk, method="pearson"), # 0.4502824
        correl_spearman = cor(likelihood_singlecell*tfmode_singlecell, likelihood_bulk*tfmode_bulk, method="spearman") # 0.455378
    )
    
    ## sign of common interactions
    networks %>% drop_na() %>% count(tfmode_singlecell) # -1: 1062; 1: 1549
    networks %>% drop_na() %>% count(tfmode_bulk) # -1: 1389; 1: 1222
    networks %>% drop_na() %>% filter((tfmode_singlecell==tfmode_bulk) & (tfmode_bulk==1)) %>% nrow() # 1091
    networks %>% drop_na() %>% filter((tfmode_singlecell==tfmode_bulk) & (tfmode_bulk==-1)) %>% nrow() # 931
    
    # prep
    splicing_factors = splicing_factors %>%
        mutate(
            in_viper_alone = ENSEMBL %in% networks_viper_alone[["regulator"]],
            in_new_w_viper = ENSEMBL %in% networks_new_w_viper[["regulator"]]
        )
    
    # plot
    plts = make_plots(splicing_factors, signatures_encore, metadata_encore)
    
    # make figdata
    figdata = make_figdata(splicing_factors, signatures_encore, metadata_encore)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

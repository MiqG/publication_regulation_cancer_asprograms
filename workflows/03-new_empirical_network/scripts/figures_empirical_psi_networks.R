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
# RESULTS_DIR = file.path(ROOT,"results","new_empirical_network")
# SUPPORT_DIR = file.path(ROOT,"support")
# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
# networks_viper_alone_dir = file.path(RESULTS_DIR,"files","viper_networks-EX")
# networks_new_w_viper_dir = file.path(RESULTS_DIR,"files","experimentally_derived_regulons_pruned_w_viper_networks-EX")
# encore_metadata_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","benchmark","metadata_encorekd.tsv.gz")
# encore_ex_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_ex_encorekd.tsv.gz")
# encore_alta_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_alta_encorekd.tsv.gz")
# encore_altd_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_altd_encorekd.tsv.gz")
# encore_int_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","benchmark","psi_int_encorekd.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","empirical_psi_networks")

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

plot_encore = function(signatures_encore, metadata_encore){
    plts = list()

    X = signatures_encore
    
    cell_lines = c("K562","HepG2")
    for (cell_line_oi in cell_lines){
        reps = encore_metadata %>% filter(cell_line==cell_line_oi & PERT_GENE=="SF3B1") %>% pull(sampleID)
        rep1 = reps[1]
        rep2 = reps[2]
        
        plt = X %>%
            drop_na(.data[[rep1]], .data[[rep2]]) %>%
            group_by(event_type) %>%
            mutate(
                n_obs = n(),
                event_type = sprintf("%s (n=%s)", event_type, n_obs)
            ) %>%
            ungroup() %>%
            ggplot(aes(x=.data[[rep1]], y=.data[[rep2]])) +
            geom_scattermore(pixels=c(1000,1000), pointsize=5, color="black", alpha=0.3) +
            stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
            theme_pubr() +
            facet_wrap(~event_type) +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(
                x="DeltaPSI Replicate 1", y="DeltaPSI Replicate 2", 
                subtitle=sprintf("ENCORE KD of SF3B1 in %s by Event Type", cell_line_oi)
            )
        
        plts[[sprintf("encore-sf3b1_validation_by_event_type-%s-scatter", cell_line_oi)]] = plt
        
        plt = X %>%
            drop_na(.data[[rep1]], .data[[rep2]]) %>%
            mutate(
                n_obs = n(),
                event_type = sprintf("(n=%s)", n_obs)
            ) %>%
            ggplot(aes(x=.data[[rep1]], y=.data[[rep2]])) +
            geom_scattermore(pixels=c(1000,1000), pointsize=5, color="black", alpha=0.3) +
            stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
            theme_pubr() +
            facet_wrap(~event_type) +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(
                x="DeltaPSI Replicate 1", y="DeltaPSI Replicate 2", 
                subtitle=sprintf("ENCORE KD of SF3B1 in %s across all Event Types", cell_line_oi)
            )
        
        plts[[sprintf("encore-sf3b1_validation_all_event_types-%s-scatter", cell_line_oi)]] = plt
    }
    
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
    
    save_plt(plts, "encore-sf3b1_validation_by_event_type-K562-scatter", '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, "encore-sf3b1_validation_all_event_types-K562-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "encore-sf3b1_validation_by_event_type-HepG2-scatter", '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, "encore-sf3b1_validation_all_event_types-HepG2-scatter", '.pdf', figs_dir, width=4, height=4)
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
        make_option("--encore_metadata_file", type="character"),
        make_option("--encore_ex_file", type="character"),
        make_option("--encore_alta_file", type="character"),
        make_option("--encore_altd_file", type="character"),
        make_option("--encore_int_file", type="character"),
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
    encore_metadata_file = args[["encore_metadata_file"]]
    encore_ex_file = args[["encore_ex_file"]]
    encore_alta_file = args[["encore_alta_file"]]
    encore_altd_file = args[["encore_altd_file"]]
    encore_int_file = args[["encore_int_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    splicing_factors = read_tsv(splicing_factors_file)
    networks_viper_alone = load_network_tsv(networks_viper_alone_dir)  
    networks_new_w_viper = load_network_tsv(networks_new_w_viper_dir)
    
    encore_metadata = read_tsv(encore_metadata_file)
    encore_ex = read_tsv(encore_ex_file)
    encore_alta = read_tsv(encore_alta_file)
    encore_altd = read_tsv(encore_altd_file)
    encore_int = read_tsv(encore_int_file)
    encore = list(
        "EX" = encore_ex,
        "ALTA" = encore_alta,
        "ALTD" = encore_altd,
        "INT" = encore_int
    )
    
    # prep
    splicing_factors = splicing_factors %>%
        mutate(
            in_viper_alone = ENSEMBL %in% networks_viper_alone[["regulator"]],
            in_new_w_viper = ENSEMBL %in% networks_new_w_viper[["regulator"]]
        )
    
    ## ENCORE: get signatures for SF3B1 KDs
    perts = encore_metadata %>% filter(PERT_GENE=="SF3B1") %>% pull(sampleID) %>% unique()
    signatures_encore = lapply(event_types, function(event_type_oi){
        
        result = lapply(perts, function(pert_oi){
            ctls = encore_metadata %>%
                filter(sampleID==pert_oi) %>%
                separate_rows(control_samples, sep=",") %>%
                pull(control_samples)

            X = encore[[event_type_oi]] %>%
                column_to_rownames("EVENT") %>%
                as.data.frame()
            ctl_psi = rowMeans(X[ctls], na.rm=TRUE)
            pert_psi = X[pert_oi]
            dpsi = pert_psi - ctl_psi

            return(dpsi)
        }) %>% 
        do.call(cbind, .) %>% 
        filter(!if_all(everything(), is.na)) %>%
        rownames_to_column("EVENT") %>%
        mutate(event_type = event_type_oi)
        
        return(result)
    }) %>% bind_rows()
    
    # stats
    ## overview
    n_total_sfs = splicing_factors %>% nrow()
    n_sfs_screened = splicing_factors %>% filter(in_viper_alone) %>% nrow()
    print(sprintf("VIPER splicing alone | Screened %s SFs out of %s", n_sfs_screened, n_total_sfs))
    n_sfs_screened = splicing_factors %>% filter(in_new_w_viper) %>% nrow()
    print(sprintf("New + VIPER splicing | Screened %s SFs out of %s", n_sfs_screened, n_total_sfs))
    
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
require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(scattermore)

# variables

# formatting
LINE_SIZE = 0.25
FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","target_selection")
# SUPPORT_DIR = file.path(ROOT,"support")

# devel_junctions_counts_file = file.path(PREP_DIR,"junction_counts","CardosoMoreira2019-ERP109002.tsv.gz")
# devel_gene_counts_file = file.path(RAW_DIR,"recount3","gene_counts","CardosoMoreira2019-ERP109002.tsv.gz")
# devel_metadata_file = file.path(RAW_DIR,"recount3","metadata","CardosoMoreira2019-ERP109002.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,'figures','junction_prioritization')

##### FUNCTIONS #####
plot_devel_analysis = function(devel_analysis){
    plts = list()
    
    X = devel_analysis %>% 
        pivot_longer(c(in_tumor,in_common,in_healthy), names_to="junction_set", values_to="in_junction_set") %>%
        filter(in_junction_set)
    
    plts[["devel-dpc_vs_junction_count_vs_junction_set-scatter"]] = X %>%
        ggplot(aes(x=log10(time_norm+1), y=log2(junction_count+1))) +
        geom_scattermore(aes(color=junction_set), alpha=0.5, pixels=c(1000,1000), pointsize=5) +
        geom_smooth(method="lm", linetype="dashed", linewidth=LINE_SIZE, color="black") +
        color_palette("Dark2") +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~junction_set, scales="free_y") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="log10(Days Post Conception + 1)", y="log2(Junction Count + 1)")
    
    X = devel_analysis %>% filter(is_significant)
    
    plts[["devel-dpc_vs_junction_count_vs_ccle_junction_type-scatter"]] = X %>%
        ggplot(aes(x=log10(time_norm+1), y=log2(junction_count+1))) +
        geom_scattermore(aes(color=junction_type), alpha=0.8, pixels=c(1000,1000), pointsize=8) +
        geom_smooth(method="lm", linetype="dashed", linewidth=LINE_SIZE, color="black") +
        color_palette(PAL_JUNCTION_TYPES) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~junction_type, scales="free_y") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="log10(Days Post Conception + 1)", y="log2(Junction Count + 1)")
        
    return(plts)
}

make_plots = function(devel_analysis){
    plts = list(
        plot_devel_analysis(devel_analysis)
    )
    plts = do.call(c,plts)
    return(plts)
}

make_figdata = function(devel_analysis){
    
    figdata = list(
        "junction_selection" = list(
            "devel_analysis"= devel_analysis,
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension=".pdf", 
                    directory="", dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    # devel analysis
    save_plt(plts, "devel-dpc_vs_junction_count_vs_junction_set-scatter", ".pdf", figs_dir, width=12, height=4)
    save_plt(plts, "devel-dpc_vs_junction_count_vs_ccle_junction_type-scatter", ".pdf", figs_dir, width=12, height=4)
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
        make_option("--devel_junctions_counts", type="character"),
        make_option("--devel_metadata", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    devel_junctions_counts = args[["devel_junctions_counts"]]
    devel_metadata = args[["devel_metadata"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    devel_junctions_counts = read_tsv(devel_junctions_counts_file)
    devel_metadata = read_tsv(devel_metadata_file)
    
    gc()
    
    # prep
    ## presence of potential cancer driver neojunctions in development vs gtex
    ## consider as detected junctions with at least 10 read counts
    devel_metadata = devel_metadata %>% 
        distinct(external_id, sra.sample_title) %>%
        separate(
            sra.sample_title, 
            c("sampleID","organism","tissue","time","sex"), 
            sep="\\."
        ) %>%
        mutate(
            # fix carnegie stages
            time = case_when(
                time=="CS13" ~ "32dpc",
                time=="CS14" ~ "33dpc",
                time=="CS16" ~ "39dpc",
                time=="CS17" ~ "41dpc",
                time=="CS18" ~ "44dpc",
                time=="CS19" ~ "46dpc",
                time=="CS20" ~ "49dpc",
                time=="CS21" ~ "51dpc",
                time=="CS22" ~ "53dpc",
                time=="CS23" ~ "56dpc",
                TRUE ~ time
            ),
            time_units = case_when(
                str_detect(time, "dpc") ~ "dpc", # days post conception
                str_detect(time, "w") ~ "w", # weeks post conception
                str_detect(time, "dpb") ~ "dpb", # days post birth
                str_detect(time, "wpb") ~ "wpb", # weeks post birth
                str_detect(time, "mpb") ~ "mpb", # months post birth
                str_detect(time, "ypb") ~ "ypb", # years post birth
            )
        ) %>%
        # normalized time to days post conception
        rowwise() %>% 
        mutate(
            time_norm = as.numeric(gsub(time_units,"",time)),
            time_norm = case_when(
                time_units=="dpc" ~ time_norm,
                time_units=="w" ~ time_norm*7,
                time_units=="dpb" ~ (9*30)+time_norm,
                time_units=="wpb" ~ (9*30)+7*time_norm,
                time_units=="mpb" ~ (9*30)+30*time_norm,
                time_units=="ypb" ~ (9*30)+365*time_norm
            )
        ) %>%
        ungroup()
    
    devel_analysis = devel_junctions_counts %>%
        filter(JUNCTION %in% model_summaries[["EVENT"]] | JUNCTION %in% summary_junctions_gtex[["JUNCTION"]]) %>%
        pivot_longer(-JUNCTION, names_to="external_id", values_to="junction_count") %>% 
        left_join(devel_metadata, by="external_id") %>%
        mutate(
            condition = sprintf("%s-%s", tissue, time),
            in_tumor_in_ccle = JUNCTION %in% junctions_tumor_in_ccle[["JUNCTION"]],
            in_tumor = JUNCTION %in% junctions_tumor[["JUNCTION"]],
            in_common = JUNCTION %in% junctions_common[["JUNCTION"]],
            in_healthy = JUNCTION %in% junctions_healthy[["JUNCTION"]],
            is_significant = JUNCTION %in% (prioritization_all %>% filter(is_significant) %>% pull(EVENT))
        ) %>%
        # add junction type
        left_join(junction_types %>% distinct(JUNCTION, junction_type), by="JUNCTION")
    
    # plot
    plts = make_plots(devel_analysis)

    # make figdata
    figdata = make_figdata(devel_analysis)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
    save_tracks_for_igv(tracks, figs_dir)
}

##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

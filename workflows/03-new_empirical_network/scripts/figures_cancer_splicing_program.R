require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)

# variables
THRESH_FDR = 0.05
THRESH_N_SUM = 5.5

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DRIVER_TYPE = c(
    "Random Genes"="darkgreen",
    "Not Driver-like"="grey",
    "Oncogenic"="#F6AE2D",
    "Tumor suppressor"="#6C98B3"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","new_empirical_network")
# diff_activity_file = file.path(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# splicing_factors_file = file.path(RESULTS_DIR,"figures","empirical_psi_networks","figdata","eda","splicing_factors.tsv.gz")
# gene_annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","cancer_splicing_program")

##### FUNCTIONS #####
plot_driver_selection = function(driver_activity){
    plts = list()
    
    # SF activity
    X = driver_activity
    x = X %>%
        count(GENE, driver_type, only_new_network) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup()
    
    sf_order = x %>%
        pivot_wider(id_cols=c("GENE","only_new_network"), names_from="driver_type", values_from="n_sign", values_fill=0) %>%
        rowwise() %>%
        mutate(diff=sum(`Tumor suppressor`,`Oncogenic`)) %>%
        ungroup() %>%
        arrange(diff, only_new_network,`Tumor suppressor`,`Oncogenic`) %>%
        pull(GENE)
    
    sfs_oi = x %>% 
        slice_max(n_sum, n=5) %>% 
        filter(driver_type=="Oncogenic") %>%
        bind_rows(
            x %>% 
            slice_min(n_sum, n=5) %>% 
            filter(driver_type=="Tumor suppressor")
        )

    plts[["driver_selection-n_signif_vs_driver_type-activity-bar"]] = x %>%
        mutate(GENE = factor(GENE, levels=sf_order)) %>%
        ggbarplot(x="GENE", y="n", fill="only_new_network", color=NA) +
        geom_text(
            aes(label=GENE, color=only_new_network),
            sfs_oi %>% filter(driver_type=="Oncogenic"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=1, vjust=1, nudge_y=0.25
        ) +
        geom_text(
            aes(label=GENE, color=only_new_network),
            sfs_oi %>% filter(driver_type=="Tumor suppressor"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=45, hjust=0, vjust=0.5, nudge_y=0.25
        ) +
        geom_hline(yintercept=THRESH_N_SUM, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_text(
            aes(x=x, label=label),
            . %>% 
            filter(abs(n_sum)>THRESH_N_SUM) %>%
            group_by(GENE) %>%
            slice_max(n, n=1) %>%
            ungroup() %>%
            count(driver_type) %>% 
                mutate(
                    label=paste0("n=",n),
                    n=c(15,15),
                    x=c(120,40)
                ),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette("Dark2") +
        color_palette("Dark2") +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~driver_type, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Splicing Factor", y="Count", fill="SF from new Version", color="SF from new Version")
    
    return(plts)
}

make_plots = function(diff_activity){
    plts = list(
        plot_driver_selection(diff_activity)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(diff_activity){
    figdata = list(
        "cancer_program" = list(
            "enrichments_reactome" = diff_activity
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
    # cancer splicing program definition
    save_plt(plts, "driver_selection-n_signif_vs_driver_type-activity-bar", '.pdf', figs_dir, width=5, height=7)
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
        make_option("--diff_activity_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--gene_annotation_file", type="character"),
        make_option("--figs_dir", type="character"),
        make_option("--random_seed", type="integer", default=1234)
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    diff_activity_file = args[["diff_activity_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    gene_annotation_file = args[["gene_annotation_file"]]
    figs_dir = args[["figs_dir"]]
    random_seed = args[["random_seed"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_activity = read_tsv(diff_activity_file)
    splicing_factors = read_tsv(splicing_factors_file)
    gene_annotation = read_tsv(gene_annotation_file) %>%
        dplyr::rename(
            GENE = `Approved symbol`,
            ENSEMBL = `Ensembl gene ID`
        )
    
    # prep
    only_new_network = splicing_factors %>% filter(!in_viper_alone & in_new_w_viper) %>% pull(ENSEMBL)
    diff_activity = diff_activity %>%
        group_by(cancer_type) %>% # correct p-values for each type of cancer
        mutate(
            padj = p.adjust(pvalue, method="fdr"),
            is_significant = padj < THRESH_FDR,
            only_new_network = regulator %in% only_new_network
        ) %>%
        ungroup() %>%
        dplyr::rename("ENSEMBL"="regulator") %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = "ENSEMBL"
        )
    
    driver_activity = diff_activity %>%
        mutate(
            driver_type = ifelse(
                `condition_a-median`>0, "Oncogenic", "Tumor suppressor"
            )
        ) %>%
        filter(is_significant)    
    
    # plot
    plts = make_plots(driver_activity)
    
    # make figdata
    figdata = make_figdata(driver_activity)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
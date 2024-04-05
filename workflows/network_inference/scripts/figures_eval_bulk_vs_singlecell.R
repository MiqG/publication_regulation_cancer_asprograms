#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------

Sys.setenv(VROOM_CONNECTION_SIZE='1000000')
require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_EVAL_TYPE = c(
    "random" = "lightgrey",
    "real" = "orange"
)
PAL_DARK = "brown"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","network_inference")
# protein_activity_bulk_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKD_K562-genexpr.tsv.gz")
# protein_activity_singlecell_file = file.path(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_K562_essential-genexpr.tsv.gz")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","eval_bulk_vs_singlecell")

##### FUNCTIONS #####
plot_evaluation = function(protein_activity_example, cancer_program_activity){
    plts = list()
    
    # correlation activity Bulk vs SC within sample
    X = protein_activity_example %>%
        group_by(PERT_ENSEMBL) %>%
        summarize(
            correlation_pearson = cor(activity_bulk, activity_singlecell, method="pearson", use="pairwise.complete.obs")
        ) %>%
        ungroup()
    
    plts[["evaluation-activity_bulk_vs_singlecell-within-violin"]] = X %>%
        pivot_longer(-PERT_ENSEMBL, names_to="correlation_type", values_to="correlation") %>%
        ggviolin(x="correlation_type", y="correlation", fill="orange", color=NA) + 
        geom_boxplot(width=0.5, outlier.size=0.1, fill=NA) + 
        ylim(-1, 1) +
        geom_text(
            aes(y = -1, label=label), 
            . %>% 
            count(correlation_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Correlation Type", y="Correlation")
    
    sample_oi = X %>% slice_min(correlation_pearson, n=1) %>% pull(PERT_ENSEMBL)
    plts[["evaluation-activity_bulk_vs_singlecell-worst-scatter"]] = protein_activity_example %>%
        filter(PERT_ENSEMBL == sample_oi) %>%
        drop_na() %>%
        ggscatter(x="activity_bulk", y="activity_singlecell", size=1, alpha=0.5, color=PAL_DARK) +
        stat_cor(size=FONT_SIZE, family=FONT_FAMILY, method="pearson") +
        theme(aspect.ratio=1) +
        labs(x="Bulk Protein Activity", y="SC Protein Activity", subtitle=sample_oi)
        
    sample_oi = X %>% slice_max(correlation_pearson, n=1) %>% pull(PERT_ENSEMBL)
    plts[["evaluation-activity_bulk_vs_singlecell-best-scatter"]] = protein_activity_example %>%
        filter(PERT_ENSEMBL == sample_oi) %>%
        drop_na() %>%
        ggscatter(x="activity_bulk", y="activity_singlecell", size=1, alpha=0.5, color=PAL_DARK) +
        stat_cor(size=FONT_SIZE, family=FONT_FAMILY, method="pearson") +
        theme(aspect.ratio=1) +
        labs(x="Bulk Protein Activity", y="SC Protein Activity", subtitle=sample_oi)
    
    # correltion cancer splicing program median activity
    X = cancer_program_activity %>%
        pivot_longer(c(activity_bulk, activity_singlecell), names_to="activity_type", values_to="activity") %>%
        pivot_wider(id_cols=c("PERT_ENSEMBL","activity_type"), names_from="driver_type", values_from="activity") %>%
        mutate(activity_diff = `Oncogenic` - `Tumor suppressor`) %>%
        pivot_wider(id_cols="PERT_ENSEMBL", names_from="activity_type", values_from="activity_diff")
    
    plts[["evaluation-program_activity_diff_bulk_vs_singlecell-scatter"]] = X %>%
        ggscatter(x="activity_bulk", y="activity_singlecell", color=PAL_DARK, size=1, alpha=0.5) +
        geom_smooth(method="lm", size=LINE_SIZE, linetype="dashed", alpha=0.5, color="black", fill="lightgray") +        
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        theme(aspect.ratio=1) +
        labs(x="Bulk Cancer Splicing Programs Activity Difference", y="SC Cancer Splicing Programs Activity Difference")
        
    X = cancer_program_activity
    
    plts[["evaluation-program_activity_bulk_vs_singlecell-scatter"]] = X %>%
        ggscatter(x="activity_bulk", y="activity_singlecell", color=PAL_DARK, size=1, alpha=0.5) +
        geom_smooth(method="lm", size=LINE_SIZE, linetype="dashed", alpha=0.5, color="black", fill="lightgray") +        
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~driver_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Bulk Cancer Splicing Programs Protein Activity", y="SC Cancer Splicing Programs Protein Activity")
    
    return(plts)
}


make_plots = function(protein_activity_example, cancer_program_activity){
    plts = list(
        plot_evaluation(protein_activity_example, cancer_program_activity)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity_example, cancer_program_activity){
    figdata = list(
        "eval_bulk_vs_singlecell" = list(
            "protein_activity_example" = protein_activity_example,
            "cancer_program_activity" = cancer_program_activity
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
    save_plt(plts, "evaluation-activity_bulk_vs_singlecell-within-violin", '.pdf', figs_dir, width=2.5, height=6)
    save_plt(plts, "evaluation-activity_bulk_vs_singlecell-worst-scatter", '.pdf', figs_dir, width=3.5, height=3.5)
    save_plt(plts, "evaluation-activity_bulk_vs_singlecell-best-scatter", '.pdf', figs_dir, width=3.5, height=3.5)
    save_plt(plts, "evaluation-program_activity_diff_bulk_vs_singlecell-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "evaluation-program_activity_bulk_vs_singlecell-scatter", '.pdf', figs_dir, width=6, height=3.5)
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
        make_option("--protein_activity_bulk_file", type="character"),
        make_option("--protein_activity_singlecell_file", type="character"),
        make_option("--cancer_program_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    protein_activity_bulk_file = args[["protein_activity_bulk_file"]]
    protein_activity_singlecell_file = args[["protein_activity_singlecell_file"]]
    cancer_program_file = args[["cancer_program_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity_bulk = read_tsv(protein_activity_bulk_file)
    protein_activity_singlecell = read_tsv(protein_activity_singlecell_file)
    cancer_program = read_tsv(cancer_program_file)
    
    # prep
    protein_activity_example = protein_activity_bulk %>%
        pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_bulk") %>%
        left_join(
            protein_activity_singlecell %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_singlecell"),
            by = c("PERT_ENSEMBL","regulator")
        ) %>%
        drop_na()
    
    cancer_program_activity = cancer_program %>%
        left_join(
            protein_activity_example, 
            by=c("ENSEMBL"="regulator")
        ) %>%
        pivot_longer(c(activity_bulk, activity_singlecell), names_to="activity_type", values_to="activity") %>%
        group_by(PERT_ENSEMBL, activity_type, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        drop_na() %>%
        pivot_wider(id_cols=c("PERT_ENSEMBL","driver_type"), names_from="activity_type", values_from="activity")
    
    # plot
    plts = make_plots(protein_activity_example, cancer_program_activity)
    
    # make figdata
    figdata = make_figdata(protein_activity_example, cancer_program_activity)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
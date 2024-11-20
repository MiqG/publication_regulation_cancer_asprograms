#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables
SETS_MAIN = c(
    'viper_networks',
    'experimentally_derived_regulons_pruned_w_viper_networks'
)

METHODS_ACTIVITY = c(
    "gsea",
    "correlation_spearman",
    "correlation_pearson",
    "viper"
)

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "#616161"
PAL_EVAL_TYPE = c(
    "random" = "lightgrey",
    "real" = "orange"
)

PAL_METHODS_ACTIVITY = c(
    "gsea"="#383F51",
    "correlation_spearman"="#B9BAA3",
    "correlation_pearson"="#8E8DBE",
    "viper"="#A22C29"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","network_inference")
# evaluation_file = file.path(RESULTS_DIR,"files","network_evaluation_scores","merged.tsv.gz")
# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","network_evaluation")

##### FUNCTIONS #####
plot_evaluation = function(evaluation){
    plts = list()
    
    X = evaluation
    
    # evaluation by dataset
    x = X %>%
        group_by(regulon_set, method_activity, curves_type, eval_direction, eval_type, 
                 signature_id) %>%
        summarize(auc_roc = median(auc_roc, na.rm=TRUE)) %>%    
        ungroup()
    
    plts[["evaluation-general-median_auc_roc-box"]] = x %>%
        filter(regulon_set%in%SETS_MAIN & eval_type=="real") %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_MAIN),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        ggplot(aes(x=regulon_set, y=auc_roc, group=interaction(regulon_set, method_activity))) +
        geom_boxplot(aes(color=method_activity), fill=NA, outlier.shape=NA, position=position_dodge(0.9), width=0.8) +
        geom_point(aes(color=method_activity, shape=signature_id), size=0.5, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.8)) +
        color_palette(PAL_METHODS_ACTIVITY) +
        geom_text(
            aes(y = 0.7, label=label), 
            . %>% 
            count(method_activity, regulon_set, eval_type, eval_direction) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", label="p.format", family=FONT_FAMILY, size=FONT_SIZE) + 
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~eval_direction) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="median(ROC AUC)", color="Inference Type", shape="Held-out Dataset")
    
    return(plts)
}

make_plots = function(evaluation){
    plts = list(
        plot_evaluation(evaluation)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(evaluation){
    figdata = list(
        "regulon_evaluation" = list(
            "evaluation" = evaluation
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
    
    # overall
    save_plt(plts, "evaluation-general-median_auc_roc-box", '.pdf', figs_dir, width=4, height=10)
    
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
        make_option("--evaluation_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    evaluation_file = args[["evaluation_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    evaluation = read_tsv(evaluation_file)
    splicing_factors = read_tsv(splicing_factors_file)
    
    # prep
    evaluation = evaluation %>%
        mutate(
            regulon_set = gsub("-.*","",regulon_set_id)
        ) %>% 
        # drop summarized evaluation
        filter(curves_type=="by_group" & n_tails=="two") %>%
        separate(PERT_ID, c("study_accession","cell_line","PERT_ENSEMBL","PERT_TYPE"), remove=FALSE, sep="___") %>%
        left_join(splicing_factors, by=c("PERT_ENSEMBL"="ENSEMBL")) %>%
        mutate(
            sf_class = case_when(
                !is.na(spliceosome_db_complex) ~ "Core",
                in_go_rbp ~ "RBP",
                TRUE ~ "Other"
            )
        )
    
    # plot
    plts = make_plots(evaluation)
    
    # make figdata
    figdata = make_figdata(evaluation)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
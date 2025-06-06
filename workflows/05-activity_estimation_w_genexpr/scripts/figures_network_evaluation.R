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
    'experimentally_derived_regulons_pruned_w_viper_networks-EX',
    'experimentally_derived_regulons_pruned-bulkgenexpr',
    'experimentally_derived_regulons_pruned-bulkscgenexpr',
    'experimentally_derived_regulons_pruned-scgenexpr'
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
# RESULTS_DIR = file.path(ROOT,"results","activity_estimation_w_genexpr")
# evaluation_file = file.path(RESULTS_DIR,"files","network_evaluation_scores","merged.tsv.gz")
# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","network_evaluation")

##### FUNCTIONS #####
plot_evaluation = function(evaluation){
    plts = list()
    
    X = evaluation
    
    # evaluation by dataset
    x = X %>%
        distinct(regulon_set, regulon_set_id, method_activity, curves_type, eval_direction, eval_type, 
                 signature_id, auc_roc) %>%
        group_by(regulon_set, regulon_set_id, method_activity, curves_type, eval_direction, eval_type, 
                 signature_id) %>%
        summarize(auc_roc = median(auc_roc, na.rm=TRUE)) %>%    
        ungroup() %>%
        mutate(benchmark_type = ifelse(str_detect(signature_id,"ReplogleWeissman2022"), "Single-cell", "Bulk"))
    
    plts[["evaluation-general-median_auc_roc-box"]] = x %>%
        filter(regulon_set_id%in%SETS_MAIN & eval_type=="real") %>%
        mutate(
            regulon_set_id = factor(regulon_set_id, levels=SETS_MAIN)
        ) %>%
        ggplot(aes(x=regulon_set_id, y=auc_roc, group=interaction(regulon_set_id, benchmark_type))) +
        geom_col(
            aes(x=regulon_set_id, y=auc_roc_avg, fill=benchmark_type), 
            . %>% group_by(benchmark_type, regulon_set_id, eval_type, eval_direction) %>% summarize(auc_roc_avg=mean(auc_roc)) %>% ungroup(),
            position=position_dodge(0.9)
        ) +
        geom_point(aes(color=signature_id, group=benchmark_type), size=0.25, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.3)) +
        color_palette("futurama") +
        fill_palette("simpsons") +
        geom_text(
            aes(y = 0.1, label=label), 
            . %>% 
            count(benchmark_type, regulon_set_id, eval_type, eval_direction) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~eval_direction) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="median(ROC AUC)", color="Data Type", shape="Held-out Dataset")
    
    # evaluation context independence
    x = X %>%
        filter(regulon_set_id%in%SETS_MAIN & eval_type=="real" & eval_direction=="perturbations") %>%
        # consider only cell lines that are not in networks
        filter(!(cell_line%in%c("K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","K562","HepG2"))) %>%
        distinct(cell_line, signature_id, auc_roc, method_activity, eval_type, eval_direction,
                 regulon_set_id, PERT_ENSEMBL, GENE, PERT_TYPE, study_accession) %>%
        rowwise() %>%
        mutate(
            cell_line = gsub("RPE1", "RPE1_RETINAL", cell_line),
            cell_line = gsub("K562", "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", cell_line),
            cell_line = gsub("HepG2", "HEPG2_LIVER", cell_line),
            cell_lineage = paste(unlist(strsplit(cell_line, "_"))[-1], collapse="_"),
            benchmark_type = ifelse(str_detect(signature_id,"ReplogleWeissman2022"), "Single-cell", "Bulk")
        ) %>%
        ungroup()
    
    cell_lineages_oi = x %>%
        count(signature_id, benchmark_type, eval_type, eval_direction, cell_lineage) %>%
        arrange(n) %>%
        filter(n>=10) %>%
        distinct(cell_lineage) %>%
        pull(cell_lineage)
    
    plts[["evaluation-context_indep-raw_auc_roc-box"]] = x %>%
        filter(cell_lineage%in%cell_lineages_oi) %>%
        mutate(
            cell_lineage = factor(cell_lineage, levels=rev(cell_lineages_oi)),
            regulon_set_id = factor(regulon_set_id, levels=SETS_MAIN)
        ) %>%
        ggplot(aes(x=regulon_set_id, y=auc_roc, group=interaction(regulon_set_id, benchmark_type))) +
        geom_boxplot(aes(color=benchmark_type), fill=NA, outlier.shape=NA, position=position_dodge(0.9)) +
        geom_point(aes(color=benchmark_type), size=0.25, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        color_palette("simpsons") + 
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(benchmark_type, regulon_set_id, 
                  eval_type, eval_direction, cell_lineage) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_hline(yintercept=0.5, linewidth=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~cell_lineage, nrow=2) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="ROC AUC", color="Held-out Data Type")
    
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
    save_plt(plts, "evaluation-general-median_auc_roc-box", '.pdf', figs_dir, width=6, height=10)
    save_plt(plts, "evaluation-context_indep-raw_auc_roc-box", '.pdf', figs_dir, width=13, height=12)
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
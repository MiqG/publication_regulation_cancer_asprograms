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
REGULON_SETS = c(
    'aracne_regulons_development',
    'mlr_regulons_development',
    'experimentally_derived_regulons_pruned',
    'aracne_and_experimental_regulons',
    'mlr_and_experimental_regulons',
    'aracne_and_mlr_regulons',
    'top100_experimentally_derived_regulons_pruned',
    'top90_experimentally_derived_regulons_pruned',
    'top80_experimentally_derived_regulons_pruned',
    'top70_experimentally_derived_regulons_pruned',
    'top60_experimentally_derived_regulons_pruned',
    'top50_experimentally_derived_regulons_pruned',
    'top40_experimentally_derived_regulons_pruned'
)

SETS_MAIN = c(
    'aracne_regulons_development',
    'mlr_regulons_development',
    'experimentally_derived_regulons_pruned'
)

SETS_ROBUSTNESS = c(
    'top40_experimentally_derived_regulons_pruned',
    'top50_experimentally_derived_regulons_pruned',
    'top60_experimentally_derived_regulons_pruned',
    'top70_experimentally_derived_regulons_pruned',
    'top80_experimentally_derived_regulons_pruned',
    'top90_experimentally_derived_regulons_pruned',
    'top100_experimentally_derived_regulons_pruned',
    'experimentally_derived_regulons_pruned'
)

SETS_LIKELIHOOD = c(
    'aracne_and_experimental_regulons',
    'mlr_and_experimental_regulons',
    'experimentally_derived_regulons_pruned'
)

SETS_MOR = c(
    'aracne_and_mlr_regulons'
)

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_EVAL_TYPE = c(
    "random" = "lightgrey",
    "real" = "orange"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# evaluation_ex_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz")
# evaluation_genexpr_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-genexpr.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","regulon_evaluation")

##### FUNCTIONS #####
plot_evaluation = function(evaluation, omic_type_oi){
    plts = list()
    
    X = evaluation %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, n_tails, regulon_set_id, pert_type_lab, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        filter(omic_type==omic_type_oi)
    
    # main networks
    plts[["evaluation-ranking_perc_vs_regulon_set_vs_pert_type-main-box"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(n_tails=="two") %>%
        ggplot(aes(x=pert_type_lab, y=ranking_perc, 
                   group=interaction(pert_type_lab, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~eval_direction+regulon_set_id, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(omic_type, pert_type_lab, regulon_set_id, eval_direction, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Validation Perturbation", y="Evaluation Score", fill="Inference Type")

    
    plts[["evaluation-ranking_perc_vs_regulon_set-main-box"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_MAIN)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    
    # robustness networks
    plts[["evaluation-ranking_perc_vs_regulon_set-robustness-box"]] = X %>%
        filter(regulon_set %in% SETS_ROBUSTNESS) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_ROBUSTNESS)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    # likelihood networks
    plts[["evaluation-ranking_perc_vs_regulon_set-likelihood-box"]] = X %>%
        filter(regulon_set %in% SETS_LIKELIHOOD) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_LIKELIHOOD)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    # mor networks
    plts[["evaluation-ranking_perc_vs_regulon_set-mor-box"]] = X %>%
        filter(regulon_set %in% SETS_MOR) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_MOR)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    # one-tailed networks
    plts[["evaluation-ranking_perc_vs_regulon_set-main_one_tailed-box"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(n_tails=="one") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_MAIN)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    names(plts) = sprintf("%s-%s", omic_type_oi, names(plts))
    
    return(plts)
}


make_plots = function(evaluation){
    plts = list(
        plot_evaluation(evaluation, "EX"),
        plot_evaluation(evaluation, "genexpr")
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
    omic_types = c("EX","genexpr")
    for (omic_type_oi in omic_types){
        # main
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set_vs_pert_type-main-box", omic_type_oi), '.pdf', figs_dir, width=6.5, height=12)
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-main-box", omic_type_oi), '.pdf', figs_dir, width=7, height=7)
        # robustness
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-robustness-box", omic_type_oi), '.pdf', figs_dir, width=12, height=7)
        # likelihood
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-likelihood-box", omic_type_oi), '.pdf', figs_dir, width=7, height=7)
        # mor
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-mor-box", omic_type_oi), '.pdf', figs_dir, width=3.5, height=7)
        # one-tailed
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-main_one_tailed-box", omic_type_oi), '.pdf', figs_dir, width=7, height=7)
    }
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
        make_option("--evaluation_ex_file", type="character"),
        make_option("--evaluation_genexpr_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    evaluation_ex_file = args[["evaluation_ex_file"]]
    evaluation_genexpr_file = args[["evaluation_genexpr_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    evaluation = list(
        read_tsv(evaluation_ex_file),
        read_tsv(evaluation_genexpr_file)
    ) %>%
    bind_rows()
    
    # prep
    evaluation = evaluation %>%
        mutate(regulon_id = gsub("-","_",regulon_id)) %>%
        filter(signature_id!=regulon_id) %>%
        filter(!(str_detect(regulon_id,"ENASFS") & (signature_id=="ENASFS"))) %>%
        # consider only signatures that we know activity 
        # of the splicing factor was altered
        filter(PERT_TYPE %in% c("KNOCKDOWN","KNOCKOUT","OVEREXPRESSION")) %>%
        mutate(
            pert_type_lab = case_when(
                PERT_TYPE=="KNOCKDOWN" ~ "KD",
                PERT_TYPE=="KNOCKOUT" ~ "KO",
                PERT_TYPE=="OVEREXPRESSION" ~ "OE"
            ),
            regulon_set = gsub("-.*","",regulon_set_id)
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
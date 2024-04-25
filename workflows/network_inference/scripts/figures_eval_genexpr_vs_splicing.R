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
# evaluation_ex_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz")
# evaluation_genexpr_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-genexpr.tsv.gz")
# evaluation_scgenexpr_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-scgenexpr.tsv.gz")
# protein_activity_ex_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-EX.tsv.gz")
# protein_activity_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-genexpr.tsv.gz")
# protein_activity_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-scgenexpr.tsv.gz")
# protein_activity_ew_model_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-EX_from_model_ewlayer_and_genexpr.tsv.gz")
# protein_activity_ew_model_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-EX_from_model_ewlayer_and_scgenexpr.tsv.gz")
# protein_activity_fc_model_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-EX_from_model_fclayer_and_genexpr.tsv.gz")
# protein_activity_fc_model_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","ENCOREKO_K562-EX_from_model_fclayer_and_scgenexpr.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","network_evaluation")

##### FUNCTIONS #####
plot_evaluation = function(evaluation, protein_activity_example){
    plts = list()
    
    X = evaluation %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, rnaseq_type,
                 n_tails, regulon_set_id, pert_type_lab, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() 
    
    # main networks
    plts[["evaluation-ranking_perc_vs_regulon_set_vs_pert_type-main-box"]] = X %>%
        ggplot(aes(x=pert_type_lab, y=ranking_perc, 
                   group=interaction(pert_type_lab, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~regulon_set_id+eval_direction+rnaseq_type, ncol=4) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(omic_type, pert_type_lab, regulon_set_id, eval_direction, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Validation Perturbation", y="Recall", fill="Inference Type")

    
    plts[["evaluation-ranking_perc_vs_regulon_set-main-box"]] = X %>% 
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator, rnaseq_type) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        ggplot(aes(x=omic_type, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~eval_direction+rnaseq_type, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Recall", fill="Inference Type")
    
    # correlation activity SF-exon vs SF-gene within sample
    X = protein_activity_example %>%
        drop_na() %>%
        group_by(PERT_ENSEMBL) %>%
        summarize(
            ex_vs_genexpr = cor(activity_ex, activity_genexpr, method="pearson"),
            ex_vs_scgenexpr = cor(activity_ex, activity_scgenexpr, method="pearson"),
            ex_vs_ew_model_genexpr = cor(activity_ex, activity_ew_model_genexpr, method="pearson"),
            ex_vs_ew_model_scgenexpr = cor(activity_ex, activity_ew_model_scgenexpr, method="pearson"),
            ex_vs_fc_model_genexpr = cor(activity_ex, activity_fc_model_genexpr, method="pearson"),
            ex_vs_fc_model_scgenexpr = cor(activity_ex, activity_fc_model_scgenexpr, method="pearson")
        ) %>%
        ungroup()
    
    plts[["evaluation-activity_ex_vs_genexpr-within-violin"]] = X %>%
        pivot_longer(-PERT_ENSEMBL, names_to="correlation_type", values_to="correlation") %>%
        drop_na() %>%
        ggviolin(x="correlation_type", y="correlation", fill="orange", color=NA, trim=TRUE) + 
        geom_boxplot(width=0.5, outlier.size=0.1, fill=NA) + 
        geom_text(
            aes(y = 1, label=label), 
            . %>% 
            count(correlation_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 25) + 
        labs(x="Correlation Type", y="Correlation")
    
#     # highlight correlations
#     sample_oi = X %>% slice_min(ex_vs_scgenexpr, n=1) %>% pull(PERT_ENSEMBL)
#     plts[["evaluation-activity_ex_vs_genexpr-worst-scatter"]] = protein_activity_example %>%
#         filter(PERT_ENSEMBL == sample_oi) %>%
#         drop_na() %>%
#         ggscatter(x="activity_ex", y="activity_scgenexpr", size=1, alpha=0.5, color=PAL_DARK) +
#         stat_cor(size=FONT_SIZE+2, family=FONT_FAMILY, method="pearson") +
#         theme(aspect.ratio=1) +
#         labs(x="SF-exon Protein Activity", y="SF-gene Protein Activity", subtitle=sample_oi)
        
#     sample_oi = X %>% slice_max(ex_vs_scgenexpr, n=1) %>% pull(PERT_ENSEMBL)
#     plts[["evaluation-activity_ex_vs_genexpr-best-scatter"]] = protein_activity_example %>%
#         filter(PERT_ENSEMBL == sample_oi) %>%
#         drop_na() %>%
#         ggscatter(x="activity_ex", y="activity_scgenexpr", size=1, alpha=0.5, color=PAL_DARK) +
#         stat_cor(size=FONT_SIZE, family=FONT_FAMILY, method="pearson") +
#         theme(aspect.ratio=1) +
#         labs(x="SF-exon Protein Activity", y="SF-gene Protein Activity", subtitle=sample_oi)
    
    return(plts)
}


make_plots = function(evaluation, protein_activity_example){
    plts = list(
        plot_evaluation(evaluation, protein_activity_example)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(evaluation, protein_activity_example){
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
    # main
    save_plt(plts, "evaluation-ranking_perc_vs_regulon_set_vs_pert_type-main-box", '.pdf', figs_dir, width=7.5, height=12)
    save_plt(plts, "evaluation-ranking_perc_vs_regulon_set-main-box", '.pdf', figs_dir, width=7, height=10)
    save_plt(plts, "evaluation-activity_ex_vs_genexpr-within-violin", '.pdf', figs_dir, width=5.5, height=5)
    save_plt(plts, "evaluation-activity_ex_vs_genexpr-worst-scatter", '.pdf', figs_dir, width=3.5, height=3.5)
    save_plt(plts, "evaluation-activity_ex_vs_genexpr-best-scatter", '.pdf', figs_dir, width=3.5, height=3.5)
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
        make_option("--protein_activity_ex_file", type="character"),
        make_option("--protein_activity_genexpr_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    evaluation_ex_file = args[["evaluation_ex_file"]]
    evaluation_genexpr_file = args[["evaluation_genexpr_file"]]
    protein_activity_ex_file = args[["protein_activity_ex_file"]]
    protein_activity_genexpr_file = args[["protein_activity_genexpr_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    evaluation = list(
        read_tsv(evaluation_ex_file),
        read_tsv(evaluation_genexpr_file),
        read_tsv(evaluation_scgenexpr_file)
    ) %>%
    bind_rows()
    protein_activity_ex = read_tsv(protein_activity_ex_file)
    protein_activity_genexpr = read_tsv(protein_activity_genexpr_file)
    protein_activity_scgenexpr = read_tsv(protein_activity_scgenexpr_file)
    protein_activity_ew_model_genexpr = read_tsv(protein_activity_ew_model_genexpr_file)
    protein_activity_ew_model_scgenexpr = read_tsv(protein_activity_ew_model_scgenexpr_file)
    protein_activity_fc_model_genexpr = read_tsv(protein_activity_fc_model_genexpr_file)
    protein_activity_fc_model_scgenexpr = read_tsv(protein_activity_fc_model_scgenexpr_file)
    
    # prep
    evaluation = evaluation %>%
        mutate(
            signature_id = gsub("-pseudobulk_across_batches","",signature_id),
            regulon_id = gsub("-","_",regulon_id),
            regulon_id = gsub("_log2fc_genexpr","",regulon_id),
            regulon_id = gsub("_delta_psi","",regulon_id)
        ) %>%
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
            rnaseq_type = case_when(
                str_detect(signature_id, "Replogle") ~ "scRNAseq",
                TRUE ~ "bulkRNAseq"
            ),
            regulon_set = gsub(".*-","",regulon_set_id)
        )
    
    protein_activity_example = protein_activity_scgenexpr %>%
        pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_scgenexpr") %>%
        left_join(
            protein_activity_genexpr %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_genexpr"),
            by = c("PERT_ENSEMBL","regulator")
        ) %>%
        left_join(
            protein_activity_ex %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_ex"),
            by = c("PERT_ENSEMBL","regulator")
        ) %>%
        left_join(
            protein_activity_ew_model_genexpr %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_ew_model_genexpr"),
            by = c("PERT_ENSEMBL","regulator")
        ) %>%
        left_join(
            protein_activity_ew_model_scgenexpr %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_ew_model_scgenexpr"),
            by = c("PERT_ENSEMBL","regulator")
        ) %>%
        left_join(
            protein_activity_fc_model_genexpr %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_fc_model_genexpr"),
            by = c("PERT_ENSEMBL","regulator")
        ) %>%
        left_join(
            protein_activity_fc_model_scgenexpr %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_fc_model_scgenexpr"),
            by = c("PERT_ENSEMBL","regulator")
        )
    
    # plot
    plts = make_plots(evaluation, protein_activity_example)
    
    # make figdata
    figdata = make_figdata(evaluation, protein_activity_example)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
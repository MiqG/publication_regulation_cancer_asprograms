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

PAL_DRIVER_TYPE = c(
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","network_inference")
# protein_activity_bulk_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-genexpr.tsv.gz")
# protein_activity_singlecell_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-scgenexpr.tsv.gz")
# metadata_singlecell_file = file.path(PREP_DIR,"singlecell","ReplogleWeissman2022_K562_essential-conditions.tsv.gz")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","eval_bulk_vs_singlecell")
PREP_VIPER_DIR = file.path(dirname(ROOT),"publication_viper_splicing","data","prep")
metadata_file = file.path(PREP_VIPER_DIR,"metadata","tumorigenesis.tsv.gz")

##### FUNCTIONS #####
plot_tumorigenesis = function(protein_activity){
    plts = list()
    
    X = protein_activity %>%
        filter(study_accession=="PRJNA193487") %>%
        drop_na(driver_type) %>%
        group_by(cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(
            activity_singlecell = median(activity_singlecell),
            activity_bulk = median(activity_bulk)
        ) %>%
        ungroup() %>%
        mutate(cell_line_name=factor(
            cell_line_name, levels=c("BJ_PRIMARY","BJ_IMMORTALIZED",
                                     "BJ_TRANSFORMED","BJ_METASTATIC")
        ))
    
    plts[["tumorigenesis-cell_line_vs_activity-violin"]] = X %>%
        pivot_longer(c(activity_singlecell,activity_bulk), names_to="activity_type", values_to="activity") %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=activity, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        facet_wrap(~activity_type, nrow=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) + 
        theme_pubr() +
        labs(x="Cell Line", y="Protein Activity", fill="Driver Type")

    return(plts)
}


plot_evaluation = function(protein_activity, cancer_program_activity){
    plts = list()
    
    # correlation activity Bulk vs SC within sample
    X = protein_activity %>%
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
    plts[["evaluation-activity_bulk_vs_singlecell-worst-scatter"]] = protein_activity %>%
        filter(PERT_ENSEMBL == sample_oi) %>%
        drop_na() %>%
        ggscatter(x="activity_bulk", y="activity_singlecell", size=1, alpha=0.5, color=PAL_DARK) +
        stat_cor(size=FONT_SIZE, family=FONT_FAMILY, method="pearson") +
        theme(aspect.ratio=1) +
        labs(x="Bulk Protein Activity", y="SC Protein Activity", subtitle=sample_oi)
        
    sample_oi = X %>% slice_max(correlation_pearson, n=1) %>% pull(PERT_ENSEMBL)
    plts[["evaluation-activity_bulk_vs_singlecell-best-scatter"]] = protein_activity %>%
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


make_plots = function(protein_activity, cancer_program_activity){
    plts = list(
        plot_evaluation(protein_activity, cancer_program_activity)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity, cancer_program_activity){
    figdata = list(
        "eval_bulk_vs_singlecell" = list(
            "protein_activity" = protein_activity,
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
    metadata = read_tsv(metadata_file)
    metadata_singlecell = read_tsv(metadata_singlecell_file)
    cancer_program = read_tsv(cancer_program_file)
    
    # prep
    protein_activity = protein_activity_singlecell %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity_singlecell") %>%
        left_join(
            protein_activity_bulk %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_bulk"),
            by = c("sampleID","regulator")
        ) %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, activity_singlecell, activity_bulk) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s%s) (%s%s) | %s | %s", condition, pert_time, pert_time_units, 
                pert_concentration, pert_concentration_units, cell_line_name, study_accession
            )
        ) %>%
        
        # summarize replicates
        group_by(condition_lab, condition, pert_time, pert_time_units, 
                 pert_concentration, pert_concentration_units, cell_line_name, study_accession,
                 PERT_ENSEMBL, PERT_GENE, regulator) %>%
        summarize(
            activity_bulk = median(activity_bulk, na.rm=TRUE),
            abs_activity_bulk = abs(activity_bulk),
            activity_singlecell = median(activity_singlecell, na.rm=TRUE),
            abs_activity_singlecell = abs(activity_singlecell)
        ) %>%
        ungroup() %>%
        
        # add activity
        group_by(condition_lab) %>%
        mutate(
            total_avail_sfs = sum(regulator %in% unlist(strsplit(PERT_ENSEMBL, ",")))
        ) %>%
        ungroup() %>%
        left_join(cancer_program, by=c("regulator"="ENSEMBL"))
    
    
    protein_activity = protein_activity_bulk %>%
        pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_bulk") %>%
        left_join(
            protein_activity_singlecell %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_singlecell"),
            by = c("PERT_ENSEMBL","regulator")
        ) %>%
        drop_na()
    
   n_cells_threshs = c(1,5,10)
   corrs_act_kd = lapply(n_cells_threshs, function(n_cells_thresh){
       corrs = protein_activity_k562 %>%
                filter(regulator == PERT_ENSEMBL) %>%
                filter(n_cells >= n_cells_thresh) %>%
                drop_na(activity_k562, pert_efficiency_fc) %>%
                group_by(PERT_ENSEMBL, PERT_GENE) %>%
                summarize(
                    correlation_pearson = cor(activity_k562, pert_efficiency_fc, method="pearson"),
                    correlation_spearman = cor(activity_k562, pert_efficiency_fc, method="spearman"),
                    n_conditions = n(),
                    pert_efficiency_fc_std = sd(pert_efficiency_fc),
                    n_cells_thresh = n_cells_thresh
                ) %>%
                ungroup() %>%
                pivot_longer(
                    c(correlation_pearson, correlation_spearman), 
                    names_to="correlation_type", values_to="correlation"
                )
       
       return(corrs)
   }) %>% bind_rows()
    # do the same but randomly
    # correlating bulk vs single cell is not trivial if we consider singlecell batches...
    
    cancer_program_activity = cancer_program %>%
        left_join(
            protein_activity, 
            by=c("ENSEMBL"="regulator")
        ) %>%
        pivot_longer(c(activity_bulk, activity_singlecell), names_to="activity_type", values_to="activity") %>%
        group_by(PERT_ENSEMBL, activity_type, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        drop_na() %>%
        pivot_wider(id_cols=c("PERT_ENSEMBL","driver_type"), names_from="activity_type", values_from="activity")
    
    # plot
    plts = make_plots(protein_activity, cancer_program_activity)
    
    # make figdata
    figdata = make_figdata(protein_activity, cancer_program_activity)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
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

ACTIVITY_TYPES = c(
    "activity_ex", 
    "activity_genexpr", "activity_ew_model_genexpr", "activity_fc_model_genexpr",  
    "activity_scgenexpr", "activity_ew_model_scgenexpr", "activity_fc_model_scgenexpr"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","network_inference")
# protein_activity_ex_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX.tsv.gz")
# protein_activity_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-genexpr.tsv.gz")
# protein_activity_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-scgenexpr.tsv.gz")
# protein_activity_ew_model_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX_from_model_ewlayer_and_genexpr.tsv.gz")
# protein_activity_ew_model_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX_from_model_ewlayer_and_scgenexpr.tsv.gz")
# protein_activity_fc_model_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX_from_model_fclayer_and_genexpr.tsv.gz")
# protein_activity_fc_model_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX_from_model_fclayer_and_scgenexpr.tsv.gz")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")

# PREP_VIPER_DIR = file.path(dirname(ROOT),"publication_viper_splicing","data","prep")
# metadata_file = file.path(PREP_VIPER_DIR,"metadata","tumorigenesis.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","eval_tumorigenesis")

##### FUNCTIONS #####
plot_tumorigenesis = function(protein_activity){
    plts = list()
    
    X = protein_activity %>%
        filter(study_accession=="PRJNA193487") %>%
        drop_na(driver_type) %>%
        group_by(cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(
            activity_ex = median(activity_ex),
            activity_scgenexpr = median(activity_scgenexpr),
            activity_genexpr = median(activity_genexpr),
            activity_ew_model_genexpr = median(activity_ew_model_genexpr),
            activity_ew_model_scgenexpr = median(activity_ew_model_scgenexpr),
            activity_fc_model_genexpr = median(activity_fc_model_genexpr),
            activity_fc_model_scgenexpr = median(activity_fc_model_scgenexpr)
        ) %>%
        ungroup() %>%
        mutate(cell_line_name=factor(
            cell_line_name, levels=c("BJ_PRIMARY","BJ_IMMORTALIZED",
                                     "BJ_TRANSFORMED","BJ_METASTATIC")
        ))
    
    x = X %>%
        pivot_longer(
            c(activity_scgenexpr, activity_genexpr, activity_ex,
              activity_ew_model_genexpr, activity_ew_model_scgenexpr,
              activity_fc_model_genexpr, activity_fc_model_scgenexpr), 
            names_to="activity_type", 
            values_to="activity"
        ) %>%
        mutate(
            activity_type = factor(activity_type, levels=ACTIVITY_TYPES)
        )
    
    plts[["tumorigenesis-cell_line_vs_activity-violin"]] = x %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=activity, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr() + 
        facet_wrap(~activity_type, ncol=4) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Cell Line", y="Protein Activity", fill="Driver Type")

    return(plts)
}


make_plots = function(protein_activity){
    plts = list(
        plot_tumorigenesis(protein_activity)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity){
    figdata = list(
        "eval_tumorigenesis" = list(
            "protein_activity" = protein_activity,
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
    save_plt(plts, "tumorigenesis-cell_line_vs_activity-violin", '.pdf', figs_dir, width=12, height=10)
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
        make_option("--protein_activity_genexpr_file", type="character"),
        make_option("--protein_activity_scgenexpr_file", type="character"),
        make_option("--cancer_program_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    protein_activity_genexpr_file = args[["protein_activity_genexpr_file"]]
    protein_activity_scgenexpr_file = args[["protein_activity_scgenexpr_file"]]
    cancer_program_file = args[["cancer_program_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity_ex = read_tsv(protein_activity_ex_file)
    protein_activity_genexpr = read_tsv(protein_activity_genexpr_file)
    protein_activity_scgenexpr = read_tsv(protein_activity_scgenexpr_file)
    protein_activity_ew_model_genexpr = read_tsv(protein_activity_ew_model_genexpr_file)
    protein_activity_ew_model_scgenexpr = read_tsv(protein_activity_ew_model_scgenexpr_file)
    protein_activity_fc_model_genexpr = read_tsv(protein_activity_fc_model_genexpr_file)
    protein_activity_fc_model_scgenexpr = read_tsv(protein_activity_fc_model_scgenexpr_file)
    metadata = read_tsv(metadata_file)
    cancer_program = read_tsv(cancer_program_file)
    
    # prep
    protein_activity = protein_activity_ex %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity_ex") %>%
        left_join(
            protein_activity_scgenexpr %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_scgenexpr"),
            by = c("sampleID","regulator")
        ) %>%
        left_join(
            protein_activity_genexpr %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_genexpr"),
            by = c("sampleID","regulator")
        ) %>%
        left_join(
            protein_activity_ew_model_genexpr %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_ew_model_genexpr"),
            by = c("sampleID","regulator")
        ) %>%
        left_join(
            protein_activity_ew_model_scgenexpr %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_ew_model_scgenexpr"),
            by = c("sampleID","regulator")
        ) %>%
        left_join(
            protein_activity_fc_model_genexpr %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_fc_model_genexpr"),
            by = c("sampleID","regulator")
        ) %>%
        left_join(
            protein_activity_fc_model_scgenexpr %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_fc_model_scgenexpr"),
            by = c("sampleID","regulator")
        ) %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, activity_scgenexpr, activity_genexpr) %>%
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
            activity_ex = median(activity_ex, na.rm=TRUE),
            activity_genexpr = median(activity_genexpr, na.rm=TRUE),
            activity_scgenexpr = median(activity_scgenexpr, na.rm=TRUE),
            activity_ew_model_genexpr = median(activity_ew_model_genexpr, na.rm=TRUE),
            activity_ew_model_scgenexpr = median(activity_ew_model_scgenexpr, na.rm=TRUE),
            activity_fc_model_genexpr = median(activity_fc_model_genexpr, na.rm=TRUE),
            activity_fc_model_scgenexpr = median(activity_fc_model_scgenexpr, na.rm=TRUE),
        ) %>%
        ungroup() %>%
        
        # add activity
        group_by(condition_lab) %>%
        mutate(
            total_avail_sfs = sum(regulator %in% unlist(strsplit(PERT_ENSEMBL, ",")))
        ) %>%
        ungroup() %>%
        left_join(cancer_program, by=c("regulator"="ENSEMBL"))
    
    # plot
    plts = make_plots(protein_activity)
    
    # make figdata
    figdata = make_figdata(protein_activity)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
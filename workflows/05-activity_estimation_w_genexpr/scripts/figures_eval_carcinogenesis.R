require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables
FIBROBLASTS = c("BJ_PRIMARY","BJ_IMMORTALIZED","BJ_TRANSFORMED","BJ_METASTATIC")
LAB_ORDER = c('WT','C','CB','CBT_228','CBT3','CBTA','CBTP','CBTP3','CBTPA')

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
    'EX',
    'bulkgenexpr',
    'bulkgenexpr-adjusted_ewlayer',
    'bulkgenexpr-adjusted_fclayer',
    'bulkscgenexpr',
    'bulkscgenexpr-adjusted_ewlayer',
    'bulkscgenexpr-adjusted_fclayer',
    'scgenexpr',
    'scgenexpr-adjusted_ewlayer',
    'scgenexpr-adjusted_fclayer'
)
PAL_ACTIVITY_TYPES = setNames(get_palette("Dark2", length(ACTIVITY_TYPES)), ACTIVITY_TYPES)


# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","activity_estimation_w_genexpr")
# protein_activity_bulk_file = file.path(RESULTS_DIR,"files","protein_activity","carcinogenesis-merged.tsv.gz")
# protein_activity_singlecell_file = file.path(RESULTS_DIR,"files","protein_activity","Hodis2022-invitro_eng_melanoc-merged.tsv.gz")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
# metadata_bulk_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets","metadata","tumorigenesis.tsv.gz")
# metadata_singlecell_file = file.path(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-conditions.tsv.gz")
# driver_types_file = file.path(ROOT,"results","new_empirical_network",'files','PANCAN','cancer_program.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","eval_carcinogenesis")

##### FUNCTIONS #####
plot_eval_bulk = function(protein_activity_bulk){
    plts = list()
    
    X = protein_activity_bulk %>%
        drop_na(driver_type) %>%
        group_by(dataset_id, omic_type, model_type, cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        mutate(
            cell_line_name = factor(cell_line_name, levels=FIBROBLASTS),
            dataset_id = gsub("carcinogenesis-","",dataset_id)
        )
    
    plts[["eval_bulk-cell_line_vs_activity_diff-raw-line"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY" & is.na(model_type)) %>%
        mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
        group_by(dataset_id, omic_type, model_type, cell_line_name, study_accession, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(dataset_id, omic_type, model_type, cell_line_name, study_accession) %>%
        summarize(activity_diff = sum(activity)) %>%
        ungroup() %>%
        ggline(
            x="cell_line_name", y="activity_diff", color="dataset_id", numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.05, palette=PAL_ACTIVITY_TYPES
        ) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        labs(x="Cell Line", y="Protein Activity Diff.", color="Network")

    plts[["eval_bulk-cell_line_vs_activity_diff-adjusted-line"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY" & !is.na(model_type) & str_detect(dataset_id,"fclayer")) %>%
        mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
        group_by(dataset_id, omic_type, model_type, cell_line_name, study_accession, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(dataset_id, omic_type, model_type, cell_line_name, study_accession) %>%
        summarize(activity_diff = sum(activity)) %>%
        ungroup() %>%
        ggline(
            x="cell_line_name", y="activity_diff", color="dataset_id", numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.05, palette=PAL_ACTIVITY_TYPES
        ) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        labs(x="Cell Line", y="Protein Activity Diff.", color="Network & Adjustment")
    
    return(plts)
}

plot_eval_singlecell = function(protein_activity_singlecell){
    plts = list()
    
    X = protein_activity_singlecell %>%
        drop_na(driver_type) %>%
        group_by(dataset_id, omic_type, model_type, treatment, driver_type, GENE) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        mutate(
            treatment = factor(treatment, levels=LAB_ORDER),
            dataset_id = gsub("Hodis2022_invitro_eng_melanoc-","",dataset_id)
        )
    
    plts[["eval_singlecell-stage_vs_activity_diff-raw-line"]] = X %>%
        drop_na(driver_type) %>%
        filter(treatment!="WT" & !str_detect(dataset_id,"ewlayer") & is.na(model_type)) %>%
        mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
        group_by(dataset_id, treatment, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(dataset_id, treatment) %>%
        summarize(activity_diff = sum(activity)) %>%
        ungroup() %>%
        ggline(
            x="treatment", y="activity_diff", color="dataset_id", numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.05, palette=PAL_ACTIVITY_TYPES
        ) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        labs(x="Carcinogenic Stage", y="Protein Activity Diff.", color="Network")
    
    plts[["eval_singlecell-stage_vs_activity_diff-adjusted-line"]] = X %>%
        drop_na(driver_type) %>%
        filter(treatment!="WT" & !str_detect(dataset_id,"ewlayer") & !is.na(model_type)) %>%
        mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
        group_by(dataset_id, treatment, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(dataset_id, treatment) %>%
        summarize(activity_diff = sum(activity)) %>%
        ungroup() %>%
        ggline(
            x="treatment", y="activity_diff", color="dataset_id", numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.05, palette=PAL_ACTIVITY_TYPES
        ) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        labs(x="Carcinogenic Stage", y="Protein Activity Diff.", color="Network & Adjustment")
    
    plts[["eval_singlecell-stage_vs_activity_diff-best-line"]] = X %>%
        drop_na(driver_type) %>%
        filter(treatment!="WT" & str_detect(dataset_id,"bulkgenexpr") & !str_detect(dataset_id,"ewlayer") & !is.na(model_type)) %>%
        mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
        group_by(dataset_id, treatment, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(dataset_id, treatment) %>%
        summarize(activity_diff = sum(activity)) %>%
        ungroup() %>%
        ggline(
            x="treatment", y="activity_diff", color="dataset_id", numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.05, palette=PAL_ACTIVITY_TYPES
        ) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        labs(x="Carcinogenic Stage", y="Protein Activity Diff.", color="Network & Adjustment")

    return(plts)
}

make_plots = function(protein_activity_bulk, protein_activity_singlecell){
    plts = list(
        plot_eval_bulk(protein_activity_bulk),
        plot_eval_singlecell(protein_activity_singlecell)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity_bulk, protein_activity_singlecell){
    figdata = list(
        "eval_carcinogenesis" = list(
            "protein_activity_bulk" = protein_activity_bulk,
            "protein_activity_singlecell" = protein_activity_singlecell
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
    save_plt(plts, "eval_bulk-cell_line_vs_activity_diff-raw-line", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "eval_bulk-cell_line_vs_activity_diff-adjusted-line", '.pdf', figs_dir, width=6, height=6)
    
    save_plt(plts, "eval_singlecell-stage_vs_activity_diff-raw-line", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "eval_singlecell-stage_vs_activity_diff-adjusted-line", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "eval_singlecell-stage_vs_activity_diff-best-line", '.pdf', figs_dir, width=7, height=6)
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
        make_option("--metadata_bulk_file", type="character"),
        make_option("--metadata_singlecell_file", type="character"),
        make_option("--driver_types_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    protein_activity_bulk_file = args[["protein_activity_bulk_file"]]
    protein_activity_singlecell_file = args[["protein_activity_singlecell_file"]]
    metadata_bulk_file = args[["metadata_bulk_file"]]
    metadata_singlecell_file = args[["metadata_singlecell_file"]]
    driver_types_file = args[["driver_types_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity_bulk = read_tsv(protein_activity_bulk_file)
    protein_activity_singlecell = read_tsv(protein_activity_singlecell_file)
    metadata_bulk = read_tsv(metadata_bulk_file)
    metadata_singlecell = read_tsv(metadata_singlecell_file)
    driver_types = read_tsv(driver_types_file)
    
    # prep
    protein_activity_bulk = protein_activity_bulk %>%
        left_join(metadata_bulk, by="sampleID") %>%
        drop_na(condition, activity) %>%
        # summarize replicates
        group_by(dataset_id, condition, pert_time, pert_time_units, 
                 pert_concentration, pert_concentration_units, cell_line_name, study_accession,
                 PERT_ENSEMBL, PERT_GENE, regulator) %>%
        summarize(
            activity = median(activity, na.rm=TRUE),
            abs_activity = abs(activity),
        ) %>%
        ungroup() %>%
        separate(dataset_id, sep="-", into=c("dataset_name","omic_type","model_type"), remove=FALSE) %>%
        # add cancer splicing programs
        left_join(driver_types, by=c("regulator"="ENSEMBL"))
    
    protein_activity_singlecell = protein_activity_singlecell %>%
        mutate(dataset_id = gsub("Hodis2022-invitro_eng_melanoc","Hodis2022_invitro_eng_melanoc",dataset_id)) %>%
        left_join(metadata_singlecell, by=c("sampleID"="condition")) %>%
        drop_na(treatment, activity) %>%
        # summarize replicates
        group_by(dataset_id, treatment, cell_type, is_ctl, regulator) %>%
        summarize(
            activity = median(activity, na.rm=TRUE),
            n_cells = sum(n_cells, na.rm=TRUE),
        ) %>%
        ungroup() %>%
        separate(dataset_id, sep="-", into=c("dataset_name","omic_type","model_type"), remove=FALSE) %>%
        left_join(driver_types, by=c("regulator"="ENSEMBL"))
    
    # plot
    plts = make_plots(protein_activity_bulk, protein_activity_singlecell)
    
    # make figdata
    figdata = make_figdata(protein_activity_bulk, protein_activity_singlecell)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
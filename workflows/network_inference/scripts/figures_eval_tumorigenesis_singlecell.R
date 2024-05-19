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

PAL_DARK = "brown"

PAL_DRIVER_TYPE = c(
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

ACTIVITY_TYPES = c(
    "activity_genexpr", "activity_ew_model_genexpr", "activity_fc_model_genexpr",  
    "activity_scgenexpr", "activity_ew_model_scgenexpr", "activity_fc_model_scgenexpr"
)

LAB_ORDER = list(
    "Hodis2022-invitro_eng_melanoc" = c('WT','C','CB','CBT_228','CBT3','CBTA','CBTP','CBTP3','CBTPA'),
    "Becker2021-adenoma" = c("Normal","Unaffected","Polyp","Adenocarcinoma"),
    "Boiarsky2022-myeloma" = c("NBM","MGUS","SMM","MM")
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","network_inference")
# protein_activity_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","Boiarsky2022-myeloma-genexpr.tsv.gz")
# protein_activity_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","Boiarsky2022-myeloma-scgenexpr.tsv.gz")
# protein_activity_ew_model_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","Boiarsky2022-myeloma-EX_from_model_ewlayer_and_genexpr.tsv.gz")
# protein_activity_ew_model_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","Boiarsky2022-myeloma-EX_from_model_ewlayer_and_scgenexpr.tsv.gz")
# protein_activity_fc_model_genexpr_file = file.path(RESULTS_DIR,"files","protein_activity","Boiarsky2022-myeloma-EX_from_model_fclayer_and_genexpr.tsv.gz")
# protein_activity_fc_model_scgenexpr_file = file.path(RESULTS_DIR,"files","protein_activity","Boiarsky2022-myeloma-EX_from_model_fclayer_and_scgenexpr.tsv.gz")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
# metadata_file = file.path(PREP_DIR,"singlecell","Boiarsky2022-myeloma-conditions.tsv.gz")
# dataset = "Boiarsky2022-myeloma"
# figs_dir = file.path(RESULTS_DIR,"figures","eval_tumorigenesis_singlecell-Boiarsky2022-myeloma")

##### FUNCTIONS #####
plot_tumorigenesis = function(protein_activity, dataset){
    plts = list()
    
    X = protein_activity %>%
        drop_na(driver_type) %>%
        group_by(treatment, cell_type, driver_type, GENE) %>%
        summarize(
            activity_scgenexpr = median(activity_scgenexpr),
            activity_genexpr = median(activity_genexpr),
            activity_ew_model_genexpr = median(activity_ew_model_genexpr),
            activity_ew_model_scgenexpr = median(activity_ew_model_scgenexpr),
            activity_fc_model_genexpr = median(activity_fc_model_genexpr),
            activity_fc_model_scgenexpr = median(activity_fc_model_scgenexpr)
        ) %>%
        ungroup() %>%
        mutate(treatment=factor(
            treatment, levels=LAB_ORDER[[dataset]]
        ))
    
    x = X %>%
        pivot_longer(
            c(activity_scgenexpr, activity_genexpr,
              activity_ew_model_genexpr, activity_ew_model_scgenexpr,
              activity_fc_model_genexpr, activity_fc_model_scgenexpr), 
            names_to="activity_type", 
            values_to="activity"
        ) %>%
        mutate(
            activity_type = factor(activity_type, levels=ACTIVITY_TYPES)
        )
    
    plts[["tumorigenesis-treatment_vs_activity-violin"]] = x %>%
        ggplot(aes(x=treatment, y=activity, group=interaction(treatment,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 45) + 
        facet_wrap(~activity_type+cell_type, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Treatment", y="Protein Activity", fill="Driver Type")
    
    x = x %>%
        group_by(treatment, driver_type, activity_type, cell_type) %>%
        summarize(
            activity = median(activity)
        ) %>%
        ungroup() %>%
        pivot_wider(id_cols=c("cell_type","treatment","activity_type"), names_from="driver_type", values_from="activity") %>%
        mutate(program_fc = `Oncogenic` - `Tumor suppressor`)
    
    plts[["tumorigenesis-treatment_vs_program_fc-line"]] = x %>%
        ggline(
            x="treatment", y="program_fc", color="activity_type", 
            palette="Paired", size=LINE_SIZE, point.size=0.05
        ) +
        theme_pubr(x.text.angle = 45) + 
        facet_wrap(~cell_type, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Treatment", y="Oncogenic vs Tumor Suppressor Activity", color="Activity Type")

    return(plts)
}


make_plots = function(protein_activity, dataset){
    plts = list(
        plot_tumorigenesis(protein_activity, dataset)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity){
    figdata = list(
        "eval_tumorigenesis_singlecell" = list(
            "protein_activity" = protein_activity
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


save_plots = function(plts, figs_dir, dataset){
    
    if (dataset=="Hodis2022-invitro_eng_melanoc"){
        
        save_plt(plts, "tumorigenesis-treatment_vs_activity-violin", '.pdf', figs_dir, width=12, height=10)
        save_plt(plts, "tumorigenesis-treatment_vs_program_fc-line", '.pdf', figs_dir, width=5, height=7)
        
    } else if (dataset=="Becker2021-adenoma"){
        
        save_plt(plts, "tumorigenesis-treatment_vs_activity-violin", '.pdf', figs_dir, width=30, height=30)
        save_plt(plts, "tumorigenesis-treatment_vs_program_fc-line", '.pdf', figs_dir, width=12, height=12)   

    } else if (dataset=="Boiarsky2022-myeloma"){
        
        save_plt(plts, "tumorigenesis-treatment_vs_activity-violin", '.pdf', figs_dir, width=30, height=30)
        save_plt(plts, "tumorigenesis-treatment_vs_program_fc-line", '.pdf', figs_dir, width=12, height=7)

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
        make_option("--protein_activity_genexpr_file", type="character"),
        make_option("--protein_activity_scgenexpr_file", type="character"),
        make_option("--protein_activity_ew_model_genexpr_file", type="character"),
        make_option("--protein_activity_ew_model_scgenexpr_file", type="character"),
        make_option("--protein_activity_fc_model_genexpr_file", type="character"),
        make_option("--protein_activity_fc_model_scgenexpr_file", type="character"),
        make_option("--cancer_program_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--dataset", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    protein_activity_genexpr_file = args[["protein_activity_genexpr_file"]]
    protein_activity_scgenexpr_file = args[["protein_activity_scgenexpr_file"]]
    protein_activity_ew_model_genexpr_file = args[["protein_activity_ew_model_genexpr_file"]]
    protein_activity_ew_model_scgenexpr_file = args[["protein_activity_ew_model_scgenexpr_file"]]
    protein_activity_fc_model_genexpr_file = args[["protein_activity_fc_model_genexpr_file"]]
    protein_activity_fc_model_scgenexpr_file = args[["protein_activity_fc_model_scgenexpr_file"]]
    metadata_file = args[["metadata_file"]]
    cancer_program_file = args[["cancer_program_file"]]
    dataset = args[["dataset"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity_genexpr = read_tsv(protein_activity_genexpr_file)
    protein_activity_scgenexpr = read_tsv(protein_activity_scgenexpr_file)
    protein_activity_ew_model_genexpr = read_tsv(protein_activity_ew_model_genexpr_file)
    protein_activity_ew_model_scgenexpr = read_tsv(protein_activity_ew_model_scgenexpr_file)
    protein_activity_fc_model_genexpr = read_tsv(protein_activity_fc_model_genexpr_file)
    protein_activity_fc_model_scgenexpr = read_tsv(protein_activity_fc_model_scgenexpr_file)
    metadata = read_tsv(metadata_file)
    cancer_program = read_tsv(cancer_program_file)
    
    # prep
    protein_activity = protein_activity_scgenexpr %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity_scgenexpr") %>%
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
        left_join(metadata, by=c("sampleID"="condition")) %>%
        drop_na(treatment, activity_scgenexpr, activity_genexpr) %>%
        # summarize replicates
        group_by(treatment, cell_type, is_ctl, regulator) %>%
        summarize(
            activity_genexpr = median(activity_genexpr, na.rm=TRUE),
            activity_scgenexpr = median(activity_scgenexpr, na.rm=TRUE),
            activity_ew_model_genexpr = median(activity_ew_model_genexpr, na.rm=TRUE),
            activity_ew_model_scgenexpr = median(activity_ew_model_scgenexpr, na.rm=TRUE),
            activity_fc_model_genexpr = median(activity_fc_model_genexpr, na.rm=TRUE),
            activity_fc_model_scgenexpr = median(activity_fc_model_scgenexpr, na.rm=TRUE),
            n_cells = sum(n_cells, na.rm=TRUE),
        ) %>%
        ungroup() %>%

        mutate(
            condition_lab = sprintf(
                "%s (n=%s)", treatment, n_cells
            )
        ) %>%
        
        # add activity
        group_by(condition_lab) %>%
        mutate(
            total_avail_sfs = n()
        ) %>%
        ungroup() %>%
        left_join(cancer_program, by=c("regulator"="ENSEMBL"))
    
    # plot
    plts = make_plots(protein_activity, dataset)
    
    # make figdata
    figdata = make_figdata(protein_activity)
    
    # save
    save_plots(plts, figs_dir, dataset)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

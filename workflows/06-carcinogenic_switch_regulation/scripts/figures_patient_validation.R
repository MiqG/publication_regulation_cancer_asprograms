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
require(ggrepel)
require(clusterProfiler)
require(scattermore)
require(ComplexHeatmap)
require(ggplotify)
require(ggvenn)

# variables
LABS_BULK = c("BJ_PRIMARY","BJ_IMMORTALIZED","BJ_TRANSFORMED","BJ_METASTATIC")
LABS_SC = c('WT','C','CB','CBT_228','CBT3','CBTA','CBTP','CBTP3','CBTPA')

RANDOM_SEED = 1234

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_EVAL_TYPE = c(
    "random" = "lightgrey",
    "real" = "orange"
)
PAL_DARK = "brown"
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

PAL_DRIVER_TYPE = c(
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

PAL_GENE_TYPE = c(
    "Not SF"="darkgreen",
    "Non-driver SF"="darkred",
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)
PAL_DATASETS = c(
    "Danielsson2013-fibroblasts"="#1729AD",
    "Hodis2022-invitro_eng_melanoc"="#FE5F55",
    "ReplogleWeissman2022_rpe1"="#12664F"
)

LABS_ORDER = c("Normal","Unaffected","Polyp","Adenocarcinoma")

GENES_OI = c(
    "ENSG00000122565", # CBX3
    "ENSG00000088205" # DDX18
)

HALLMARKS_OI = c("HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2")

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","program_regulation")
# NETWORKS_DIR = file.path(ROOT,"results","network_inference")

# carcinogenesis_genexpr_file = file.path(PREP_DIR,"singlecell","Becker2021-adenoma-pseudobulk.tsv.gz")
# carcinogenesis_activity_file = file.path(NETWORKS_DIR,"files","protein_activity","Becker2021-adenoma-EX_from_model_fclayer_and_scgenexpr.tsv.gz")
# carcinogenesis_hallmarks_file = file.path(NETWORKS_DIR,"files","gsea","Becker2021-adenoma-hallmarks.tsv.gz")
# carcinogenesis_metadata_file = file.path(PREP_DIR,"singlecell","Becker2021-adenoma-conditions.tsv.gz")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,"figures","patient_validation")

##### FUNCTIONS #####
plot_carcinogenesis = function(genexpr, activity, hallmarks){
    plts = list()
    
    plts[["carcinogenesis-activity-box"]] = activity %>%
        filter(cell_type=="Stem") %>%
        drop_na(driver_type) %>%
        group_by(treatment, cell_type, replicate, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        ggstripchart(x="treatment", y="activity", color="driver_type", position=position_jitterdodge(0.2, dodge.width=0.9)) +
        geom_boxplot(aes(color=driver_type), fill=NA, outlier.shape=NA) + 
        color_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle=45) +
        facet_wrap(~cell_type) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Carcinogenic Stage", y="Splicing Program Activity")
    
    plts[["carcinogenesis-activity_diff-box"]] = activity %>%
        filter(cell_type=="Stem") %>%
        drop_na(driver_type) %>%
        group_by(treatment, cell_type, replicate, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        pivot_wider(names_from="driver_type", values_from="activity") %>%
        mutate(activity_diff=`Oncogenic` - `Tumor suppressor`) %>%
        ggstripchart(x="treatment", y="activity_diff", color=PAL_DARK, size=1, position=position_jitter(0.1)) +
        geom_boxplot(width=0.2, color=PAL_DARK, fill=NA, outlier.shape=NA) +
        geom_hline(yintercept=0, size=LINE_SIZE , linetype="dashed", color="black") +
        geom_text(
            aes(y=-0.6, label=label),
            . %>% count(treatment) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        theme_pubr(x.text.angle=45) +
        facet_wrap(~cell_type) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Carcinogenic Stage", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")

    
    plts[["carcinogenesis-hallmarks-box"]] = hallmarks %>%
        filter(cell_type=="Stem" & Description%in%HALLMARKS_OI) %>%
        ggstripchart(x="treatment", y="NES", color="Description", size=1, position=position_jitter(0.1)) +
        geom_boxplot(aes(color=Description), width=0.2, fill=NA, outlier.shape=NA) +
        geom_text(
            aes(y=-2, label=label),
            . %>% count(treatment) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        color_palette(c("blue","lightblue")) +
        theme_pubr(x.text.angle=45) +
        facet_wrap(~cell_type+Description, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Carcinogenic Stage", y="NES")
    
    plts[["carcinogenesis-genexpr-box"]] = genexpr %>%
        filter(ENSEMBL%in%GENES_OI & cell_type=="Stem") %>%
        ggstripchart(x="treatment", y="genexpr", color="ENSEMBL", position=position_jitterdodge(0.2, dodge.width=0.9)) +
        geom_boxplot(aes(color=ENSEMBL), fill=NA, outlier.shape=NA) + 
        theme_pubr(x.text.angle=45) +
        facet_wrap(~cell_type) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Carcinogenic Stage", y="log2(TPM + 1)")
    
    terms_oi = corrs %>%
        group_by(Description) %>%
        summarize(avg_correlation = median(correlation_pearson)) %>%
        ungroup() 
    terms_oi = terms_oi %>%
        slice_max(avg_correlation, n=5) %>%
        bind_rows(
            terms_oi %>% slice_min(avg_correlation, n=5)
        ) %>%
        pull(Description)
    
    plts[["carcinogenesis-corrs-activity_diff_vs_nes-box"]] = corrs %>%
        filter(Description%in%terms_oi) %>%
        mutate(Description = fct_reorder(Description, correlation_pearson)) %>%
        ggstripchart(x="Description", y="correlation_pearson", color="cell_type") +
        geom_boxplot(color="black", fill=NA, outlier.shape=NA) + 
        color_palette("Paired") +
        labs(x="Top Correlating Hallmark", y="Pearson Correlation", color="Cell Type") +
        coord_flip()
    
    return(plts)
}


make_plots = function(genexpr, activity, hallmarks){
    plts = list(
        plot_carcinogenesis(genexpr, activity, hallmarks)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(corrs){
    figdata = list(
        "gsea_carcinogenesis" = list(
            "correlation_analysis" = corrs
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
    save_plt(plts, "corrs-nes_vs_activity_diff-violin", '.pdf', figs_dir, width=4, height=4)
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
        make_option("--carcinogenesis_genexpr_file", type="character"),
        make_option("--carcinogenesis_activity_file", type="character"),
        make_option("--carcinogenesis_hallmarks_file", type="character"),
        make_option("--carcinogenesis_metadata_file", type="character"),
        make_option("--cancer_program_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    carcinogenesis_genexpr_file = args[["carcinogenesis_genexpr_file"]]
    carcinogenesis_activity_file = args[["carcinogenesis_activity_file"]]
    carcinogenesis_hallmarks_file = args[["carcinogenesis_hallmarks_file"]]
    carcinogenesis_metadata_file = args[["carcinogenesis_metadata_file"]]
    cancer_program_file = args[["cancer_program_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    set.seed(RANDOM_SEED)
    
    # load
    carcinogenesis_genexpr = read_tsv(carcinogenesis_genexpr_file)
    carcinogenesis_activity = read_tsv(carcinogenesis_activity_file)
    carcinogenesis_hallmarks = read_tsv(carcinogenesis_hallmarks_file)
    carcinogenesis_metadata = read_tsv(carcinogenesis_metadata_file)
    cancer_program = read_tsv(cancer_program_file)
    
    # prep
    carcinogenesis_metadata = carcinogenesis_metadata %>%
        mutate(treatment = factor(treatment, levels=LABS_ORDER))
    
    genexpr = carcinogenesis_genexpr %>%
        pivot_longer(-ENSEMBL, names_to="sampleID", values_to="genexpr") %>%
        left_join(carcinogenesis_metadata, by=c("sampleID"="condition"))
    
    activity = carcinogenesis_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(carcinogenesis_metadata, by=c("sampleID"="condition")) %>%
        left_join(cancer_program %>% distinct(ENSEMBL, GENE, driver_type), by=c("regulator"="ENSEMBL"))
    
    hallmarks = carcinogenesis_hallmarks %>%
        left_join(carcinogenesis_metadata, by=c("sampleID"="condition"))    
    
    # correlations hallmarks NES vs Cancer Splicing Program activity difference
    corrs = activity %>%
        #filter(cell_type=="Stem") %>%
        drop_na(driver_type) %>%
        group_by(treatment, cell_type, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        pivot_wider(names_from="driver_type", values_from="activity") %>%
        mutate(activity_diff=`Oncogenic` - `Tumor suppressor`) %>%
        left_join(
            hallmarks %>%
                group_by(Description, treatment, cell_type) %>%
                summarize(NES = median(NES)) %>%
                ungroup(), 
            by=c("treatment","cell_type")
        ) %>%
        group_by(Description, cell_type) %>%
        summarize(
            correlation_pearson = cor(NES, activity_diff, method="pearson", use="pairwise.complete.obs"),
            n_obs = n()
        ) %>%
        ungroup() %>%
        arrange(correlation_pearson)
    
    # plot
    plts = make_plots(genexpr, activity, hallmarks)
    
    # make figdata
    figdata = make_figdata(corrs)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}

##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

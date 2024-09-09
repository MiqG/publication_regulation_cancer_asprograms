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

# variables
LABS_BULK = c("BJ_PRIMARY","BJ_IMMORTALIZED","BJ_TRANSFORMED","BJ_METASTATIC")
LABS_SC = c('WT','C','CB','CBT_228','CBT3','CBTA','CBTP','CBTP3','CBTPA')

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
    "Oncogenic"="#F6AE2D",
    "Tumor suppressor"="#6C98B3"
)
PAL_DATASETS = c(
    "Danielsson2013-fibroblasts"="#1729AD",
    "Hodis2022-invitro_eng_melanoc"="#FE5F55",
    "ReplogleWeissman2022_rpe1"="#12664F"
)
# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# NETWORKS_DIR = file.path(ROOT,"results","network_inference") 
# RESULTS_DIR = file.path(ROOT,"results","program_regulation")

# carcinogenesis_bulk_activity_file = file.path(NETWORKS_DIR,"figures","eval_tumorigenesis","figdata","eval_tumorigenesis","protein_activity.tsv.gz")
# carcinogenesis_bulk_hallmarks_file = file.path(NETWORKS_DIR,"files","gsea","tumorigenesis-genexpr-hallmarks.tsv.gz")

# carcinogenesis_singlecell_activity_file = file.path(NETWORKS_DIR,"figures","eval_tumorigenesis_singlecell-Hodis2022-invitro_eng_melanoc","figdata","eval_tumorigenesis_singlecell","protein_activity.tsv.gz")
# carcinogenesis_singlecell_hallmarks_file = file.path(NETWORKS_DIR,"files","gsea","Hodis2022-invitro_eng_melanoc-hallmarks.tsv.gz")

# pertseq_activity_file = file.path(RESULTS_DIR,"figures","upstream_regulators","figdata","upstream_regulators","cancer_program_activity.tsv.gz")
# pertseq_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","ReplogleWeissman2022_rpe1-hallmarks.tsv.gz")

# PREP_VIPER_DIR = file.path(dirname(ROOT),"publication_viper_splicing","data","prep")
# carcinogenesis_bulk_metadata_file = file.path(PREP_VIPER_DIR,"metadata","tumorigenesis.tsv.gz")

# carcinogenesis_singlecell_metadata_file = file.path(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-conditions.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,"figures","gsea_carcinogenesis")

##### FUNCTIONS #####
plot_corrs = function(corrs, experiments){
    plts = list()
    
    X = corrs

    # distributions of correlations
    plts[["corrs-nes_vs_activity_diff-violin"]] = X %>%
        ggviolin(x="dataset", y="correlation_diff_activity", fill="dataset", 
                 color=NA, palette=PAL_DATASETS, trim=TRUE) +
        geom_text(
            aes(y=-1.05, label=label),
            . %>% count(dataset) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        facet_wrap(~dataset, scales="free_x") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Dataset", y="Pearson Correlation")
    
    # line plots top correlations during carcinogenesis
    ## bulk
    pathways_oi = X %>%
        filter(dataset=="Danielsson2013-fibroblasts") %>%
        slice_max(correlation_diff_activity, n=5)
    x = experiments %>%
        filter(dataset=="Danielsson2013-fibroblasts") %>%
        filter(Description %in% pathways_oi[["Description"]]) %>%
        pivot_longer(cols=c("NES","activity_diff")) %>%
        distinct(cell_line_name, Description, value, name) %>%
        mutate(Description = ifelse(name=="activity_diff", "switch", Description)) %>%
        distinct(cell_line_name, Description, value, name) %>%
        left_join(pathways_oi %>% distinct(Description, correlation_diff_activity), by="Description") %>%
        mutate(
            correlation_diff_activity = ifelse(name=="activity_diff", 1, correlation_diff_activity),
            corr_class = cut(correlation_diff_activity, breaks = 10)
        )
    plts[["corrs-stage_vs_act_diff_and_enrichment-bulk-line"]] = x %>%
        ggplot(aes(x=cell_line_name, y=value, group=Description)) +
        geom_line(aes(color=corr_class)) +
        geom_text_repel(
            aes(label=Description, color=corr_class),
            . %>% group_by(Description) %>% slice_max(cell_line_name),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1,
            max.overlaps=50
        ) +
        theme_pubr() +
        facet_wrap(~name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        color_palette(get_palette(c("grey","#1729AD"), 4)) +
        labs(x="Carcinogenic Stage", y="NES or DeltaActivity", color="Pearson Correlation")
    
    ## single cell
    pathways_oi = X %>%
        filter(dataset=="Hodis2022-invitro_eng_melanoc") %>%
        slice_max(correlation_diff_activity, n=5)
    x = experiments %>%
        filter(dataset=="Hodis2022-invitro_eng_melanoc") %>%
        filter(Description %in% pathways_oi[["Description"]]) %>%
        pivot_longer(cols=c("NES","activity_diff")) %>%
        distinct(treatment, Description, value, name) %>%
        mutate(Description = ifelse(name=="activity_diff", "switch", Description)) %>%
        distinct(treatment, Description, value, name) %>%
        left_join(pathways_oi %>% distinct(Description, correlation_diff_activity), by="Description") %>%
        mutate(
            correlation_diff_activity = ifelse(name=="activity_diff", 1, correlation_diff_activity),
            corr_class = cut(correlation_diff_activity, breaks = 10)
        )
    plts[["corrs-stage_vs_act_diff_and_enrichment-singlecell-line"]] = x %>%
        ggplot(aes(x=treatment, y=value, group=Description)) +
        geom_line(aes(color=corr_class)) +
        geom_text_repel(
            aes(label=Description, color=corr_class),
            . %>% group_by(Description) %>% slice_max(treatment),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1,
            max.overlaps=50
        ) +
        theme_pubr() +
        facet_wrap(~name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        color_palette(get_palette(c("grey","#FE5F55"), 3)) +
        labs(x="Carcinogenic Stage", y="NES or DeltaActivity", color="Pearson Correlation")
           
    # top correlations with splicing switch
    ## bulk
    plts[["corrs-nes_vs_activity_diff-bulk-bar"]] = X %>%
        filter(dataset=="Danielsson2013-fibroblasts") %>%
        slice_max(correlation_diff_activity, n=10) %>%
        ggbarplot(x="Description", y="correlation_diff_activity", fill="correlation_diff_activity", color=NA) +
        labs(x="MSigDB Hallmark", y="Pearson Correlation") +
        guides(fill="none") +
        coord_flip() +
        scale_x_discrete(limits=rev) +
        scale_fill_gradient(low="black", high="#1729AD")
    
    ## singlecell
    plts[["corrs-nes_vs_activity_diff-siglecell-bar"]] = X %>%
        filter(dataset=="Hodis2022-invitro_eng_melanoc") %>%
        slice_max(correlation_diff_activity, n=10) %>%
        ggbarplot(x="Description", y="correlation_diff_activity", fill="correlation_diff_activity", color=NA) +
        labs(x="MSigDB Hallmark", y="Pearson Correlation") +
        guides(fill="none") +
        coord_flip() +
        scale_x_discrete(limits=rev) +
        scale_fill_gradient(low="black", high="#FE5F55")
    
    ## pert seq
    plts[["corrs-nes_vs_activity_diff-pertseq-bar"]] = X %>%
        filter(dataset=="ReplogleWeissman2022_rpe1") %>%
        slice_max(correlation_diff_activity, n=5) %>%
        ggbarplot(x="Description", y="correlation_diff_activity", fill="correlation_diff_activity", color=NA) +
        labs(x="MSigDB Hallmark", y="Pearson Correlation") +
        guides(fill="none") +
        coord_flip() +
        scale_x_discrete(limits=rev) +
        scale_fill_gradient(low="grey", high="#12664F")
    
    # relationship delta activity with NES
    pathways_oi = X %>%
        filter(dataset=="ReplogleWeissman2022_rpe1") %>%
        slice_max(correlation_diff_activity, n=5)
    x = experiments %>%
        filter(dataset=="ReplogleWeissman2022_rpe1") %>%
        filter(Description %in% pathways_oi[["Description"]]) %>% 
        left_join(pathways_oi %>% distinct(Description, correlation_diff_activity), by="Description")
    
    plts[["corrs-nes_vs_activity_diff-pertseq-scatter"]] = x %>%
        mutate(Description = fct_reorder(Description, -correlation_diff_activity)) %>%
        ggscatter(x="NES", y="activity_diff", color="correlation_diff_activity", alpha=0.2) +
        facet_wrap(~Description, nrow=1) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(color="none") +
        scale_color_gradient(low="grey", high="#12664F")
    
    
    
    
    
    
    
    
    %>%
        pivot_longer(cols=c("NES","activity_diff")) %>%
        distinct(treatment, Description, value, name) %>%
        mutate(Description = ifelse(name=="activity_diff", "switch", Description)) %>%
        distinct(treatment, Description, value, name) %>%
        left_join(pathways_oi %>% distinct(Description, correlation_diff_activity), by="Description") %>%
        mutate(
            correlation_diff_activity = ifelse(name=="activity_diff", 1, correlation_diff_activity),
            corr_class = cut(correlation_diff_activity, breaks = 10)
        )
        
    
    return(plts)
}


make_plots = function(corrs, experiments){
    plts = list(
        plot_corrs(corrs, experiments)
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
    save_plt(plts, "corrs-nes_vs_activity_diff-bulk", '.pdf', figs_dir, width=7.5, height=4)
    save_plt(plts, "corrs-nes_vs_activity_diff-singlecell", '.pdf', figs_dir, width=7.5, height=4)
    save_plt(plts, "corrs-nes_vs_activity_diff-pertseq", '.pdf', figs_dir, width=7.5, height=4)
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
        make_option("--carcinogenesis_bulk_activity_file", type="character"),
        make_option("--carcinogenesis_bulk_hallmarks_file", type="character"),
        make_option("--carcinogenesis_singlecell_activity_file", type="character"),
        make_option("--carcinogenesis_singlecell_hallmarks_file", type="character"),
        make_option("--pertseq_activity_file", type="character"),
        make_option("--pertseq_hallmarks_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    carcinogenesis_bulk_activity_file = args[["carcinogenesis_bulk_activity_file"]]
    carcinogenesis_bulk_metadata_file = args[["carcinogenesis_bulk_metadata_file"]]
    carcinogenesis_bulk_hallmarks_file = args[["carcinogenesis_bulk_hallmarks_file"]]
    carcinogenesis_singlecell_activity_file = args[["carcinogenesis_singlecell_activity_file"]]
    carcinogenesis_singlecell_hallmarks_file = args[["carcinogenesis_singlecell_hallmarks_file"]]
    pertseq_activity_file = args[["pertseq_activity_file"]]
    pertseq_hallmarks_file = args[["pertseq_hallmarks_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    carcinogenesis_bulk_activity = read_tsv(carcinogenesis_bulk_activity_file)
    carcinogenesis_bulk_metadata = read_tsv(carcinogenesis_bulk_metadata_file)
    carcinogenesis_bulk_hallmarks = read_tsv(carcinogenesis_bulk_hallmarks_file)
    carcinogenesis_singlecell_activity = read_tsv(carcinogenesis_singlecell_activity_file)
    carcinogenesis_singlecell_metadata = read_tsv(carcinogenesis_singlecell_metadata_file)
    carcinogenesis_singlecell_hallmarks = read_tsv(carcinogenesis_singlecell_hallmarks_file)
    pertseq_activity = read_tsv(pertseq_activity_file)
    pertseq_hallmarks = read_tsv(pertseq_hallmarks_file)
    
    # prep: combine activity and hallmarks
    ## bulk carcinogenesis
    carcinogenesis_bulk = carcinogenesis_bulk_hallmarks %>%
        left_join(carcinogenesis_bulk_metadata, by="sampleID") %>%
        filter(cell_line_name %in% LABS_BULK) %>%
        drop_na(NES) %>%
        group_by(Description, cell_line_name) %>%
        summarize(
            NES = mean(NES),
        ) %>%
        ungroup() %>%
        left_join(
            carcinogenesis_bulk_activity %>%
                group_by(cell_line_name, driver_type) %>%
                summarize(
                    activity = median(activity_fc_model_scgenexpr)
                ) %>%
                ungroup() %>%
                drop_na(driver_type) %>%
                pivot_wider(id_cols="cell_line_name", names_from="driver_type", values_from="activity") %>%
                mutate(activity_diff = `Oncogenic` - `Tumor suppressor`),
            by=c("cell_line_name")
        ) %>%
        mutate(
            cell_line_name = factor(cell_line_name, levels=LABS_BULK),
            dataset = "Danielsson2013-fibroblasts"
        )
    
    ## single cell carcinogenesis
    carcinogenesis_singlecell = carcinogenesis_singlecell_hallmarks %>%
        left_join(carcinogenesis_singlecell_metadata, by=c("sampleID"="condition")) %>%
        filter(treatment %in% LABS_SC) %>%
        drop_na(NES) %>%
        group_by(Description, treatment, cell_type) %>%
        summarize(
            NES = mean(NES),
        ) %>%
        left_join(
            carcinogenesis_singlecell_activity %>%
                group_by(treatment, driver_type, cell_type) %>%
                summarize(
                    activity = median(activity_fc_model_scgenexpr)
                ) %>%
                ungroup() %>%
                drop_na(driver_type) %>%
                pivot_wider(id_cols=c("cell_type","treatment"), names_from="driver_type", values_from="activity") %>%
                mutate(activity_diff = `Oncogenic` - `Tumor suppressor`),
            by=c("treatment","cell_type")
        ) %>%
        mutate(
            treatment = factor(treatment, levels=LABS_SC),
            dataset = "Hodis2022-invitro_eng_melanoc"
        )
    
    
    ## Perturb seq
    pertseq = pertseq_hallmarks %>%
        drop_na(NES) %>%
        left_join(
            pertseq_activity %>% filter(activity_type=="activity_rpe1"), 
            by=c("sampleID"="PERT_ENSEMBL")
        ) %>%
        mutate(dataset = "ReplogleWeissman2022_rpe1")
    
    # compute correlations between switch activation and hallmark enrichments
    ## bulk
    corrs_bulk = carcinogenesis_bulk %>%
        group_by(dataset, Description) %>%
        summarize(
            correlation_stage = cor(NES, as.numeric(cell_line_name), method="spearman"),
            correlation_diff_activity = cor(NES, activity_diff, method="pearson"),
            n_obs = n()
        ) %>%
        ungroup()
    
    ## single cell
    corrs_singlecell = carcinogenesis_singlecell %>%
        group_by(dataset, Description) %>%
        summarize(
            correlation_stage = cor(NES, as.numeric(treatment), method="spearman"),
            correlation_diff_activity = cor(NES, activity_diff, method="pearson"),
            n_obs = n()
        ) %>%
        ungroup()
    
    ## perturb seq
    corrs_pertseq = pertseq %>%
        group_by(dataset, Description) %>%
        summarize(
            correlation_diff_activity = cor(NES, activity_diff, method="pearson"),
            n_obs = n()
        ) %>%
        ungroup()
    
    # merge
    experiments = bind_rows(carcinogenesis_bulk, carcinogenesis_singlecell, pertseq)
    corrs = bind_rows(corrs_bulk, corrs_singlecell, corrs_pertseq)
    
    # plot
    plts = make_plots(corrs, experiments)
    
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

#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Notes
# --------------
# - MYC turns the carcinogenic splicing switch on
# - MYC upregulates DDX18 transcription
# - overexpression of DDX18 decreases its activity
#
# TODO
# ----
# - plot everything as fold changes w.r.t. MCF10A
# - Check gene expression/activity changes of all splicing factors upon MYC activation that are known targets of MYC (CHEA database) --> volcano plot --> we are interested in those that both change their expression at the same time that the switch is turned on 
#     - consider whether they are annotated as suppressor or oncogenic?
#     - correlate their FC gene expression or activity with the activation of the switch
#     - we are interested in those whose change in expression (and activity) correlates with the switch turning on or not upon knockdown in the perturb seq (if available)

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
# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# PREP_VIPER_DIR = file.path(dirname(ROOT),"publication_viper_splicing","data","prep")
# NETWORKS_DIR = file.path(ROOT,"results","network_inference") 
# RESULTS_DIR = file.path(ROOT,"results","program_regulation")

# carcinogenesis_bulk_genexpr_file = file.path(PREP_VIPER_DIR,'genexpr_tpm','tumorigenesis.tsv.gz')
# carcinogenesis_bulk_activity_file = file.path(NETWORKS_DIR,"figures","eval_tumorigenesis","figdata","eval_tumorigenesis","protein_activity.tsv.gz")
# carcinogenesis_bulk_hallmarks_file = file.path(NETWORKS_DIR,"files","gsea","tumorigenesis-genexpr-hallmarks.tsv.gz")

# carcinogenesis_singlecell_genexpr_file = file.path(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-pseudobulk.tsv.gz")
# carcinogenesis_singlecell_activity_file = file.path(NETWORKS_DIR,"figures","eval_tumorigenesis_singlecell-Hodis2022-invitro_eng_melanoc","figdata","eval_tumorigenesis_singlecell","protein_activity.tsv.gz")
# carcinogenesis_singlecell_hallmarks_file = file.path(NETWORKS_DIR,"files","gsea","Hodis2022-invitro_eng_melanoc-hallmarks.tsv.gz")
# pertseq_activity_file = file.path(RESULTS_DIR,"figures","upstream_regulators","figdata","upstream_regulators","cancer_program_activity.tsv.gz")
# pertseq_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","ReplogleWeissman2022_rpe1-hallmarks.tsv.gz")
# PREP_VIPER_DIR = file.path(dirname(ROOT),"publication_viper_splicing","data","prep")
# carcinogenesis_bulk_metadata_file = file.path(PREP_VIPER_DIR,"metadata","tumorigenesis.tsv.gz")
# carcinogenesis_singlecell_metadata_file = file.path(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-conditions.tsv.gz")
# msigdb_dir = file.path(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")

# urbanski_metadata_file = file.path(PREP_DIR,"metadata","Urbanski2022.tsv.gz")
# urbanski_genexpr_file = file.path(PREP_DIR,'genexpr_tpm',"Urbanski2022.tsv.gz")
# urbanski_ex_file = file.path(PREP_DIR,'event_psi',"Urbanski2022-EX.tsv.gz")
# urbanski_activiy_file = file.path(RESULTS_DIR,"files","protein_activity","Urbanski2022-EX.tsv.gz")
# urbanski_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","Urbanski2022-hallmarks.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,"figures","gsea_carcinogenesis")

##### FUNCTIONS #####
load_ontologies = function(msigdb_dir){
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt"))
    )
    return(ontologies)
}

plot_corrs = function(corrs, experiments, ontologies, cancer_program){
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
    
    # bulk vs single cell
    plts[["corrs-bulk_vs_singlecell-scatter"]] = X %>%
        filter(dataset %in% c("Danielsson2013-fibroblasts","Hodis2022-invitro_eng_melanoc")) %>%
        pivot_wider(id_cols="Description", names_from="dataset", values_from="correlation_diff_activity") %>%
        rowwise() %>%
        mutate(avg_corr = mean(c(`Danielsson2013-fibroblasts`,`Hodis2022-invitro_eng_melanoc`)) ) %>%
        ungroup() %>%
        ggscatter(
            x="Danielsson2013-fibroblasts", y="Hodis2022-invitro_eng_melanoc", 
            alpha=0.5, size=1, color=PAL_DARK
        ) +
        # geom_text_repel(
        #     aes(label=Description),
        #     . %>% slice_max(abs(avg_corr), n=10),
        #     size=FONT_SIZE+2, family=FONT_FAMILY, segment.size=0.1,
        #     max.overlaps=50
        # ) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_abline(intercept=0, slope=1, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        theme(aspect.ratio=1)
    
    x = X %>%
        group_by(Description) %>%
        filter(dataset %in% c("Danielsson2013-fibroblasts","Hodis2022-invitro_eng_melanoc")) %>%
        mutate(avg_corr = mean(correlation_diff_activity)) %>%
        ungroup() %>%
        slice_max(abs(avg_corr), n=20) %>%
        mutate(Description = fct_reorder(Description, avg_corr))
    
    plts[["corrs-bulk_vs_singlecell-bar"]] = x %>%
        ggbarplot(x="Description", y="correlation_diff_activity", color=NA, fill="dataset", 
                  palette=PAL_DATASETS, position=position_dodge(0.9)) +
        labs(x="Top Correlating Hallmark", y="Pearson Correlation", fill="") +
        coord_flip()
    
    # Upset plot showing the overlaps between shortlisted pathways
    onto_list = ontologies[["hallmarks"]] %>%
        filter(term %in% x[["Description"]]) %>%
        mutate(term = as.character(term)) %>%
        with(., split(gene, term))
    m = onto_list %>% 
        list_to_matrix() %>% 
        make_comb_mat()
    plts[["corrs-hallmark_sets_overlaps-upset"]] = m %>%
        UpSet(comb_order = order(comb_size(m)), 
              comb_col = PAL_DARK,
              top_annotation = upset_top_annotation(m, gp = gpar(fill = PAL_DARK, col=NA)),
              right_annotation = upset_right_annotation(m, gp = gpar(fill = PAL_DARK, col=NA))) %>%
            draw() %>%
            grid.grabExpr() %>%
            as.ggplot()
    
    # ranges of enrichments in perturb seq for selected pathways
    plts[["corrs-term_vs_nes-pertseq-violin"]] = experiments %>%
        filter(dataset=="ReplogleWeissman2022_rpe1" & Description %in% x[["Description"]]) %>%
        mutate(Description = factor(Description, levels=levels(x[["Description"]]))) %>%
        distinct(Description, PERT_GENE, NES, enrichmentScore, dataset) %>%
        drop_na() %>%
        ggviolin(x="Description", y="NES", trim=TRUE, color=NA, fill="dataset", palette=PAL_DATASETS) +
        geom_text(
            aes(y=5, label=label),
            . %>% count(Description) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        labs(x="Top Correlating Hallmark", y="Perturb-seq NES") +
        coord_flip()
    
    # correlations (co-linearity) between NES across perturb seq for selected pathways
    nes_oi = experiments %>%
        filter(dataset=="ReplogleWeissman2022_rpe1" & Description%in%x[["Description"]]) %>%
        distinct(Description, PERT_GENE, NES, enrichmentScore, dataset, activity_diff) %>%
        drop_na()
    mat = nes_oi %>%
        pivot_wider(id_cols="PERT_GENE", names_from="Description", values_from="NES") %>%
        column_to_rownames("PERT_GENE") %>%
        cor(method="pearson", use="pairwise.complete.obs")
    
    plts[["corrs-nes_corrmat-heatmap"]] = mat %>% 
        Heatmap(
            name="Pearson Corr.",
            row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            heatmap_legend_param = list(legend_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY)),
            cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f", mat[i, j]), x, y, 
                              gp = gpar(fontsize=6, fontfamily=FONT_FAMILY))
        }) %>% 
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    # linear model: activity_diff ~ pathways
    pathways_oi = nes_oi %>% distinct(Description) %>% pull(Description)
    f = sprintf("%s ~ %s", "activity_diff", paste(pathways_oi, collapse="+"))
    data = nes_oi %>%
        pivot_wider(id_cols=c("PERT_GENE","activity_diff"), names_from="Description", values_from="NES")
    fit = lm(formula=as.formula(f), data=data)
    result = broom::tidy(fit)
    plts[["corrs-nes_lm_coefs-bar"]] = result %>% 
        mutate(
            term = factor(term, levels=levels(x[["Description"]])),
            dataset = "ReplogleWeissman2022_rpe1",
            p.value = format.pval(p.value, 2)
        ) %>%
        ggbarplot(x="term", y="estimate", color=NA, fill="dataset", palette=PAL_DATASETS) +
        geom_text(aes(y=0.4, label=p.value), size=FONT_SIZE, family=FONT_FAMILY) +
        guides(fill="none") +
        labs(x="Model Variable (Top Hallmarks)", y="Fitted Coefficient") +
        coord_flip()
    
    # venn diagram with overlaps between MYC v1 and v2
    plts[["corrs-overlap_myc_sets-venn"]] = ontologies[["hallmarks"]] %>%
        filter(term %in% c("HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2")) %>%
        bind_rows(
            splicing_factors %>%
                mutate(
                    term = "SFs",
                    gene = GENE
                ) %>%
                distinct(term,gene)
        ) %>%
        mutate(in_pathway=TRUE) %>%
        pivot_wider(id_cols="gene", names_from="term", values_from="in_pathway", values_fill=FALSE) %>%
        ggplot(aes(A=HALLMARK_MYC_TARGETS_V1, B=HALLMARK_MYC_TARGETS_V2, C=SFs)) +
        geom_venn(stroke_color=NA, 
                  fill_color=c("blue","lightblue","darkblue"),
                  set_name_size = FONT_SIZE+0.5, text_size = FONT_SIZE) +
        coord_fixed() +
        theme_void()
    
    # distribution of sfs_oi in perturb-seq
    sfs_oi = c("CBX3","DDX18")
    experiments %>%
        filter(dataset=="ReplogleWeissman2022_rpe1") %>%
        distinct(PERT_GENE, activity_diff) %>%
        arrange(-activity_diff) %>%
        mutate(
            ranking = row_number(),
            is_sf_oi = PERT_GENE %in% sfs_oi,
            is_oncogenic = PERT_GENE %in% (cancer_program %>% filter(driver_type=="Oncogenic") %>% pull(GENE))
        ) %>%
        filter(is_oncogenic) %>%
        ggscatter(x="ranking", y="activity_diff", color="is_sf_oi") +
        geom_point(data=.%>%filter(is_sf_oi))
    
    
    # expression of sfs_oi along carcinogenesis
    genes_oi = c(
        "ENSG00000122565", # CBX3
        "ENSG00000088205" # DDX18
    )
    ## bulk
    plts[["corrs-myc_target_sfs-genexpr_bulk-strip"]] = carcinogenesis_genexpr %>%
        filter(dataset=="Danielsson2013-fibroblasts" & ENSEMBL%in%genes_oi) %>%
        drop_na(cell_line_name) %>%
        ggstripchart(x="cell_line_name", y="genexpr", color="ENSEMBL", size=1) +
        labs(x="Carcinogenic Stage", y="log2(TPM+1)", subtitle="Danielsson2013-fibroblasts")
        
        
    ## singlecell
    plts[["corrs-myc_target_sfs-genexpr_singlecell-strip"]] = carcinogenesis_genexpr %>%
        filter(dataset=="Hodis2022-invitro_eng_melanoc" & ENSEMBL%in%genes_oi) %>%
        ggstripchart(x="treatment", y="genexpr", color="ENSEMBL", size=1) +
        labs(x="Carcinogenic Stage", y="log2(TPM+1)", subtitle="Hodis2022-invitro_eng_melanoc")
    
    return(plts)
}


plot_urbanski = function(urbanski_genexpr, urbanski_activity, urbanski_hallmarks){
    plts = list()
    
    X = urbanski_activity
    
    genes_oi = c(
        "ENSG00000122565", # CBX3
        "ENSG00000088205" # DDX18
    )
    
    # carcinogenic splicing switch
    plts[["urbanski-time_vs_activity-box"]] = X %>%
        drop_na(driver_type) %>%
        group_by(driver_type, condition, replicate, pert_time) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        ggstripchart(x="pert_time", y="activity", color="driver_type", size=1,
                     palette=PAL_DRIVER_TYPE, position=position_dodge(0.9)) +
        geom_boxplot(aes(color=driver_type), fill=NA) +
        facet_wrap(~condition, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Perturbation Time", y="Protein Activity", color="Driver Type")
    
    # activity of genes of interest
    plts[["urbanski-time_vs_activity_genes_oi-box"]] = X %>%
        filter(regulator %in% genes_oi) %>%
        ggstripchart(x="pert_time", y="activity", color="condition", size=1,
                     palette="Dark2", position=position_dodge(0.9)) +
        geom_boxplot(aes(color=condition), fill=NA) +
        facet_wrap(~regulator, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Perturbation Time", y="Protein Activity", color="Splicing Factor")
    
    # expression of genes of interest
    X = urbanski_genexpr
    plts[["urbanski-time_vs_genexpr-box"]] = X %>%
        filter(ID %in% genes_oi) %>%
        ggstripchart(x="pert_time", y="genexpr_logtpm", color="condition", size=1, position=position_jitterdodge(0.2, dodge.width=0.9), palette="Dark2") +
        geom_boxplot(aes(color=condition), fill=NA, position=position_dodge(0.9)) +
        facet_wrap(~ID, ncol=1, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Perturbation Time", y="log2(TPM + 1)")
    
    plts[["urbanski-genexpr_fc_distr_sfs-violin"]] = X %>%
        filter(ID %in% splicing_factors[["ENSEMBL"]]) %>%
        group_by(ID, condition, pert_time) %>%
        summarize(genexpr_logtpm=mean(genexpr_logtpm)) %>%
        ungroup() %>%
        pivot_wider(values_from="genexpr_logtpm", names_from="pert_time") %>%
        mutate(
            genexpr_fc = `8h` - `0h`,
            is_oi = ID %in% genes_oi
        ) %>%
        ggviolin(x="condition", y="genexpr_fc", color=NA, fill="condition") +
        geom_text_repel(
            aes(label=label),
            . %>% filter(is_oi) %>% mutate(label=ID),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Condition", y="logFC Gene Expression 8h vs 0h")
    
    # enrichment of pathways
    pathways_oi = c("HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2")
    X = urbanski_hallmarks
    plts[["urbanski-time_vs_hallmarks-box"]] = X %>%
        filter(Description %in% pathways_oi) %>%
        ggstripchart(x="pert_time", y="NES", color="condition", size=1, position=position_jitterdodge(0.2, dodge.width=0.9), palette="Dark2") +
        geom_boxplot(aes(color=condition), fill=NA, position=position_dodge(0.9)) +
        facet_wrap(~Description, ncol=1, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Perturbation Time", y="NES")
    
    return(plts)
}

make_plots = function(corrs, experiments, ontologies, cancer_program, 
                      urbanski_genexpr, urbanski_activity, urbanski_hallmarks){
    plts = list(
        plot_corrs(corrs, experiments, ontologies, cancer_program),
        plot_urbanski(urbanski_genexpr, urbanski_activity, urbanski_hallmarks)
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
    save_plt(plts, "corrs-bulk_vs_singlecell-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "corrs-bulk_vs_singlecell-bar", '.pdf', figs_dir, width=7, height=5)
    save_plt(plts, "corrs-hallmark_sets_overlaps-upset", '.pdf', figs_dir, width=14, height=10)
    save_plt(plts, "corrs-term_vs_nes-pertseq-violin", '.pdf', figs_dir, width=7, height=5)
    save_plt(plts, "corrs-nes_lm_coefs-bar", '.pdf', figs_dir, width=7, height=5.5)
    save_plt(plts, "corrs-overlap_myc_sets-venn", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "corrs-myc_target_sfs-genexpr_bulk-strip", '.pdf', figs_dir, width=5, height=4.5)
    save_plt(plts, "corrs-myc_target_sfs-genexpr_singlecell-strip", '.pdf', figs_dir, width=5, height=4.5)

    save_plt(plts, "urbanski-time_vs_activity-box", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "urbanski-time_vs_activity_genes_oi-box", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "urbanski-time_vs_genexpr-box", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "urbanski-genexpr_fc_distr_sfs-violin", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "urbanski-time_vs_hallmarks-box", '.pdf', figs_dir, width=5, height=7)
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
        make_option("--msigdb_dir", type="character"),
        make_option("--splicing_factors_file", type="character"),
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
    msigdb_dir = args[["msigdb_dir"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    set.seed(RANDOM_SEED)
    
    # load
    carcinogenesis_bulk_genexpr = read_tsv(carcinogenesis_bulk_genexpr_file)
    carcinogenesis_bulk_activity = read_tsv(carcinogenesis_bulk_activity_file)
    carcinogenesis_bulk_metadata = read_tsv(carcinogenesis_bulk_metadata_file)
    carcinogenesis_bulk_hallmarks = read_tsv(carcinogenesis_bulk_hallmarks_file)
    
    carcinogenesis_singlecell_genexpr = read_tsv(carcinogenesis_singlecell_genexpr_file)
    carcinogenesis_singlecell_activity = read_tsv(carcinogenesis_singlecell_activity_file)
    carcinogenesis_singlecell_metadata = read_tsv(carcinogenesis_singlecell_metadata_file)
    carcinogenesis_singlecell_hallmarks = read_tsv(carcinogenesis_singlecell_hallmarks_file)
    
    pertseq_activity = read_tsv(pertseq_activity_file)
    pertseq_hallmarks = read_tsv(pertseq_hallmarks_file)
    ontologies = load_ontologies(msigdb_dir)
    splicing_factors = read_tsv(splicing_factors_file)
    cancer_program = read_tsv(cancer_program_file) %>%
        mutate(driver_type = factor(driver_type, levels=names(PAL_DRIVER_TYPE)))
    
    urbanski_metadata = read_tsv(urbanski_metadata_file)
    urbanski_genexpr = read_tsv(urbanski_genexpr_file)
    urbanski_ex = read_tsv(urbanski_ex_file)
    urbanski_activity = read_tsv(urbanski_activiy_file)     
    urbanski_hallmarks = read_tsv(urbanski_hallmarks_file)
    
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
    experiments = bind_rows(carcinogenesis_bulk, carcinogenesis_singlecell, pertseq) %>%
        mutate(
            is_sf = PERT_GENE %in% splicing_factors[["GENE"]],
            gene_type = case_when(
                is_sf & is.na(driver_type) ~ "Non-driver SF",
                is_sf & driver_type=="Oncogenic" ~ "Oncogenic",
                is_sf & driver_type=="Tumor suppressor" ~ "Tumor suppressor",
                !is_sf & !is.na(PERT_GENE) ~ "Not SF"
            )
        )
    corrs = bind_rows(corrs_bulk, corrs_singlecell, corrs_pertseq)
    carcinogenesis_genexpr = carcinogenesis_bulk_genexpr %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr") %>%
        rename(ENSEMBL=ID) %>%
        left_join(carcinogenesis_bulk_metadata %>% distinct(run_accession, cell_line_name) %>% drop_na(run_accession), by=c("sampleID"="run_accession")) %>%
        mutate(
            cell_line_name = factor(cell_line_name, levels=LABS_BULK),
            dataset = "Danielsson2013-fibroblasts"
        ) %>%
        bind_rows(
            carcinogenesis_singlecell_genexpr %>%
            pivot_longer(-ENSEMBL, names_to="sampleID", values_to="genexpr") %>%
            left_join(carcinogenesis_singlecell_metadata, by=c("sampleID"="condition")) %>%
            mutate(
                treatment = factor(treatment, levels=LABS_SC),
                dataset = "Hodis2022-invitro_eng_melanoc"
            )
        )
    
    # enrichment
    hallmarks_oi = c("HALLMARK_MYC_TARGETS_V2","HALLMARK_MYC_TARGETS_V1")
    genes = experiments %>% 
        filter(dataset=="ReplogleWeissman2022_rpe1" & Description%in%hallmarks_oi) %>% 
        distinct(NES,PERT_GENE,Description) %>%
        group_by(PERT_GENE) %>%
        summarize(NES=mean(NES)) %>%
        ungroup() %>%
        arrange(-NES) %>% 
        distinct(PERT_GENE, NES) %>%
        drop_na(PERT_GENE) %>% 
        deframe() 
    enrichment = GSEA(geneList = genes, TERM2GENE=ontologies[["GO_BP"]], seed=RANDOM_SEED)
    ## perturbing cell cycle genes activates MYC pathways
    ## perturbing splicing factors does not specifically activate MYC pathways
    ## but, across perturbations the activations of the switch and MYC pathways also correlate
    
    # MYC knockdown does not alter the switch
    experiments %>% filter(PERT_GENE=="MYC") %>% distinct(PERT_GENE, activity_diff)
    
    # sort PERT_GENES by distances from perfect correlation between average MYC NES and switch
    x = experiments %>% 
        filter(dataset=="ReplogleWeissman2022_rpe1" & Description%in%hallmarks_oi) %>% 
        distinct(NES,PERT_GENE,Description, activity_diff) %>%
        group_by(PERT_GENE, activity_diff) %>%
        summarize(NES=mean(NES)) %>%
        ungroup() %>%
        mutate(
            NES = (NES - mean(NES))/sd(NES),
            activity_diff = (activity_diff - mean(activity_diff))/sd(activity_diff)
        ) %>%
        drop_na(PERT_GENE)
    fit = lm(activity_diff ~ NES, data=x)
    genes = setNames(1/abs(fit[["residuals"]]), x[["PERT_GENE"]])
    genes = sort(genes, decreasing=TRUE)
    enrichment = GSEA(geneList = genes, TERM2GENE=ontologies[["GO_BP"]], seed=RANDOM_SEED, scoreType="pos")
    
    # prep Urbanski2022
    urbanski_metadata = urbanski_metadata %>%
        separate(sample_title, into=c("condition","pert_time","replicate")) %>%
        distinct(run_accession, condition, pert_time, replicate) %>%
        mutate(pert_time = factor(pert_time, levels=c("0h","8h","24h","48h")))
    urbanski_genexpr = urbanski_genexpr %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr_logtpm") %>%
        left_join(urbanski_metadata, by=c("sampleID"="run_accession"))
    urbanski_activity = urbanski_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(urbanski_metadata, by=c("sampleID"="run_accession")) %>%
        left_join(cancer_program %>% distinct(ENSEMBL,GENE,driver_type), by=c("regulator"="ENSEMBL"))
    urbanski_hallmarks = urbanski_hallmarks %>%
        left_join(urbanski_metadata, by=c("sampleID"="run_accession"))
    
    # plot
    plts = make_plots(corrs, experiments, ontologies, cancer_program,
                      urbanski_genexpr, urbanski_activity, urbanski_hallmarks)
    
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

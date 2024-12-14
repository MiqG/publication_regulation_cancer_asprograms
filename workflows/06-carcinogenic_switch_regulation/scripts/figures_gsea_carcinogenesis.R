Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)
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
THRESH_FDR = 0.05
LABS_BULK = c("BJ_PRIMARY","BJ_IMMORTALIZED","BJ_TRANSFORMED","BJ_METASTATIC")
LABS_SC = c('WT','C','CB','CBT_228','CBT3','CBTA','CBTP','CBTP3','CBTPA')

RANDOM_SEED = 1234

LABS_HALLMARKS = c(
    'HALLMARK_INFLAMMATORY_RESPONSE',
    'HALLMARK_COMPLEMENT',
    'HALLMARK_IL6_JAK_STAT3_SIGNALING',
    'HALLMARK_P53_PATHWAY',
    'HALLMARK_UV_RESPONSE_DN',
    'HALLMARK_SPERMATOGENESIS',
    'HALLMARK_E2F_TARGETS',
    'HALLMARK_G2M_CHECKPOINT',
    'HALLMARK_MYC_TARGETS_V1',
    'HALLMARK_MYC_TARGETS_V2'
)

TISSUES = c(
    'Heart',
    'Forebrain',
    'Hindbrain',
    'Liver',
    'Kidney',
    'Ovary',
    'Testis',
    'Brain',
    # 'KidneyTestis',
    'Cerebellum'
)

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
# VIPER_SPLICING_DIR = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets")
# CARCINOGENESIS_BULK_DIR = file.path(ROOT,"results","new_empirical_network")
# CARCINOGENESIS_SC_DIR = file.path(ROOT,"results","activity_estimation_w_genexpr")
# RESULTS_DIR = file.path(ROOT,"results","carcinogenic_switch_regulation")

# carcinogenesis_bulk_genexpr_file = file.path(VIPER_SPLICING_DIR,'genexpr_tpm','tumorigenesis.tsv.gz')
# carcinogenesis_bulk_activity_file = file.path(CARCINOGENESIS_BULK_DIR,"figures","carcinogenesis","figdata","carcinogenesis","protein_activity.tsv.gz")
# carcinogenesis_bulk_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","tumorigenesis-genexpr-hallmarks.tsv.gz")
# carcinogenesis_bulk_metadata_file = file.path(VIPER_SPLICING_DIR,"metadata","tumorigenesis.tsv.gz")

# carcinogenesis_singlecell_genexpr_file = file.path(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-pseudobulk.tsv.gz")
# carcinogenesis_singlecell_activity_file = file.path(CARCINOGENESIS_SC_DIR,"figures","eval_carcinogenesis","figdata","eval_carcinogenesis","protein_activity_singlecell.tsv.gz")
# carcinogenesis_singlecell_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","Hodis2022-invitro_eng_melanoc-hallmarks.tsv.gz")
# carcinogenesis_singlecell_metadata_file = file.path(PREP_DIR,"singlecell","Hodis2022-invitro_eng_melanoc-conditions.tsv.gz")

# pertseq_activity_file = file.path(RESULTS_DIR,"figures","upstream_regulators","figdata","upstream_regulators","cancer_program_activity.tsv.gz")
# pertseq_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","ReplogleWeissman2022_rpe1-hallmarks.tsv.gz")
# pertseq_genexpr_file = file.path(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz")

# urbanski_metadata_file = file.path(PREP_DIR,"metadata","Urbanski2022.tsv.gz")
# urbanski_genexpr_file = file.path(PREP_DIR,'genexpr_tpm',"Urbanski2022.tsv.gz")
# urbanski_ex_file = file.path(PREP_DIR,'event_psi',"Urbanski2022-EX.tsv.gz")
# urbanski_activity_file = file.path(RESULTS_DIR,"files","protein_activity","Urbanski2022-EX.tsv.gz")
# urbanski_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","Urbanski2022-hallmarks.tsv.gz")

# cardoso_metadata_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets","metadata","CardosoMoreira2020.tsv.gz")
# cardoso_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","CardosoMoreira2020-hallmarks.tsv.gz")

# msigdb_dir = file.path(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# chea_file = file.path(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz")
# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")
# cancer_program_file = file.path(CARCINOGENESIS_BULK_DIR,'files','PANCAN','cancer_program.tsv.gz')
# gene_annot_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,"figures","gsea_carcinogenesis")

##### FUNCTIONS #####
load_ontologies = function(msigdb_dir, chea_file){
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "CHEA" = read.gmt(chea_file)
    )
    return(ontologies)
}

plot_corrs = function(corrs, experiments, ontologies){
    plts = list()
    
    X = corrs

    # distributions of correlations
    plts[["corrs-nes_vs_activity_diff-violin"]] = X %>%
        filter(dataset %in% c("Danielsson2013-fibroblasts","Hodis2022-invitro_eng_melanoc")) %>%
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
    
    plts[["corrs-pertseq-bar"]] = X %>%
        filter(dataset=="ReplogleWeissman2022_rpe1" & Description %in% x[["Description"]]) %>%
        mutate(Description = factor(Description, levels=LABS_HALLMARKS)) %>%
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
        mutate(Description = factor(Description, levels=LABS_HALLMARKS)) %>%
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
    
    return(plts)
}

plot_urbanski = function(urbanski_activity, urbanski_hallmarks){
    plts = list()
    
    # cancer splicing programs
    X = urbanski_activity  %>%
            drop_na(driver_type) %>%
            # median program activity in each replicate
            group_by(condition, driver_type, pert_time, replicate) %>%
            summarize(activity_norm = median(activity_norm)) %>%
            ungroup()
    plts[["urbanski-time_vs_activity_norm-box"]] = X %>%
        ggstripchart(x="pert_time", y="activity_norm", color="driver_type", size=1,
                     palette=PAL_DRIVER_TYPE, position=position_dodge(0.9)) +
        geom_boxplot(aes(color=driver_type), fill=NA) +
        facet_wrap(~condition, ncol=3) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Perturbation Time", y="Protein Activity Norm.", color="Driver Type")
    
    plts[["urbanski-time_vs_activity_norm_diff-line"]] = X %>%
        group_by(condition, driver_type, pert_time) %>%
        summarize(activity_norm = mean(activity_norm)) %>%
        ungroup() %>%
        mutate(
            activity_norm = ifelse(driver_type=="Tumor suppressor", -activity_norm, activity_norm)
        ) %>%
        group_by(condition, pert_time) %>%
        summarize(activity_diff = sum(activity_norm)) %>%
        ungroup() %>%
        filter(condition=="MYCER") %>%
        ggline(
            x="pert_time", y="activity_diff", color=PAL_DARK, numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.05
        ) +
        geom_hline(yintercept=0, linetype="dashed", linewidth=LINE_SIZE) +
        labs(x="Perturbation Time", y="Protein Activity Norm. Diff.")
    
    # enrichment of pathways
    pathways_oi = LABS_HALLMARKS
    X = urbanski_hallmarks
    plts[["urbanski-time_vs_hallmarks-box"]] = X %>%
        filter(Description %in% pathways_oi) %>%
        mutate(
            Description = factor(Description, levels=rev(pathways_oi))
        ) %>%
        ggstripchart(x="pert_time", y="NES", color="condition", size=1, position=position_jitterdodge(0.2, dodge.width=0.9), palette="Dark2") +
        geom_boxplot(aes(color=condition), fill=NA, position=position_dodge(0.9)) +
        facet_wrap(~Description, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Perturbation Time", y="NES")

    ctls = X %>%
        filter(Description %in% pathways_oi & condition=="MCF10A") %>%
        group_by(Description, pert_time) %>%
        summarize(ctl_nes = mean(NES)) %>%
        ungroup()
    plts[["urbanski-time_vs_hallmarks_norm-box"]] = X %>%
        filter(Description%in%pathways_oi & condition!="MCF10A") %>%
        left_join(ctls, by=c("Description","pert_time")) %>%
        mutate(
            nes_diff = NES - ctl_nes,
            Description = factor(Description, levels=rev(pathways_oi)) 
        ) %>%
        ggstripchart(x="pert_time", y="nes_diff", color="condition", size=1, position=position_jitterdodge(0.2, dodge.width=0.9), palette="Dark2") +
        geom_boxplot(aes(color=condition), fill=NA, position=position_dodge(0.9)) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        facet_wrap(~Description) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Perturbation Time", y="NES Norm.")
    
    return(plts)
}

plot_myc = function(enrichments, gene_corrs, activity_corrs, pertseq_activity){
    plts = list()
    
    # Are MYC targets enriched in splicing factors? YES
    plts[["myc-enrichments-go_bp-dot"]] = enrichments[["GO_BP"]] %>%
        as.data.frame() %>%
        filter(p.adjust < THRESH_FDR) %>%
        rowwise() %>%
        mutate(GeneRatio = eval(parse(text=GeneRatio))) %>%
        ungroup() %>%
        arrange(GeneRatio, Description) %>%
        slice_max(GeneRatio, n=10) %>%
        mutate(Description = factor(Description, levels=rev(unique(Description)))) %>%
        ggscatter(x="GeneRatio", y="Description", size="Count", color="p.adjust") +
        scale_size(range=c(0.5,3)) + 
        scale_color_continuous(
            low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
            name="FDR", guide=guide_colorbar(reverse=TRUE)
        )
    
    # distributions of gene correlations
    X = gene_corrs %>%
        drop_na(correlation) %>%
        group_by(ENSEMBL, `Approved symbol`, gene_type) %>%
        summarize(
            correlation = median(correlation, na.rm=TRUE),
            n_obs = n()
        ) %>%
        ungroup() %>%
        arrange(correlation) %>%
        mutate(
            ranking=row_number(),
            gene_type = factor(gene_type, levels=names(PAL_GENE_TYPE))
        )
    
    plts[["myc-gene_type_vs_genexpr_pearson-violin"]] = X %>%
        ggviolin(x="gene_type", y="correlation", color=NA, fill="gene_type", palette=PAL_GENE_TYPE, trim=TRUE) +
        geom_boxplot(width=0.25, color="black", fill=NA, outlier.size = 0.1) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        stat_compare_means(ref.group="Not SF", method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY, label="p.format") +
        geom_text(
            aes(y=-1.05, label=label),
            . %>% count(gene_type) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        labs(x="Gene Type", y="Aggregated Pearson Corr. Coef.")
    
    plts[["myc-consensus_ranking_vs_genexpr_pearson-scatter"]] = X %>%
        ggscatter(x="ranking", y="correlation", color="gene_type", palette=PAL_GENE_TYPE, size=1) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        facet_wrap(~gene_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Ranking", y="Aggregated Pearson Corr. Coef.")
    
    # distributions of SF activity correlations
    X = activity_corrs %>%
        filter(dataset=="Urbanski2022") %>%
        distinct(dataset, ENSEMBL, correlation, `Approved symbol`, gene_type) %>%
        drop_na(correlation) %>%
        arrange(correlation) %>%
        mutate(
            ranking=row_number(),
            gene_type = factor(gene_type, levels=names(PAL_GENE_TYPE)),
            label = sprintf("%s (%s)", `Approved symbol`, ranking)
        ) %>%
        left_join(
            pertseq_activity %>% distinct(PERT_ENSEMBL, PERT_GENE, activity_diff, activity_type),
            by = c("ENSEMBL"="PERT_ENSEMBL")
        )
    
    plts[["myc-gene_type_vs_activity_pearson-bar"]] = X %>%
        ggbarplot(x="ranking", y="correlation", color=NA, fill="gene_type", palette=PAL_GENE_TYPE) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        geom_text_repel(
            aes(label=label, color=gene_type),
            . %>% slice_max(correlation, n=5),
            size=FONT_SIZE, family=FONT_FAMILY, 
            segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=label, color=gene_type),
            . %>% slice_min(correlation, n=5),
            size=FONT_SIZE, family=FONT_FAMILY, 
            segment.size=0.1, max.overlaps=50
        ) +
        labs(x="Ranking", y="Pearson Corr. Coef.", fill="Gene Type", color="Gene Type")
    
    plts[["myc-urbanski_correlation_vs_pertseq_activity_diff-scatter"]] = X %>%
        ggplot(aes(x=activity_diff, y=correlation)) +
        geom_scattermore(aes(color=gene_type), pixels = c(1000,1000), pointsize=8, alpha=0.5) +
        color_palette(PAL_GENE_TYPE) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        geom_vline(xintercept=0, linetype="dashed", color="black") +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_text_repel(
            aes(label=PERT_GENE, color=gene_type),
            . %>% slice_max(correlation+(-activity_diff), n=3),
            size=FONT_SIZE, family=FONT_FAMILY, 
            segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=PERT_GENE, color=gene_type),
            . %>% slice_min(correlation+(-activity_diff), n=3),
            size=FONT_SIZE, family=FONT_FAMILY, 
            segment.size=0.1, max.overlaps=50
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Protein Activity Norm. Diff. RPE1 (Replogle2020)",
             y="Pearson Corr. Coef. (Urbanski2020)")
    
    return(plts)
}

plot_cardoso = function(cardoso_hallmarks){
    plts = list()
    
    X = cardoso_hallmarks %>%
        drop_na(time_norm, Description) %>%
        filter(Description%in%LABS_HALLMARKS) %>%
        group_by(tissue, Description) %>%
        mutate(time = as.numeric(cut(log10(time_norm+1), breaks=10))) %>%
        ungroup() %>%
        group_by(tissue, time, Description) %>%
        summarize(NES = median(NES)) %>%
        ungroup() %>%
        mutate(
            tissue = factor(tissue, levels=TISSUES),
            Description = factor(Description, levels=rev(LABS_HALLMARKS))
        )
    
    plts[["cardoso-differentiation_vs_hallmarks-line"]] = X  %>%
        ggline(
            x="time", y="NES", color="tissue", numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.01
        ) +
        color_palette(get_palette("Dark2", length(TISSUES))) +
        geom_hline(yintercept=0, color="black", linetype="dashed", linewidth=LINE_SIZE) +
        stat_cor(aes(color=tissue), method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY, label.y = -1.5) +
        facet_wrap(~Description) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="log10(Days Post Conception + 1) Binned", y="Median NES", color="Tissue")
    
    hallmarks_oi = c(
        'HALLMARK_MYC_TARGETS_V1',
        'HALLMARK_MYC_TARGETS_V2'
    )
    plts[["cardoso-differentiation_vs_hallmarks-myc_only-line"]] = X %>%
        filter(Description %in% hallmarks_oi) %>%
        ggscatter(
            x="time", y="NES", color="Description", size=0.1, alpha=0.5, position=position_jitter(0.1)
        ) +
        geom_smooth(aes(color=Description), linewidth=LINE_SIZE, linetype="dashed", fill="lightgrey") +
        color_palette("Dark2") +
        geom_hline(yintercept=0, color="black", linetype="dashed", linewidth=LINE_SIZE) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY, label.y = -1.5) +
        facet_wrap(~Description, nrow=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        scale_x_continuous(breaks = unique(X[["time"]])) +
        labs(x="log10(Days Post Conception + 1) Binned", y="Median NES", color="Tissue")
    
    return(plts)
}

make_plots = function(corrs, experiments, ontologies, urbanski_activity, urbanski_hallmarks, 
                      enrichments, gene_corrs, activity_corrs, pertseq_activity, cardoso_hallmarks){
    plts = list(
        plot_corrs(corrs, experiments, ontologies),
        plot_urbanski(urbanski_activity, urbanski_hallmarks),
        plot_myc(enrichments, gene_corrs, activity_corrs, pertseq_activity),
        plot_cardoso(cardoso_hallmarks)
    )
    plts = do.call(c,plts)
    return(plts)
}

make_figdata = function(corrs, experiments, ontologies, urbanski_activity, urbanski_hallmarks, 
                        enrichments, gene_corrs, activity_corrs, pertseq_activity, cardoso_hallmarks){
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
    save_plt(plts, "corrs-nes_vs_activity_diff-violin", '.pdf', figs_dir, width=3, height=4)
    save_plt(plts, "corrs-bulk_vs_singlecell-bar", '.pdf', figs_dir, width=7, height=5)
    save_plt(plts, "corrs-pertseq-bar", '.pdf', figs_dir, width=7, height=5)
    save_plt(plts, "corrs-hallmark_sets_overlaps-upset", '.pdf', figs_dir, width=14, height=10)
    save_plt(plts, "corrs-term_vs_nes-pertseq-violin", '.pdf', figs_dir, width=7, height=5)

    save_plt(plts, "urbanski-time_vs_activity_norm-box", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "urbanski-time_vs_activity_norm_diff-line", '.pdf', figs_dir, width=5, height=4)
    save_plt(plts, "urbanski-time_vs_hallmarks-box", '.pdf', figs_dir, width=12, height=12)
    save_plt(plts, "urbanski-time_vs_hallmarks_norm-box", '.pdf', figs_dir, width=8, height=9)
    
    save_plt(plts, "myc-enrichments-go_bp-dot", '.pdf', figs_dir, width=10, height=5)
    save_plt(plts, "myc-gene_type_vs_genexpr_pearson-violin", '.pdf', figs_dir, width=5, height=4)
    save_plt(plts, "myc-consensus_ranking_vs_genexpr_pearson-scatter", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "myc-gene_type_vs_activity_pearson-bar", '.pdf', figs_dir, width=9, height=5)
    save_plt(plts, "myc-urbanski_correlation_vs_pertseq_activity_diff-scatter", '.pdf', figs_dir, width=4, height=6)
    
    save_plt(plts, "cardoso-differentiation_vs_hallmarks-line", '.pdf', figs_dir, width=15, height=15)
    save_plt(plts, "cardoso-differentiation_vs_hallmarks-myc_only-line", '.pdf', figs_dir, width=11, height=6)
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
        make_option("--carcinogenesis_bulk_genexpr_file", type="character"),
        make_option("--carcinogenesis_bulk_activity_file", type="character"),
        make_option("--carcinogenesis_bulk_metadata_file", type="character"),
        make_option("--carcinogenesis_bulk_hallmarks_file", type="character"),
        make_option("--carcinogenesis_singlecell_genexpr_file", type="character"),
        make_option("--carcinogenesis_singlecell_activity_file", type="character"),
        make_option("--carcinogenesis_singlecell_hallmarks_file", type="character"),
        make_option("--carcinogenesis_singlecell_metadata_file", type="character"),
        make_option("--pertseq_activity_file", type="character"),
        make_option("--pertseq_hallmarks_file", type="character"),
        make_option("--pertseq_genexpr_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--chea_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--cancer_program_file", type="character"),
        make_option("--gene_annot_file", type="character"),
        make_option("--urbanski_metadata_file", type="character"),
        make_option("--urbanski_genexpr_file", type="character"),
        make_option("--urbanski_ex_file", type="character"),
        make_option("--urbanski_activity_file", type="character"),
        make_option("--urbanski_hallmarks_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    carcinogenesis_bulk_genexpr_file = args[["carcinogenesis_bulk_genexpr_file"]]
    carcinogenesis_bulk_activity_file = args[["carcinogenesis_bulk_activity_file"]]
    carcinogenesis_bulk_metadata_file = args[["carcinogenesis_bulk_metadata_file"]]
    carcinogenesis_bulk_hallmarks_file = args[["carcinogenesis_bulk_hallmarks_file"]]
    carcinogenesis_singlecell_genexpr_file = args[["carcinogenesis_singlecell_genexpr_file"]]
    carcinogenesis_singlecell_activity_file = args[["carcinogenesis_singlecell_activity_file"]]
    carcinogenesis_singlecell_hallmarks_file = args[["carcinogenesis_singlecell_hallmarks_file"]]
    carcinogenesis_singlecell_metadata_file = args[["carcinogenesis_singlecell_metadata_file"]]
    pertseq_activity_file = args[["pertseq_activity_file"]]
    pertseq_hallmarks_file = args[["pertseq_hallmarks_file"]]
    pertseq_genexpr_file = args[["pertseq_genexpr_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    chea_file = args[["chea_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    cancer_program_file = args[["cancer_program_file"]]
    gene_annot_file = args[["gene_annot_file"]]
    urbanski_metadata_file = args[["urbanski_metadata_file"]]
    urbanski_genexpr_file = args[["urbanski_genexpr_file"]]
    urbanski_ex_file = args[["urbanski_ex_file"]]
    urbanski_activity_file = args[["urbanski_activity_file"]]
    urbanski_hallmarks_file = args[["urbanski_hallmarks_file"]]
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
    pertseq_genexpr = read_tsv(pertseq_genexpr_file)
    
    ontologies = load_ontologies(msigdb_dir, chea_file)
    splicing_factors = read_tsv(splicing_factors_file)
    cancer_program = read_tsv(cancer_program_file) %>%
        mutate(driver_type = factor(driver_type, levels=names(PAL_DRIVER_TYPE)))
    gene_annot = read_tsv(gene_annot_file)
    
    urbanski_metadata = read_tsv(urbanski_metadata_file)
    urbanski_genexpr = read_tsv(urbanski_genexpr_file)
    urbanski_ex = read_tsv(urbanski_ex_file)
    urbanski_activity = read_tsv(urbanski_activity_file)     
    urbanski_hallmarks = read_tsv(urbanski_hallmarks_file)
    
    cardoso_metadata = read_tsv(cardoso_metadata_file)
    cardoso_hallmarks = read_tsv(cardoso_hallmarks_file)
    
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
                    activity = median(activity)
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
        ungroup() %>%
        left_join(
            carcinogenesis_singlecell_activity %>%
                filter(dataset_id=="Hodis2022_invitro_eng_melanoc-bulkgenexpr-adjusted_fclayer") %>%
                group_by(treatment, driver_type, cell_type) %>%
                summarize(
                    activity = median(activity)
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
        separate(sampleID, sep="___", into=c("study","cell_line","PERT_ENSEMBL","PERT_TYPE"), remove=FALSE) %>%
        left_join(
            pertseq_activity %>% filter(activity_type=="activity_rpe1"), 
            by="PERT_ENSEMBL"
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
        ungroup() %>%
        arrange(correlation_diff_activity)
    
    # merge
    ## activities
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
    ## correlations
    corrs = bind_rows(corrs_bulk, corrs_singlecell, corrs_pertseq)
    
    ## gene expression
    carcinogenesis_genexpr = carcinogenesis_bulk_genexpr %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr") %>%
        rename(ENSEMBL=ID) %>%
        left_join(
            carcinogenesis_bulk_metadata %>% 
                distinct(run_accession, cell_line_name) %>% 
                drop_na(run_accession), 
            by=c("sampleID"="run_accession")
        ) %>%
        drop_na(cell_line_name) %>%
        mutate(
            cell_line_name = factor(cell_line_name, levels=LABS_BULK),
            dataset = "Danielsson2013-fibroblasts"
        ) %>%
        left_join(
            carcinogenesis_bulk %>% distinct(dataset, cell_line_name, activity_diff),
            by=c("dataset","cell_line_name")
        ) %>%
        bind_rows(
            carcinogenesis_singlecell_genexpr %>%
            pivot_longer(-ENSEMBL, names_to="sampleID", values_to="genexpr") %>%
            left_join(carcinogenesis_singlecell_metadata, by=c("sampleID"="condition")) %>%
            mutate(
                treatment = factor(treatment, levels=LABS_SC),
                dataset = "Hodis2022-invitro_eng_melanoc"
            ) %>%
            left_join(
                    carcinogenesis_singlecell %>% distinct(dataset, treatment, activity_diff),
                by=c("dataset","treatment")
            )
        )
    
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

    # normalize Urbanski by condition and time
    ## activity
    urbanski_activity = urbanski_activity %>%
        group_by(regulator, GENE, driver_type, condition, replicate, pert_time) %>%
        summarize(activity = median(activity)) %>%
        ungroup()
    ctls_cond = urbanski_activity %>%
        filter(condition=="MCF10A") %>%
        group_by(regulator, GENE, driver_type, pert_time) %>%
        summarize(ctl_activity = mean(activity)) %>%
        ungroup()
    urbanski_activity = urbanski_activity %>%
        left_join(ctls_cond, by=c("regulator","GENE","driver_type","pert_time")) %>%
        mutate(activity_norm = activity - ctl_activity)   
    
    ## gene expression
    ctls_cond = urbanski_genexpr %>%
        filter(condition=="MCF10A") %>%
        group_by(pert_time, ID) %>%
        summarize(ctl_genexpr = mean(genexpr_logtpm, na.rm=TRUE)) %>%
        ungroup()
    urbanski_genexpr = urbanski_genexpr %>%
        left_join(ctls_cond, by=c("pert_time","ID")) %>%
        mutate(genexpr_norm = genexpr_logtpm - ctl_genexpr)
    ctls_time = urbanski_genexpr %>%
        filter(pert_time=="0h") %>%
        group_by(condition, ID) %>%
        summarize(ctl_genexpr_time = mean(genexpr_norm, na.rm=TRUE)) %>%
        ungroup()
    urbanski_genexpr = urbanski_genexpr %>%
        left_join(ctls_time, by=c("condition","ID")) %>%
        mutate(genexpr_norm_time = genexpr_norm - ctl_genexpr_time)
    
    # gene correlations: correlate gene expression changes of MYC targets with changes in switch activation
    genes_oi = ontologies[["CHEA"]] %>%
        filter(term=="MYC") %>%
        pull(gene)
    ensembl_oi = gene_annot %>%
        filter(`Approved symbol` %in% genes_oi) %>%
        pull(`Ensembl gene ID`)
    
    ## carcinogenesis datasets
    gene_corrs_carcinogenesis = carcinogenesis_genexpr %>%
        # only MYC targets
        filter(ENSEMBL %in% ensembl_oi) %>%
        # average gene expression across condition replicates
        group_by(dataset, ENSEMBL, activity_diff, cell_line_name, treatment) %>%
        summarize(genexpr = mean(genexpr, na.rm=TRUE)) %>%
        ungroup() %>%
        drop_na(genexpr, activity_diff) %>%
        # correlate gene expression with program activity differences
        group_by(dataset, ENSEMBL) %>%
        summarize(correlation = cor(genexpr, activity_diff, method="pearson")) %>%
        ungroup()
    
    ## perturb seq dataset
    pertseq_genexpr_mat = pertseq_genexpr %>%
        # only MYC targets
        filter(ensembl_id %in% ensembl_oi) %>%
        column_to_rownames("ensembl_id") %>%
        as.data.frame()
    colnames(pertseq_genexpr_mat) = gsub("ReplogleWeissman2022_rpe1___RPE1___","",colnames(pertseq_genexpr_mat))
    colnames(pertseq_genexpr_mat) = gsub("___KNOCKDOWN","",colnames(pertseq_genexpr_mat))
    pertseq_activity_vec = pertseq_activity %>% 
        filter(activity_type=="activity_rpe1") %>%
        distinct(PERT_ENSEMBL, activity_diff) %>%
        deframe()
    gene_corrs_pertseq = cor(
            t(pertseq_genexpr_mat), pertseq_activity_vec[colnames(pertseq_genexpr_mat)], 
            method="pearson"
        )[,1] %>%
        enframe("ENSEMBL", "correlation") %>%
        mutate(dataset = "ReplogleWeissman2022_rpe1")
    
    ## MYC activation dataset (Urbanski2022)
    gene_corrs_urbanski = urbanski_genexpr %>%
        rename(ENSEMBL=ID) %>%
        filter(ENSEMBL%in%ensembl_oi & condition=="MYCER") %>%
        # average normalized gene expression across replicates
        group_by(condition, ENSEMBL, pert_time) %>%
        summarize(genexpr_norm_time = mean(genexpr_norm_time, na.rm=TRUE)) %>%
        ungroup() %>%
        # add program activity difference
        left_join(
            urbanski_activity %>%
                drop_na(driver_type) %>%
                # median program activity in each replicate
                group_by(condition, driver_type, pert_time, replicate) %>%
                summarize(activity_norm = median(activity_norm)) %>%
                ungroup() %>%
                # average normalized activity across replicates
                group_by(condition, driver_type, pert_time) %>%
                summarize(activity_norm = mean(activity_norm, na.rm=TRUE)) %>%
                ungroup() %>%
                # compute activity difference
                mutate(activity_norm = ifelse(driver_type=="Tumor suppressor", -activity_norm, activity_norm)) %>%
                group_by(condition, pert_time) %>%
                summarize(activity_diff = sum(activity_norm)) %>%
                ungroup(),
            by = c("condition","pert_time")
        ) %>%
        # compute correlation
        drop_na(genexpr_norm_time, activity_diff) %>%
        group_by(ENSEMBL) %>%
        summarize(correlation = cor(genexpr_norm_time, activity_diff, method="pearson")) %>%
        ungroup() %>%
        drop_na(correlation) %>%
        mutate(dataset = "Urbanski2022")
    
    ## combine gene correlations
    gene_corrs = bind_rows(
        gene_corrs_carcinogenesis,
        gene_corrs_pertseq,
        gene_corrs_urbanski
    )
    
    ## annotate
    gene_corrs = gene_corrs %>%
        left_join(cancer_program %>% distinct(ENSEMBL, driver_type), by="ENSEMBL") %>%
        left_join(gene_annot %>% distinct(`Ensembl gene ID`,`Approved symbol`), by=c("ENSEMBL"="Ensembl gene ID")) %>%
        mutate(
            is_sf = ENSEMBL %in% splicing_factors[["ENSEMBL"]],
            gene_type = case_when(
                is_sf & is.na(driver_type) ~ "Non-driver SF",
                is_sf & driver_type=="Oncogenic" ~ "Oncogenic",
                is_sf & driver_type=="Tumor suppressor" ~ "Tumor suppressor",
                !is_sf ~ "Not SF",
                !is_sf & !is.na(ENSEMBL) ~ "Not SF"
            )
        )
    
    # MYC targets enrichments
    myc_targets = ontologies[["CHEA"]] %>% 
        filter(term=="MYC") %>%
        mutate(is_sf = gene %in% splicing_factors[["GENE"]])
    
    print("MYC CHEA targets that are splicing factors:")
    myc_targets %>% count(is_sf) %>% print()
    
    genes = myc_targets %>% pull(gene)
    enrichments = list(
        "GO_BP" = enricher(gene=genes, TERM2GENE=ontologies[["GO_BP"]])
    )
    
    # SF activity correlations
    ## bulk
    activity_corrs_bulk = carcinogenesis_bulk_activity %>%
        distinct(cell_line_name, regulator, activity)
    
    diffs = activity_corrs_bulk %>%
        # program activity
        left_join(cancer_program %>% distinct(ENSEMBL, driver_type), by=c("regulator"="ENSEMBL")) %>%
        drop_na(driver_type) %>%
        group_by(cell_line_name, driver_type) %>%
        summarize(activity_program = median(activity)) %>%
        ungroup() %>%
        group_by(cell_line_name) %>%
        summarize(
            activity_diff = ifelse(driver_type=="Tumor suppressor", -activity_program, activity_program),
            activity_diff = sum(activity_diff)
        ) %>%
        ungroup() %>%
        distinct()
    
    activity_corrs_bulk = activity_corrs_bulk %>%
        left_join(diffs, by="cell_line_name")  %>%
        # only MYC targets
        filter(regulator %in% ensembl_oi) %>%
        mutate(dataset = "Danielsson2013-fibroblasts") %>%
        group_by(cell_line_name, regulator) %>%
        mutate(correlation = cor(activity, activity_diff, method="pearson")) %>%
        ungroup()
    
    ## singlecell
    activity_corrs_singlecell = carcinogenesis_singlecell_activity %>%
        filter(dataset_id=="Hodis2022_invitro_eng_melanoc-bulkgenexpr-adjusted_fclayer") %>%
        distinct(treatment, regulator, activity)
    
    diffs = activity_corrs_singlecell %>%
        # program activity
        left_join(cancer_program %>% distinct(ENSEMBL, driver_type), by=c("regulator"="ENSEMBL")) %>%
        drop_na(driver_type) %>%
        group_by(treatment, driver_type) %>%
        summarize(activity_program = median(activity)) %>%
        ungroup() %>%
        group_by(treatment) %>%
        summarize(
            activity_diff = ifelse(driver_type=="Tumor suppressor", -activity_program, activity_program),
            activity_diff = sum(activity_diff)
        ) %>%
        ungroup() %>%
        distinct()
    
    activity_corrs_singlecell = activity_corrs_singlecell %>%
        left_join(diffs, by="treatment")  %>%
        # only MYC targets
        filter(regulator %in% ensembl_oi) %>%
        mutate(dataset = "Hodis2022-invitro_eng_melanoc") %>%
        group_by(treatment, regulator) %>%
        mutate(correlation = cor(activity, activity_diff, method="pearson")) %>%
        ungroup()
    
    ## MYC (Urbanski2022)
    activity_corrs_urbanski = urbanski_activity %>%
        filter(condition=="MYCER") %>%
        distinct(condition, pert_time, replicate, regulator, activity_norm)
    
    diffs = activity_corrs_urbanski %>%
        # program activity
        left_join(cancer_program %>% distinct(ENSEMBL, driver_type), by=c("regulator"="ENSEMBL")) %>%
        drop_na(driver_type) %>%
        # median program activity in each replicate
        group_by(condition, driver_type, pert_time, replicate) %>%
        summarize(activity_norm = median(activity_norm)) %>%
        ungroup() %>%
        # average normalized activity across replicates
        group_by(condition, driver_type, pert_time) %>%
        summarize(activity_norm = mean(activity_norm, na.rm=TRUE)) %>%
        ungroup() %>%
        # compute activity difference
        mutate() %>%
        group_by(condition, pert_time) %>%
        summarize(
            activity_diff = ifelse(driver_type=="Tumor suppressor", -activity_norm, activity_norm),
            activity_diff = sum(activity_diff)) %>%
        ungroup() %>%
        distinct()
    
    activity_corrs_urbanski = activity_corrs_urbanski %>%
        # average activity by splicing factor
        group_by(condition, pert_time, regulator) %>%
        summarize(activity = mean(activity_norm, na.rm=TRUE)) %>%
        ungroup() %>%
        left_join(diffs, by=c("pert_time","condition")) %>%
        # only MYC targets
        filter(regulator %in% ensembl_oi) %>%
        mutate(dataset = "Urbanski2022") %>%
        group_by(condition, regulator) %>%
        mutate(correlation = cor(activity, activity_diff, method="pearson")) %>%
        ungroup()
    
    ## combine and annotate
    activity_corrs = bind_rows(
            activity_corrs_bulk,
            activity_corrs_singlecell,
            activity_corrs_urbanski
        ) %>%
        rename(ENSEMBL=regulator) %>%
        left_join(cancer_program %>% distinct(ENSEMBL, driver_type), by="ENSEMBL") %>%
        left_join(gene_annot %>% distinct(`Ensembl gene ID`,`Approved symbol`), by=c("ENSEMBL"="Ensembl gene ID")) %>%
        mutate(
            is_sf = ENSEMBL %in% splicing_factors[["ENSEMBL"]],
            gene_type = case_when(
                is_sf & is.na(driver_type) ~ "Non-driver SF",
                is_sf & driver_type=="Oncogenic" ~ "Oncogenic",
                is_sf & driver_type=="Tumor suppressor" ~ "Tumor suppressor",
                !is_sf ~ "Not SF",
                !is_sf & !is.na(ENSEMBL) ~ "Not SF"
            )
        )
    
    # add time to CardosoMoreira2020
    cardoso_metadata = cardoso_metadata %>% 
        distinct(study_accession, sampleID, sample_title) %>%
        separate(
            sample_title, 
            c("code","organism","tissue","time","sex"), 
            sep="\\."
        ) %>%
        mutate(
            # fix carnegie stages
            time = case_when(
                time=="CS13" ~ "32dpc",
                time=="CS14" ~ "33dpc",
                time=="CS16" ~ "39dpc",
                time=="CS17" ~ "41dpc",
                time=="CS18" ~ "44dpc",
                time=="CS19" ~ "46dpc",
                time=="CS20" ~ "49dpc",
                time=="CS21" ~ "51dpc",
                time=="CS22" ~ "53dpc",
                time=="CS23" ~ "56dpc",
                TRUE ~ time
            ),
            time_units = case_when(
                str_detect(time, "dpc") ~ "dpc", # days post conception
                str_detect(time, "w") ~ "w", # weeks post conception
                str_detect(time, "dpb") ~ "dpb", # days post birth
                str_detect(time, "wpb") ~ "wpb", # weeks post birth
                str_detect(time, "mpb") ~ "mpb", # months post birth
                str_detect(time, "ypb") ~ "ypb", # years post birth
            )
        ) %>%
        # normalized time to days post conception
        rowwise() %>% 
        mutate(
            time_norm = as.numeric(gsub(time_units,"",time)),
            time_norm = case_when(
                time_units=="dpc" ~ time_norm,
                time_units=="w" ~ time_norm*7,
                time_units=="dpb" ~ (9*30)+time_norm,
                time_units=="wpb" ~ (9*30)+7*time_norm,
                time_units=="mpb" ~ (9*30)+30*time_norm,
                time_units=="ypb" ~ (9*30)+365*time_norm
            )
        ) %>%
        ungroup() %>%
        filter(tissue %in% TISSUES)
    
    cardoso_hallmarks = cardoso_hallmarks %>%
        left_join(cardoso_metadata, by="sampleID")
    
    # plot
    plts = make_plots(corrs, experiments, ontologies, urbanski_activity, urbanski_hallmarks, 
                      enrichments, gene_corrs, activity_corrs, pertseq_activity, cardoso_hallmarks)
    
    # make figdata
    figdata = make_figdata(corrs, experiments, ontologies, urbanski_activity, urbanski_hallmarks, 
                           enrichments, gene_corrs, activity_corrs, pertseq_activity, cardoso_hallmarks)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}

##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)
require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)
require(clusterProfiler)
require(ggraph)
require(tidygraph)

# variables
COSMIC_DRIVER_TYPES = c(
    "Not in COSMIC",
    "COSMIC Suppressor",
    "COSMIC Oncogenic",
    "COSMIC Dual"
)

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
    "Oncogenic"="#F6AE2D",
    "Tumor suppressor"="#6C98B3"
)
PAL_GENE_TYPE = c(
    "Not SF"="darkgreen",
    "Non-driver SF"="darkred",
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","carcinogenic_switch_regulation")

# protein_activity_rpe1_file = file.path(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_rpe1-bulkgenexpr-adjusted_fclayer.tsv.gz")
# metadata_rpe1_file = file.path(PREP_DIR,"singlecell","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-conditions.tsv.gz")
# gene_info_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# cancer_program_file = file.path(ROOT,"results","new_empirical_network",'files','PANCAN','cancer_program.tsv.gz')
# msigdb_dir = file.path(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# cosmic_genes_file = file.path(RAW_DIR,"COSMIC","cancer_gene_census.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","upstream_regulators")

# regulons_dir = file.path(ROOT,"results","new_empirical_network","files","experimentally_derived_regulons_pruned_w_viper_networks-EX")
# event_info_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# splicing_factors_file = file.path(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv")

# shortest_paths_pert_sfs_file = file.path(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_pert_splicing_factors.tsv.gz')
# shortest_paths_random_file = file.path(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_pert_splicing_factors_random.tsv.gz')
# ppi_network_file = file.path(PREP_DIR,'ppi','STRINGDB.tsv.gz')

# gsea_hallmarks_file = file.path(RESULTS_DIR,"files","gsea","ReplogleWeissman2022_rpe1-hallmarks.tsv.gz")
# gsea_reactome_file = file.path(RESULTS_DIR,"files","gsea","ReplogleWeissman2022_rpe1-reactome.tsv.gz")

##### FUNCTIONS #####
load_ontologies = function(msigdb_dir, cosmic_genes_file){
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "GO_CC" = read.gmt(file.path(msigdb_dir,"c5.go.cc.v7.4.symbols.gmt")),
        "cosmic" = read_tsv(cosmic_genes_file) %>%
            dplyr::select("Gene Symbol","Role in Cancer") %>%
            rename(gene = `Gene Symbol`, cosmic_driver_type = `Role in Cancer`) %>%
            mutate(
                cosmic_driver_type = case_when(
                    str_detect(cosmic_driver_type, "oncogene") & !str_detect(cosmic_driver_type, "TSG") ~ "COSMIC Oncogenic",
                    !str_detect(cosmic_driver_type, "oncogene") & str_detect(cosmic_driver_type, "TSG") ~ "COSMIC Suppressor",
                    str_detect(cosmic_driver_type, "oncogene") & str_detect(cosmic_driver_type, "TSG") ~ "COSMIC Dual"
                ),
                term = "COSMIC_CENSUS"
            ) %>%
        dplyr::select(term,gene,cosmic_driver_type)
    )
    return(ontologies)
}

plot_program_activity = function(cancer_program_activity){
    plts = list()
    
    X = cancer_program_activity
    x = X %>%
        pivot_wider(
            id_cols=c("PERT_ENSEMBL","PERT_GENE","cosmic_driver_type","target_in_cosmic","driver_type","gene_type"), 
            names_from="activity_type", values_from="activity_diff"
        )
    
    plts[["program_activity-diff_vs_pert_efficiency-scatter"]] = X %>%
        filter(activity_type=="activity_rpe1") %>%
        drop_na(pert_efficiency_fc, activity_diff) %>%
        ggplot(aes(x=pert_efficiency_fc, y=activity_diff)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="KD Efficiency", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
        
    
    plts[["program_activity-ts_vs_onco_by_cell_line-scatter"]] = X %>%
        ggplot(aes(x=`Tumor suppressor`, y=`Oncogenic`)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        geom_text(
            aes(label=PERT_GENE),
            . %>% 
                group_by(activity_type) %>% 
                slice_max(abs(`Tumor suppressor`*`Oncogenic`), n=1) %>%
                ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(label=PERT_GENE),
            . %>% 
                group_by(activity_type) %>% 
                slice_max(`Tumor suppressor`+`Oncogenic`, n=1) %>%
                ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() + 
        geom_vline(xintercept=0, size=LINE_SIZE , linetype="dashed", color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~activity_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Tumor Suppressor Splicing Program Activity", y="Oncogenic Splicing Program Activity")

    plts[["program_activity-diff_ranking_overview-scatter"]] = X %>%
        drop_na(activity_diff) %>%
        group_by(activity_type) %>%
        arrange(activity_diff) %>%
        mutate(
            ranking = row_number(),
            rel_ranking = ranking / n()
        ) %>%
        ungroup() %>%
        ggplot(aes(x=rel_ranking, y=activity_diff)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        theme_pubr() + 
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type) %>% 
            slice_max(activity_diff, n=3) %>% 
            ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type) %>% 
            slice_min(activity_diff, n=3) %>% 
            ungroup(),
            size=FONT_SIZE,family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        facet_wrap(~activity_type, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Relative Ranking", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
    
    
    plts[["program_activity-diff_ranking_by_cell_line_and_cosmic-scatter"]] = X %>%
        drop_na(activity_diff) %>%
        group_by(activity_type) %>%
        arrange(-activity_diff) %>%
        mutate(
            ranking = row_number(),
            rel_ranking = ranking / n()
        ) %>%
        ungroup() %>%
        ggplot(aes(x=rel_ranking, y=activity_diff)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        theme_pubr() + 
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type, cosmic_driver_type) %>% 
            slice_max(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type, cosmic_driver_type) %>% 
            slice_min(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE,family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        facet_wrap(~cosmic_driver_type+activity_type, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Relative Ranking", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
    
    plts[["program_activity-diff_ranking_by_cell_line_and_program-scatter"]] = X %>%
        drop_na(activity_diff) %>%
        group_by(activity_type) %>%
        arrange(-activity_diff) %>%
        mutate(
            ranking = row_number(),
            rel_ranking = ranking / n()
        ) %>%
        ungroup() %>%
        ggplot(aes(x=rel_ranking, y=activity_diff)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        theme_pubr() + 
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type, driver_type) %>% 
            slice_max(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type, driver_type) %>% 
            slice_min(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE,family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        facet_wrap(~driver_type+activity_type, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Relative Ranking", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
    
    plts[["program_activity-diff_ranking_by_cell_line_and_sfornot-scatter"]] = X %>%
        drop_na(activity_diff) %>%
        group_by(activity_type) %>%
        arrange(-activity_diff) %>%
        mutate(
            ranking = row_number(),
            rel_ranking = ranking / n()
        ) %>%
        ungroup() %>%
        ggplot(aes(x=rel_ranking, y=activity_diff)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        theme_pubr() + 
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type, pert_is_sf) %>% 
            slice_max(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type, pert_is_sf) %>% 
            slice_min(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE,family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        facet_wrap(~pert_is_sf+activity_type, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Relative Ranking", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
    
    comparisons = list(c("Not in COSMIC", "COSMIC Suppressor"), c("Not in COSMIC", "COSMIC Oncogenic"), c("Not in COSMIC", "COSMIC Dual"))
    plts[["program_activity-diff_vs_cosmic-violin"]] = X %>%
        mutate(
            cosmic_driver_type = replace_na(cosmic_driver_type, "Not in COSMIC"),
            cosmic_driver_type = factor(cosmic_driver_type, levels=COSMIC_DRIVER_TYPES)
        ) %>%
        ggviolin(x="cosmic_driver_type", y="activity_diff", fill="cosmic_driver_type", palette="Paired", color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        geom_text(
            aes(y = -4.5, label=label), 
            . %>% 
            count(activity_type, cosmic_driver_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(comparisons=comparisons, method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr(x.text.angle = 70) +
        facet_wrap(~activity_type, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="COSMIC Driver Type", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")

    
    comparisons = list(
        c("Not SF","Non-driver SF"),c("Not SF","Tumor suppressor"),c("Not SF","Oncogenic"),
        c("Non-driver SF","Tumor suppressor"),c("Non-driver SF","Oncogenic"),
        c("Tumor suppressor", "Oncogenic")
    )
    plts[["program_activity-diff_vs_splicing_program-violin"]] = X %>%
        mutate(gene_type = factor(gene_type, levels=names(PAL_GENE_TYPE))) %>%
        filter(!is.na(gene_type)) %>%
        ggviolin(x="gene_type", y="activity_diff", fill="gene_type", palette=PAL_GENE_TYPE, color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black") +
        geom_text(
            aes(y = -2.25, label=label), 
            . %>% 
            count(activity_type, gene_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        #stat_compare_means(comparisons=comparisons, method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr(x.text.angle = 70) +
        facet_wrap(~activity_type, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Cancer Splicing Program", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
    
    return(plts)
}

plot_gsea = function(df){
        
    plt = df %>%
        mutate(
            Count = str_count(core_enrichment, "/")+1, 
            GeneRatio=Count/setSize
        ) %>%
        arrange(NES, Description) %>%
        mutate(Description = factor(Description, levels=unique(Description))) %>%
        ggscatter(x="Description", y="NES", size="GeneRatio", color="p.adjust") +
        scale_size(range=c(0.5,3)) + 
        scale_color_continuous(
            low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
            name="FDR", guide=guide_colorbar(reverse=TRUE)) +
        theme_pubr(x.text.angle = 70)
    
    return(plt)
}


plot_enrichments = function(enrichments){
    plts = list()
    
    # GO BP
    X = enrichments[["rpe1_activity_diff_gobp"]] %>%
        as.data.frame() %>%
        slice_max(NES, n=10)
    
        # bind_rows(
        #     enrichments[["rpe1_activity_diff_gobp"]] %>%
        #         as.data.frame() %>%
        #         slice_min(NES, n=5)
        # )
        
    plts[["enrichments-activity_diff-dot"]] = X %>%
        mutate(
            Count = str_count(core_enrichment, "/")+1, 
            GeneRatio=Count/setSize
        ) %>%
        arrange(NES, Description) %>%
        mutate(Description = factor(Description, levels=unique(Description))) %>%
        ggscatter(x="NES", y="Description", size="GeneRatio", color="p.adjust") +
        scale_size(range=c(0.5,3)) + 
        scale_color_continuous(
            low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
            name="FDR", guide=guide_colorbar(reverse=TRUE))
    
    return(plts)
}

plot_sf_network_analysis = function(networks_sf_ex, splicing_factors, cancer_program){
    plts = list()
    
    X = networks_sf_ex
    
    x = X %>% 
        distinct(regulator,PERT_GENE,activity_diff,target_type,regulator_type,GENE) %>% 
        count(regulator,PERT_GENE,activity_diff,target_type, regulator_type) %>%
        group_by(regulator) %>%
        mutate(
            n_total = sum(n),
            perc = 100 * n / n_total
        ) %>%
        ungroup() %>%
        filter(n_total>=25)
    
    # overview target gene percentages for each splicing factor (barplot)
    sf_order = x %>%
        filter(target_type=="Not SF") %>%
        arrange(perc) %>%
        pull(regulator)
    plts[["sf_network_analysis-sf_vs_target_type_freq-bar"]] = x %>%
        mutate(
            target_type = factor(target_type, levels=names(PAL_GENE_TYPE)),
            regulator = factor(regulator, levels=sf_order)
        ) %>%
        ggbarplot(x="regulator", y="perc", fill="target_type", color=NA) +
        fill_palette(PAL_GENE_TYPE) +
        labs(x="Splicing Factor", y="% Target Genes")
    
    # correlation KD switch and oncogenic/suppressor %?
    plts[["sf_network_analysis-target_type_freq_vs_activity_diff-non_driver-scatter"]] = x %>%
        filter(regulator_type=="Non-driver SF") %>%
        mutate(target_type = factor(target_type, levels=names(PAL_GENE_TYPE))) %>%
        drop_na(perc, activity_diff) %>%
        ggscatter(x="perc", y="activity_diff", size=1, alpha=0.5, color="regulator_type") +
        color_palette(PAL_GENE_TYPE) +
        fill_palette(PAL_GENE_TYPE) +
        stat_cor(method="spearman987", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~target_type, scales="free") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="% Target Genes", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
    
    # how many of the possible targets are regulated by another SF?
    total_sfs = splicing_factors %>%
        left_join(cancer_program, by=c("ENSEMBL","GENE")) %>%
        mutate(
            sf_type = case_when(
                is.na(driver_type) ~ "Non-driver SF",
                driver_type=="Oncogenic" ~ "Oncogenic",
                driver_type=="Tumor suppressor" ~ "Tumor suppressor"
            )
        ) %>%
        count(sf_type, name="n_total")
    
    edges = X %>% 
        distinct(GENE, regulator_type, target_type) %>%
        filter(target_type!="Not SF") %>%
        count(regulator_type, target_type) %>%
        left_join(total_sfs, by=c("target_type"="sf_type")) %>%
        mutate(
            perc = n / n_total,
            edge_label = sprintf("%s (n=%s)", round(perc, 2), n),
            regulator_class = regulator_type,
            node_label = sprintf("%s (n=%s)", target_type, n_total)
        )
    
    nodes = edges %>%
        distinct(node_label) %>%
        mutate(regulator_class = gsub(" \\(.*", "", node_label)) %>%
        distinct(regulator_class, node_label)
    
    g = tbl_graph(nodes=nodes, edges=edges, directed=TRUE)
        
    plts[["sf_network_analysis-regulator_type_perc-graph"]] = g %>%
        ggraph() +
        geom_edge_arc(
            arrow=arrow(length=unit(3, 'mm')), 
            aes(label=edge_label, width=perc, color=regulator_class),
            start_cap=circle(3, 'mm'),
            end_cap=circle(3, 'mm'),
            strength=0.2,
            alpha=0.5,
            label_size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        geom_edge_loop(
            aes(label=edge_label, width=perc, color=regulator_class), 
            arrow=arrow(length=unit(3, 'mm')),
            alpha=0.5,
            label_size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_node_text(
            aes(label=node_label, color=regulator_class),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        color_palette(PAL_GENE_TYPE) +
        scale_edge_color_manual(values=PAL_GENE_TYPE) +
        scale_edge_width_continuous(range = c(0.2, 1)) +
        theme_void() +
        theme(aspect.ratio=1)
    
    return(plts)
}

plot_ppi_network_analysis = function(shortest_paths_pert_sfs, cancer_program_activity, splicing_factors, cancer_program, ppi_network){
    plts = list()
    
    X = shortest_paths_pert_sfs %>%
        filter(source!=target) %>%
        drop_na() %>%
        mutate(pair_type = factor(pair_type, levels=c("Both Weak","Mixed","Both Strong"))) %>%
        left_join(
            cancer_program_activity %>% 
                distinct(PERT_GENE, activity_diff) %>%
                rename(activity_diff_source=activity_diff),
            by=c("source"="PERT_GENE")
        ) %>%
        left_join(
            cancer_program_activity %>% 
                distinct(PERT_GENE, activity_diff) %>%
                rename(activity_diff_target=activity_diff),
            by=c("target"="PERT_GENE")
        ) %>%
        rowwise() %>%
        mutate(activity_diff_avg = mean(activity_diff_source,activity_diff_target)) %>%
        ungroup()
    
    comparisons = list(c("Both Weak", "Both Strong"),c("Mixed", "Both Strong"),c("Mixed","Both Weak"))
    plts[["ppi_network_analysis-pair_type_vs_path_length-violin"]] = X %>%
        ggviolin(x="pair_type", y="shortest_path_length", trim=TRUE, fill="darkkhaki", color=NA, 
                 add="median_iqr", add.params=list(color="black", size=0.1)) +
        stat_compare_means(comparisons=comparisons, method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_text(
            aes(y=0.1, label=label),
            . %>% count(pair_type) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        coord_flip() +
        labs(x="SF Pair", y="STRINGDB Shortest Path Length")
    
    corr = cor.test(as.numeric(X[["shortest_path_length"]]), X[["activity_diff_avg"]], method="spearman")
    corr_label = sprintf("R = %s, p %s", round(corr[["estimate"]],2), format.pval(corr[["p.value"]]))
    plts[["ppi_network_analysis-activity_diff_avg_vs_path_length-scatter"]] = X %>%
        ggviolin(x="shortest_path_length", y="activity_diff_avg", color=NA, fill="darkkhaki", trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black") +
        geom_text(
            aes(y = -2.25, label=label), 
            . %>% 
            count(shortest_path_length) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        annotate(
            "text", x=3, y=3.7, label=corr_label, size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        geom_hline(yintercept=0, size=LINE_SIZE , linetype="dashed", color="black") +
        labs(x="STRINGDB Shortest Path Length", y="Median Program Act. Diff.\n(Avg. by SF Pair)")
    
    # how many of the possible targets are regulated by another SF?
    sfs = splicing_factors %>%
        left_join(cancer_program, by=c("ENSEMBL","GENE")) %>%
        mutate(
            sf_type = case_when(
                is.na(driver_type) ~ "Non-driver SF",
                driver_type=="Oncogenic" ~ "Oncogenic",
                driver_type=="Tumor suppressor" ~ "Tumor suppressor"
            )
        ) %>%
        # consider only SFs available in PPI
        filter()
    
    total_sfs = sfs %>%
        count(sf_type, name="n_total")
    
    edges = ppi_network %>% 
        distinct(source, target) %>%
        left_join(
            sfs %>% distinct(GENE, sf_type) %>% rename(regulator_type=sf_type), by=c("source"="GENE")
        ) %>%
        left_join(
            sfs %>% distinct(GENE, sf_type) %>% rename(target_type=sf_type), by=c("target"="GENE")
        ) %>%
        drop_na() %>%
        distinct(target, regulator_type, target_type) %>%
        filter(target_type!="Not SF") %>%
        count(regulator_type, target_type) %>%
        left_join(total_sfs, by=c("target_type"="sf_type")) %>%
        mutate(
            perc = n / n_total,
            edge_label = sprintf("%s (n=%s)", round(perc, 2), n),
            regulator_class = regulator_type,
            node_label = sprintf("%s (n=%s)", target_type, n_total)
        )
    
    nodes = edges %>%
        distinct(node_label) %>%
        mutate(regulator_class = gsub(" \\(.*", "", node_label)) %>%
        distinct(regulator_class, node_label)
    
    g = tbl_graph(nodes=nodes, edges=edges, directed=TRUE)
        
    plts[["ppi_network_analysis-regulator_type_perc-graph"]] = g %>%
        ggraph() +
        geom_edge_arc(
            arrow=arrow(length=unit(3, 'mm')), 
            aes(label=edge_label, width=perc, color=regulator_class),
            start_cap=circle(3, 'mm'),
            end_cap=circle(3, 'mm'),
            strength=0.2,
            alpha=0.5,
            label_size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        geom_edge_loop(
            aes(label=edge_label, width=perc, color=regulator_class), 
            arrow=arrow(length=unit(3, 'mm')),
            alpha=0.5,
            label_size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_node_text(
            aes(label=node_label, color=regulator_class),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        color_palette(PAL_GENE_TYPE) +
        scale_edge_color_manual(values=PAL_GENE_TYPE) +
        scale_edge_width_continuous(range = c(0.2, 1)) +
        theme_void() +
        theme(aspect.ratio=1)
    
    return(plts)
}

make_plots = function(cancer_program_activity, enrichments, networks_sf_ex, shortest_paths_pert_sfs, ppi_network, splicing_factors, cancer_program){
    plts = list(
        plot_program_activity(cancer_program_activity),
        plot_enrichments(enrichments),
        plot_sf_network_analysis(networks_sf_ex, splicing_factors, cancer_program),
        plot_ppi_network_analysis(shortest_paths_pert_sfs, cancer_program_activity, splicing_factors, cancer_program, ppi_network)
    )
    plts = do.call(c,plts)
    return(plts)
}

make_figdata = function(cancer_program_activity, enrichments,  networks_sf_ex, shortest_paths_pert_sfs, ppi_network, splicing_factors, cancer_program){
    figdata = list(
        "upstream_regulators" = list(
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
    save_plt(plts, "program_activity-diff_vs_pert_efficiency-scatter", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "program_activity-ts_vs_onco_by_cell_line-scatter", '.pdf', figs_dir, width=6, height=4)
    save_plt(plts, "program_activity-diff_ranking_overview-scatter", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "program_activity-diff_ranking_by_cell_line_and_cosmic-scatter", '.pdf', figs_dir, width=8, height=17)
    save_plt(plts, "program_activity-diff_ranking_by_cell_line_and_program-scatter", '.pdf', figs_dir, width=8, height=17)
    save_plt(plts, "program_activity-diff_ranking_by_cell_line_and_sfornot-scatter", '.pdf', figs_dir, width=8, height=12)
    save_plt(plts, "program_activity-diff_vs_cosmic-violin", '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, "program_activity-diff_vs_splicing_program-violin", '.pdf', figs_dir, width=4.5, height=5)
    
    save_plt(plts, "enrichments-activity_diff-dot", '.pdf', figs_dir, width=11, height=5.5)
    
    save_plt(plts, "sf_network_analysis-sf_vs_target_type_freq-bar", '.pdf', figs_dir, width=8, height=9)
    save_plt(plts, "sf_network_analysis-target_type_freq_vs_activity_diff-non_driver-scatter", '.pdf', figs_dir, width=8, height=9)
    save_plt(plts, "sf_network_analysis-regulator_type_perc-graph", '.pdf', figs_dir, width=7, height=8)
    
    save_plt(plts, "ppi_network_analysis-pair_type_vs_path_length-violin", '.pdf', figs_dir, width=6.5, height=4)
    save_plt(plts, "ppi_network_analysis-activity_diff_avg_vs_path_length-scatter", '.pdf', figs_dir, width=6.5, height=4)
    save_plt(plts, "ppi_network_analysis-regulator_type_perc-graph", '.pdf', figs_dir, width=7, height=8)
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
        make_option("--protein_activity_rpe1_file", type="character"),
        make_option("--metadata_rpe1_file", type="character"),
        make_option("--cancer_program_file", type="character"),
        make_option("--gene_info_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--cosmic_genes_file", type="character"),
        make_option("--regulons_dir", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--shortest_paths_pert_sfs_file", type="character"),
        make_option("--shortest_paths_random_file", type="character"),
        make_option("--ppi_network_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    protein_activity_rpe1_file = args[["protein_activity_rpe1_file"]]
    metadata_rpe1_file = args[["metadata_rpe1_file"]]
    cancer_program_file = args[["cancer_program_file"]]
    gene_info_file = args[["gene_info_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    cosmic_genes_file = args[["cosmic_genes_file"]]
    regulons_dir = args[["regulons_dir"]]
    event_info_file = args[["event_info_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    shortest_paths_pert_sfs_file = args[["shortest_paths_pert_sfs_file"]]
    shortest_paths_random_file = args[["shortest_paths_random_file"]]
    ppi_network_file = args[["ppi_network_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    set.seed(RANDOM_SEED)
    
    # load
    protein_activity_rpe1 = read_tsv(protein_activity_rpe1_file)
    metadata_rpe1 = read_tsv(metadata_rpe1_file)
    cancer_program = read_tsv(cancer_program_file)
    gene_info = read_tsv(gene_info_file)
    ontologies = load_ontologies(msigdb_dir, cosmic_genes_file)
    networks_sf_ex = lapply(list.files(regulons_dir, full.names=TRUE), function(regulons_file){
        regulon_id = basename(regulons_file) %>% gsub("-delta_psi.tsv.gz","",.)
        regulons = read_tsv(regulons_file) %>%
            mutate(regulon_id = regulon_id)
        return(regulons)
    }) %>% bind_rows()
    event_info = read_tsv(event_info_file)
    splicing_factors = read_tsv(splicing_factors_file)
    
    shortest_paths_pert_sfs = read_tsv(shortest_paths_pert_sfs_file)
    shortest_paths_random = read_tsv(shortest_paths_random_file)
    ppi_network = read_tsv(ppi_network_file)
    
    # prep
    gene_info = gene_info %>%
        mutate(PERT_GENE = `Approved symbol`) %>%
        dplyr::select(`Ensembl gene ID`, PERT_GENE)

    protein_activity_rpe1 = protein_activity_rpe1 %>%
        pivot_longer(-regulator, names_to="PERT_ID", values_to="activity_rpe1") %>%
        separate(PERT_ID, sep="___", into=c("study","cell_line","PERT_ENSEMBL","PERT_TYPE"), remove=FALSE)
    
    protein_activity = protein_activity_rpe1 %>%
        filter(!if_all(c(activity_rpe1), is.na)) %>%
        left_join(gene_info, by=c("PERT_ENSEMBL"="Ensembl gene ID"))
    
    cancer_program_activity = cancer_program %>%
        left_join(
            protein_activity, 
            by=c("ENSEMBL"="regulator")
        ) %>%
        pivot_longer(c(activity_rpe1), names_to="activity_type", values_to="activity") %>%
        drop_na(activity) %>%
        group_by(PERT_ENSEMBL, PERT_GENE, activity_type, driver_type) %>%
        summarize(activity = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        pivot_wider(id_cols=c("PERT_ENSEMBL","PERT_GENE","activity_type"), names_from="driver_type", values_from="activity") %>%
        mutate(
            activity_diff = `Oncogenic` - `Tumor suppressor`
        )  %>%
        left_join(
            ontologies[["cosmic"]] %>% distinct(gene,cosmic_driver_type), by=c("PERT_GENE"="gene")
        ) %>%
        mutate(
            target_in_cosmic = !is.na(cosmic_driver_type)
        ) %>%
        left_join(
            cancer_program %>% distinct(GENE, ENSEMBL, driver_type), 
            by=c("PERT_ENSEMBL"="ENSEMBL", "PERT_GENE"="GENE")
        ) %>%
        mutate(
            is_sf = PERT_GENE %in% splicing_factors[["GENE"]],
            gene_type = case_when(
                is_sf & is.na(driver_type) ~ "Non-driver SF",
                is_sf & driver_type=="Oncogenic" ~ "Oncogenic",
                is_sf & driver_type=="Tumor suppressor" ~ "Tumor suppressor",
                !is_sf ~ "Not SF",
                !is_sf & !is.na(PERT_GENE) ~ "Not SF"
            )
        ) %>%
        left_join(
            metadata_rpe1 %>%
                distinct(PERT_ENSEMBL, n_cells, pert_efficiency_fc) %>%
                mutate(activity_type = "activity_rpe1"),
            by = c("activity_type","PERT_ENSEMBL")
        ) %>%
        mutate(
            pert_is_sf = PERT_ENSEMBL %in% splicing_factors[["ENSEMBL"]]
        )
    
    # enrichments
    ## full
    enrichments = list()
    genes = cancer_program_activity %>% 
        filter(activity_type=="activity_rpe1") %>% 
        arrange(-activity_diff) %>% 
        distinct(PERT_GENE,activity_diff) %>% 
        deframe()
    enrichments[["rpe1_activity_diff_gobp"]] = GSEA(geneList = genes, TERM2GENE=ontologies[["GO_BP"]], seed=RANDOM_SEED)
    enrichments[["rpe1_activity_diff_reactome"]] = GSEA(geneList = genes, TERM2GENE=ontologies[["reactome"]], seed=RANDOM_SEED)
    
    # target exon network analysis
    networks_sf_ex = networks_sf_ex %>% 
        distinct(regulator, target) %>%
        left_join(
            event_info %>% distinct(EVENT,GENE),
            by = c("target"="EVENT")
        ) %>%
        left_join(
            cancer_program %>% 
                distinct(driver_type, GENE) %>%
                rename(driver_type_target=driver_type),
            by = "GENE"
        ) %>%
        left_join(
            cancer_program %>% 
                distinct(driver_type, ENSEMBL) %>%
                rename(driver_type_regulator=driver_type),
            by = c("regulator"="ENSEMBL")
        ) %>%
        mutate(
            target_is_sf = GENE %in% splicing_factors[["GENE"]],
            target_type = case_when(
                target_is_sf & is.na(driver_type_target) ~ "Non-driver SF",
                target_is_sf & driver_type_target=="Oncogenic" ~ "Oncogenic",
                target_is_sf & driver_type_target=="Tumor suppressor" ~ "Tumor suppressor",
                !target_is_sf ~ "Not SF"
            ),
            regulator_is_sf = regulator %in% splicing_factors[["ENSEMBL"]],
            regulator_type = case_when(
                regulator_is_sf & is.na(driver_type_regulator) ~ "Non-driver SF",
                regulator_is_sf & driver_type_regulator=="Oncogenic" ~ "Oncogenic",
                regulator_is_sf & driver_type_regulator=="Tumor suppressor" ~ "Tumor suppressor",
                !regulator_is_sf ~ "Not SF"
            )
        ) %>%
        left_join(
            cancer_program_activity %>% 
                distinct(PERT_ENSEMBL,PERT_GENE,activity_type,activity_diff) %>%
                filter(activity_type=="activity_rpe1"),
            by = c("regulator"="PERT_ENSEMBL")
        ) %>%
        left_join(
            splicing_factors,
            by = c("PERT_GENE"="GENE", "regulator"="ENSEMBL")
        )
    
    # SF PPI analysis
    ## are SFs closer to each other that produce strong changes?
    sf_switch_type = cancer_program_activity %>% 
        filter(activity_type=="activity_rpe1" & pert_is_sf) %>% 
        mutate(
            switch_type = case_when(
                activity_diff >= 2 ~ "strong",
                activity_diff < 2 ~ "weak"
            ) 
        ) %>%
        distinct(PERT_GENE, switch_type)
    
    shortest_paths_pert_sfs = shortest_paths_pert_sfs %>%
        left_join(
            sf_switch_type %>% rename(switch_type_source = switch_type),
            by=c("source"="PERT_GENE")
        ) %>%
        left_join(
            sf_switch_type %>% rename(switch_type_target = switch_type),
            by=c("target"="PERT_GENE")
        ) %>%
        mutate(
            pair_type = case_when(
                switch_type_source=="weak" & switch_type_target=="weak" ~ "Both Weak",
                switch_type_source=="weak" & switch_type_target=="strong" ~ "Mixed",
                switch_type_source=="strong" & switch_type_target=="weak" ~ "Mixed",
                switch_type_source=="strong" & switch_type_target=="strong" ~ "Both Strong"
            )
        )
    
    # plot
    plts = make_plots(cancer_program_activity, enrichments, networks_sf_ex, shortest_paths_pert_sfs, ppi_network, splicing_factors, cancer_program)
    
    # make figdata
    figdata = make_figdata(cancer_program_activity, enrichments, networks_sf_ex, shortest_paths_pert_sfs, ppi_network, splicing_factors, cancer_program)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}

##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
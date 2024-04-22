#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------

Sys.setenv(VROOM_CONNECTION_SIZE='100000000')
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

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","program_regulation")
# protein_activity_rpe1_file = file.path(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_rpe1-genexpr.tsv.gz")
# protein_activity_k562_file = file.path(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_K562_essential-genexpr.tsv.gz")
# gene_info_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# cancer_program_file = file.path(SUPPORT_DIR,"supplementary_tables","cancer_program.tsv.gz")
# msigdb_dir = file.path(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# cosmic_genes_file = file.path(RAW_DIR,"COSMIC","cancer_gene_census.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","upstream_regulators")

##### FUNCTIONS #####
load_ontologies = function(msigdb_dir, cosmic_genes_file){
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "GO_CC" = read.gmt(file.path(msigdb_dir,"c5.go.cc.v7.4.symbols.gmt")),
        "cosmic" = read_tsv(cosmic_genes_file) %>%
            dplyr::select("Gene Symbol") %>%
            rename(gene = `Gene Symbol`) %>%
            mutate(term = "COSMIC_CENSUS") %>%
            dplyr::select(term,gene)
    )
    return(ontologies)
}


plot_program_activity = function(cancer_program_activity){
    plts = list()
    
    X = cancer_program_activity 
    x = X %>%
        pivot_wider(
            id_cols=c("PERT_ENSEMBL","PERT_GENE","in_cosmic"), 
            names_from="activity_type", values_from="activity_diff"
        )
    
    plts[["program_activity-ts_vs_onco_by_cell_line-scatter"]] = X %>%
        ggplot(aes(x=`Tumor suppressor`, y=`Oncogenic`)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        theme_pubr() + 
        geom_abline(intercept=0, slope=1, size=LINE_SIZE , linetype="dashed", color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~activity_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Tumor Suppressor Splicing Program Activity", y="Oncogenic Splicing Program Activity")
    
    plts[["program_activity-diff_ranking_by_cell_line-scatter"]] = X %>%
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
            group_by(activity_type, in_cosmic) %>% 
            slice_max(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(activity_type, in_cosmic) %>% 
            slice_min(activity_diff, n=5) %>% 
            ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        facet_wrap(~activity_type+in_cosmic) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Relative Ranking", y="Oncogenic vs Tumor Suppressor\nSplicing Program Activity")
    
    plts[["program_activity-diff_rpe1_vs_k562-scatter"]] = x %>%
        ggplot(aes(x=activity_rpe1, y=activity_k562)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8, alpha=0.5, color=PAL_DARK) +
        theme_pubr() + 
        geom_hline(yintercept=0, size=LINE_SIZE , linetype="dashed", color="black") +
        geom_vline(xintercept=0, size=LINE_SIZE , linetype="dashed", color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~in_cosmic) +
        geom_text_repel(
            aes(label=PERT_GENE),
            . %>% 
            group_by(in_cosmic) %>% 
            slice_max(abs(abs(activity_rpe1) - abs(activity_k562)), n=6) %>% 
            ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="RPE1 Cancer Splicing Programs Activity Difference", y="K562 Cancer Splicing Programs Activity Difference")
    
    return(plts)
}


plot_enrichments = function(enrichments){
    plts = list()
    
    X = enrichments
    for (x in names(X)){
        plt_name = sprintf("enrichments-activity_diff-%s-dot", x)
        plts[[plt_name]] = X[[x]] %>%
            dotplot() + 
            scale_size(range=c(0.5,3)) + 
            scale_color_continuous(
                low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
                name="FDR", guide=guide_colorbar(reverse=TRUE)) +
            theme_pubr()
            
    }
            
    return(plts)
}

make_plots = function(cancer_program_activity, enrichments){
    plts = list(
        plot_program_activity(cancer_program_activity),
        plot_enrichments(enrichments)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(cancer_program_activity, enrichments){
    figdata = list(
        "eval_bulk_vs_singlecell" = list(
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
    save_plt(plts, "program_activity-ts_vs_onco_by_cell_line-scatter", '.pdf', figs_dir, width=6, height=4)
    save_plt(plts, "program_activity-diff_ranking_by_cell_line-scatter", '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, "program_activity-diff_rpe1_vs_k562-scatter", '.pdf', figs_dir, width=6, height=4)
    save_plt(plts, "enrichments-activity_diff-k562-dot", '.pdf', figs_dir, width=8.5, height=6)
    save_plt(plts, "enrichments-activity_diff-rpe1-dot", '.pdf', figs_dir, width=13, height=6)
    save_plt(plts, "enrichments-activity_diff-k562_vs_rpe1-dot", '.pdf', figs_dir, width=11.5, height=6)
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
        make_option("--protein_activity_k562_file", type="character"),
        make_option("--cancer_program_file", type="character"),
        make_option("--gene_info_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--cosmic_genes_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    protein_activity_rpe1_file = args[["protein_activity_rpe1_file"]]
    protein_activity_k562_file = args[["protein_activity_k562_file"]]
    cancer_program_file = args[["cancer_program_file"]]
    gene_info_file = args[["gene_info_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    cosmic_genes_file = args[["cosmic_genes_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity_rpe1 = read_tsv(protein_activity_rpe1_file)
    protein_activity_k562 = read_tsv(protein_activity_k562_file)
    cancer_program = read_tsv(cancer_program_file)
    gene_info = read_tsv(gene_info_file)
    ontologies = load_ontologies(msigdb_dir, cosmic_genes_file)
    
    # prep
    gene_info = gene_info %>%
        mutate(PERT_GENE = `Approved symbol`) %>%
        dplyr::select(`Ensembl gene ID`, PERT_GENE)

    protein_activity_k562 = protein_activity_k562 %>%
        pivot_longer(-regulator, names_to="condition", values_to="activity_k562") %>%
        left_join(metadata_k562, by="condition") %>%
        left_join(gene_info, by=c("PERT_ENSEMBL"="Ensembl gene ID"))
    
    cancer_program_k562 = cancer_program %>%
        left_join(
            protein_activity_k562, 
            by=c("ENSEMBL"="regulator")
        ) %>%
        group_by(driver_type, PERT_ENSEMBL, PERT_GENE, pert_efficiency_fc, condition, n_cells) %>%
        summarize(activity_k562 = median(activity_k562, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(
            target_in_cosmic = ifelse(
                PERT_GENE %in% ontologies[["cosmic"]][["gene"]], "In COSMIC", "Not in COSMIC"
            )
        ) 
    
    protein_activity_rpe1 = protein_activity_rpe1 %>%
        pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_rpe1")
    
    protein_activity_k562 = protein_activity_k562 %>%
            pivot_longer(-regulator, names_to="PERT_ENSEMBL", values_to="activity_k562")
    
    protein_activity = merge(
            protein_activity_rpe1, protein_activity_k562, 
            by=c("PERT_ENSEMBL","regulator"), all.x=TRUE, all.y=TRUE
        ) %>%
        filter(!if_all(c(activity_rpe1,activity_k562), is.na)) %>%
        left_join(gene_info, by=c("PERT_ENSEMBL"="Ensembl gene ID"))
    
    cancer_program_activity = cancer_program %>%
        left_join(
            protein_activity, 
            by=c("ENSEMBL"="regulator")
        ) %>%
        pivot_longer(c(activity_rpe1, activity_k562), names_to="activity_type", values_to="activity") %>%
        drop_na(activity) %>%
        group_by(PERT_ENSEMBL, PERT_GENE, activity_type, driver_type) %>%
        summarize(activity = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        pivot_wider(id_cols=c("PERT_ENSEMBL","PERT_GENE","activity_type"), names_from="driver_type", values_from="activity") %>%
        mutate(
            activity_diff = `Oncogenic` - `Tumor suppressor`,
            in_cosmic = ifelse(PERT_GENE %in% ontologies[["cosmic"]][["gene"]], "In COSMIC", "Not in COSMIC")
        ) 
    
    # enrichments
    enrichments = list()
    genes = cancer_program_activity %>% 
        filter(activity_type=="activity_k562") %>% 
        arrange(-activity_diff) %>% 
        distinct(PERT_GENE,activity_diff) %>% 
        deframe()
    enrichments[["k562"]] = GSEA(geneList = genes, TERM2GENE=ontologies[["GO_BP"]])
    genes = cancer_program_activity %>% 
        filter(activity_type=="activity_rpe1") %>% 
        arrange(-activity_diff) %>% 
        distinct(PERT_GENE,activity_diff) %>% 
        deframe()
    enrichments[["rpe1"]] = GSEA(geneList = genes, TERM2GENE=ontologies[["GO_BP"]])
    genes = cancer_program_activity %>% 
        pivot_wider(id_cols=c("PERT_ENSEMBL","PERT_GENE"), names_from="activity_type", values_from="activity_diff") %>%
        drop_na() %>%
        mutate(activity_diff_diff = activity_k562 - activity_rpe1) %>%
        arrange(-activity_diff_diff) %>% 
        distinct(PERT_GENE,activity_diff_diff) %>% 
        deframe()
    enrichments[["k562_vs_rpe1"]] = GSEA(geneList = genes, TERM2GENE=ontologies[["GO_BP"]])
    
    # plot
    plts = make_plots(cancer_program_activity, enrichments)
    
    # make figdata
    figdata = make_figdata(cancer_program_activity, enrichments)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(ggrepel)
require(clusterProfiler)

# variables
RANDOM_SEED = 1234
THRESH_FDR = 0.05

FIBROBLASTS = c("BJ_PRIMARY","BJ_IMMORTALIZED","BJ_TRANSFORMED","BJ_METASTATIC")

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DRIVER_TYPE = c(
    #"Non-driver"="lightgrey",
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

PAL_GENE_TYPE = c(
    "Not SF"="darkgreen",
    "Non-driver SF"="darkred",
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

PAL_DARK = "darkred"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","new_empirical_network")
# genexpr_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets","genexpr_tpm","tumorigenesis.tsv.gz")
# splicing_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets","event_psi","tumorigenesis-EX.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","carcinogenesis-EX.tsv.gz")
# metadata_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets","metadata","tumorigenesis.tsv.gz")
# driver_types_file = file.path(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
# event_info_file = file.path(SUPPORT_DIR,"supplementary_tables","supdata01_event_prior_knowledge.txt")
# regulons_dir = file.path(ROOT,"results","new_empirical_network","files","experimentally_derived_regulons_pruned_w_viper_networks-EX")
# msigdb_dir = file.path(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# figs_dir = file.path(RESULTS_DIR,"figures","carcinogenesis")

##### FUNCTIONS #####
load_ontologies = function(msigdb_dir, cosmic_genes_file){
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "GO_CC" = read.gmt(file.path(msigdb_dir,"c5.go.cc.v7.4.symbols.gmt"))
    )
    return(ontologies)
}

plot_carcinogenesis = function(protein_activity, genexpr, splicing, program_networks){
    plts = list()
    
    X = protein_activity %>%
        drop_na(driver_type) %>%
        group_by(cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        mutate(
            cell_line_name = factor(cell_line_name, levels=FIBROBLASTS),
            driver_type = factor(driver_type, levels=names(PAL_DRIVER_TYPE))
        )
    
    plts[["carcinogenesis-cell_line_vs_activity-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=activity, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_text(
            aes(y=-3, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Cell Line", y="Protein Activity", fill="Driver Type")
    
    plts[["carcinogenesis-cell_line_vs_activity_diff-line"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
        group_by(cell_line_name, study_accession, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(cell_line_name, study_accession) %>%
        summarize(activity_diff = sum(activity)) %>%
        ungroup() %>%
        ggline(
            x="cell_line_name", y="activity_diff", color=PAL_DARK, numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="dashed", point.size=0.05
        ) +
        geom_hline(yintercept=0, color="black", linetype="dashed", linewidth=LINE_SIZE) +
        labs(x="Cell Line", y="Protein Activity Diff.")

    # genexpr
    X = genexpr %>%
        drop_na(driver_type) %>%
        filter(GENE %in% X[["GENE"]]) %>%
        group_by(cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(genexpr_tpm = median(genexpr_tpm)) %>%
        ungroup() %>%
        mutate(cell_line_name=factor(
            cell_line_name, levels=FIBROBLASTS
        ))
    ctl_genexpr = genexpr %>%
        filter(cell_line_name=="BJ_PRIMARY") %>%
        distinct(genexpr_tpm, GENE) %>%
        dplyr::rename(ctl_tpm = genexpr_tpm)
    X = X %>%
        left_join(ctl_genexpr, by="GENE") %>%
        mutate(genexpr_tpm_fc = genexpr_tpm - ctl_tpm)

    plts[["carcinogenesis-cell_line_vs_genexpr_fc-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=genexpr_tpm_fc, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_text(
            aes(y=-1.5, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Cell Line", y="Gene Expression log2FC", fill="Driver Type")
    
    # splicing of cancer driver exons
    X = splicing %>%
        mutate(
            cell_line_name = factor(cell_line_name, levels=FIBROBLASTS),
            is_known_driver = replace_na(is_known_driver, FALSE)
        ) %>%
        drop_na(cell_line_name)
    ctl_splicing = splicing %>%
        filter(cell_line_name=="BJ_PRIMARY") %>%
        distinct(event_psi, EVENT) %>%
        dplyr::rename(ctl_psi = event_psi)
    X = X %>%
        left_join(ctl_splicing, by="EVENT") %>%
        mutate(
            event_dpsi = event_psi - ctl_psi,
            abs_event_dpsi = abs(event_dpsi),
        )
    
    plts[["carcinogenesis-cell_line_vs_splicing_dpsi-driver_exons-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=abs_event_dpsi, group=interaction(cell_line_name,is_known_driver))) +
        geom_boxplot(aes(color=is_known_driver), outlier.size=0.1, fill=NA, position=position_dodge(0.9)) +
        color_palette(c("lightgrey",PAL_DARK)) +
        geom_text(
            aes(y=-3, label=label),
            . %>% count(cell_line_name, is_known_driver) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Carcinogenic Stage", y="|Delta PSI|", color="Cancer Driver Exon")
    
    X = program_networks %>%
        distinct(regulator, target, GENE, driver_type, is_known_driver) %>%
        left_join(X %>% distinct(EVENT, cell_line_name, event_dpsi, abs_event_dpsi), by=c("target"="EVENT")) %>%
        mutate(
            driver_type = replace_na(driver_type, "Non-driver SF"),
            is_known_driver = replace_na(is_known_driver, FALSE)            
        )
    
    plts[["carcinogenesis-cell_line_vs_splicing_dpsi-driver_programs_targets-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=abs_event_dpsi, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), trim=TRUE, color=NA) +
        geom_boxplot(width=0.2, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_GENE_TYPE) +
        geom_text(
            aes(y=-3, label=label),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Carcinogenic Stage", y="|Delta PSI|", fill="Regulator SF Type")
    
    # random control
    X = protein_activity %>%
        group_by(cell_line_name, study_accession) %>%
        mutate(driver_type = sample(driver_type)) %>%
        ungroup() %>%
        drop_na(driver_type) %>%
        group_by(cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        filter(study_accession=="PRJNA193487") %>%
        mutate(cell_line_name=factor(
            cell_line_name, levels=FIBROBLASTS
        ))
    
    plts[["carcinogenesis-cell_line_vs_activity-random-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=activity, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_text(
            aes(y=-3, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Cell Line", y="Protein Activity", fill="Driver Type")
    
    return(plts)
}

plot_enrichments = function(enrichments){
    plts = list()
    
    X = enrichments
    
    terms_oi = X %>%
        group_by(gene_set) %>%
        slice_max(GeneRatio, n=10) %>%
        ungroup() %>%
        pull(Description) %>%
        unique()
    
    plts[["enrichments-reactome-bar"]] = X %>%
        filter(Description %in% terms_oi) %>%
        group_by(Description) %>%
        mutate(ratio_sums = sum(GeneRatio)) %>%
        ungroup() %>%
        arrange(GeneRatio) %>%
        ggbarplot(x="Description", y="GeneRatio", fill="driver_type", color=NA,
                  palette=PAL_DRIVER_TYPE, position=position_dodge(0.9)) +
        geom_text(aes(label=Count, group=driver_type), 
                  size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9), hjust=-0.1) +
        labs(x="Description", y="GeneRatio", fill="Driver Type") +
        coord_flip()
    
    return(plts)
}

make_plots = function(
    protein_activity, genexpr, splicing, enrichments, program_networks
){
    plts = list(
        plot_carcinogenesis(protein_activity, genexpr, splicing, program_networks),
        plot_enrichments(enrichments)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    protein_activity
){
    figdata = list(
        "carcinogenesis" = list(
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


save_plots = function(plts, figs_dir){
    save_plt(plts, "carcinogenesis-cell_line_vs_activity-violin", '.pdf', figs_dir, width=6, height=4.5)
    save_plt(plts, "carcinogenesis-cell_line_vs_activity_diff-line", '.pdf', figs_dir, width=6, height=2.25)
    save_plt(plts, "carcinogenesis-cell_line_vs_genexpr_fc-violin", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "carcinogenesis-cell_line_vs_splicing_dpsi-driver_exons-violin", '.pdf', figs_dir, width=6, height=4.5)
    save_plt(plts, "carcinogenesis-cell_line_vs_splicing_dpsi-driver_programs_targets-violin", '.pdf', figs_dir, width=6, height=5.5)
    save_plt(plts, "carcinogenesis-cell_line_vs_activity-random-violin", '.pdf', figs_dir, width=6, height=6)
    
    save_plt(plts, "enrichments-reactome-bar", '.pdf', figs_dir, width=16, height=8)
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
        make_option("--genexpr_file", type="character"),
        make_option("--protein_activity_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--driver_types_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--regulons_dir", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    genexpr_file = args[["genexpr_file"]]
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    driver_types_file = args[["driver_types_file"]]
    event_info_file = args[["event_info_file"]]
    regulons_dir = args[["regulons_dir"]]
    msigdb_dir = args[["msigdb_dir"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    set.seed(RANDOM_SEED)    
    
    # load
    genexpr = read_tsv(genexpr_file)
    splicing = read_tsv(splicing_file)
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    driver_types = read_tsv(driver_types_file)
    event_info = read_tsv(event_info_file)
    networks_sf_ex = lapply(list.files(regulons_dir, full.names=TRUE), function(regulons_file){
        regulon_id = basename(regulons_file) %>% gsub("-delta_psi.tsv.gz","",.)
        regulons = read_tsv(regulons_file) %>%
            mutate(regulon_id = regulon_id)
        return(regulons)
    }) %>% bind_rows()
    ontologies = load_ontologies(msigdb_dir)
    
    # prep
    protein_activity = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, activity) %>%
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
            activity = median(activity, na.rm=TRUE),
            abs_activity = abs(activity),
        ) %>%
        ungroup() %>%
        
        # add activity
        group_by(condition_lab) %>%
        arrange(activity) %>%
        mutate(
            activity_ranking = row_number(),
            total_avail_sfs = sum(regulator %in% unlist(strsplit(PERT_ENSEMBL, ",")))
        ) %>%
        arrange(abs_activity) %>%
        mutate(
            abs_activity_ranking = row_number(),
        ) %>%
        ungroup() %>%
        left_join(driver_types, by=c("regulator"="ENSEMBL"))
    
    genexpr = genexpr %>%
        filter(ID%in%driver_types[["ENSEMBL"]]) %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr_tpm") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, genexpr_tpm) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s%s) (%s%s) | %s | %s", condition, pert_time, pert_time_units, 
                pert_concentration, pert_concentration_units, cell_line_name, study_accession
            )
        ) %>%
        
        # summarize replicates
        group_by(condition_lab, condition, pert_time, pert_time_units, 
                 pert_concentration, pert_concentration_units, cell_line_name, study_accession,
                 PERT_ENSEMBL, PERT_GENE, ID) %>%
        summarize(
            genexpr_tpm = median(genexpr_tpm, na.rm=TRUE),
        ) %>%
        ungroup() %>%
        left_join(driver_types, by=c("ID"="ENSEMBL")) %>%
        drop_na(driver_type)
    
    splicing = splicing %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="event_psi") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, event_psi) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s%s) (%s%s) | %s | %s", condition, pert_time, pert_time_units, 
                pert_concentration, pert_concentration_units, cell_line_name, study_accession
            )
        ) %>%
        # summarize replicates
        group_by(condition_lab, condition, pert_time, pert_time_units, 
                 pert_concentration, pert_concentration_units, cell_line_name, study_accession,
                 PERT_ENSEMBL, PERT_GENE, EVENT) %>%
        summarize(
            event_psi = median(event_psi, na.rm=TRUE),
        ) %>%
        ungroup() %>%
        left_join(event_info, by="EVENT")

    # enrichment program targets
    program_networks = networks_sf_ex %>%
        distinct(regulator, target) %>%
        left_join(event_info, by=c("target"="EVENT")) %>%
        left_join(driver_types %>% distinct(ENSEMBL,driver_type), by=c("regulator"="ENSEMBL"))
    
    enrichments = list()
    genes = program_networks %>% filter(driver_type=="Oncogenic") %>% pull(GENE) %>% unique()
    enrichments[["oncogenics"]] = enricher(gene=genes, TERM2GENE=ontologies[["reactome"]])
    genes = program_networks %>% filter(driver_type=="Tumor suppressor") %>% pull(GENE) %>% unique()
    enrichments[["suppressors"]] = enricher(gene=genes, TERM2GENE=ontologies[["reactome"]])
    enrichments = lapply(names(enrichments), function(gene_set_oi){
            x = enrichments[[gene_set_oi]] %>%
                as.data.frame() %>%
                mutate(gene_set = gene_set_oi)
            return(x)
        }) %>%
        bind_rows() %>%
        filter(p.adjust < THRESH_FDR) %>%
        rowwise() %>%
        mutate(GeneRatio = eval(parse(text=GeneRatio))) %>%
        ungroup() %>%
        mutate(
            driver_type = case_when(
                gene_set=="suppressors" ~ "Tumor suppressor",
                gene_set=="oncogenics" ~ "Oncogenic"
            )
        )
    
    # plot
    plts = make_plots(protein_activity, genexpr, splicing, enrichments, program_networks)
    
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
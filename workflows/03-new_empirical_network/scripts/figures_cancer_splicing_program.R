require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(survival)
require(survminer)

# variables
THRESH_FDR = 0.05
THRESH_N_SUM = 5.5

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DRIVER_TYPE = c(
    "Random Genes"="darkgreen",
    "Not Driver-like"="grey",
    "Oncogenic"="#F6AE2D",
    "Tumor suppressor"="#6C98B3"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","new_empirical_network")
# diff_activity_file = file.path(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# splicing_factors_file = file.path(RESULTS_DIR,"figures","empirical_psi_networks","figdata","eda","splicing_factors.tsv.gz")
# gene_annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")

# program_activity_diff_file = file.path(RESULTS_DIR,"files","protein_activity",'PANCAN-PrimaryTumor_vs_SolidTissueNormal-program_activity_diff.tsv.gz')
# metadata_file = file.path(RAW_DIR,'UCSCXena','TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.tsv')
# mutations_file = file.path(RAW_DIR,'UCSCXena','TCGA','snv','mc3.v0.2.8.PUBLIC.xena.gz')
# driver_types_file = file.path(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')

# figs_dir = file.path(RESULTS_DIR,"figures","cancer_splicing_program")

##### FUNCTIONS #####
plot_driver_selection = function(driver_activity){
    plts = list()
    
    # SF activity
    X = driver_activity
    x = X %>%
        count(GENE, driver_type, only_new_network) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup()
    
    sf_order = x %>%
        pivot_wider(id_cols=c("GENE","only_new_network"), names_from="driver_type", values_from="n_sign", values_fill=0) %>%
        rowwise() %>%
        mutate(diff=sum(`Tumor suppressor`,`Oncogenic`)) %>%
        ungroup() %>%
        arrange(diff, only_new_network,`Tumor suppressor`,`Oncogenic`) %>%
        pull(GENE)
    
    sfs_oi = x %>% 
        slice_max(n_sum, n=5) %>% 
        filter(driver_type=="Oncogenic") %>%
        bind_rows(
            x %>% 
            slice_min(n_sum, n=5) %>% 
            filter(driver_type=="Tumor suppressor")
        )

    plts[["driver_selection-n_signif_vs_driver_type-activity-bar"]] = x %>%
        mutate(GENE = factor(GENE, levels=sf_order)) %>%
        ggbarplot(x="GENE", y="n", fill="only_new_network", color=NA) +
        geom_text(
            aes(label=GENE, color=only_new_network),
            sfs_oi %>% filter(driver_type=="Oncogenic"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=1, vjust=1, nudge_y=0.25
        ) +
        geom_text(
            aes(label=GENE, color=only_new_network),
            sfs_oi %>% filter(driver_type=="Tumor suppressor"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=45, hjust=0, vjust=0.5, nudge_y=0.25
        ) +
        geom_hline(yintercept=THRESH_N_SUM, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_text(
            aes(x=x, label=label),
            . %>% 
            filter(abs(n_sum)>THRESH_N_SUM) %>%
            group_by(GENE) %>%
            slice_max(n, n=1) %>%
            ungroup() %>%
            count(driver_type) %>% 
                mutate(
                    label=paste0("n=",n),
                    n=c(15,15),
                    x=c(120,40)
                ),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette("Dark2") +
        color_palette("Dark2") +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~driver_type, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Splicing Factor", y="Count", fill="SF from new Version", color="SF from new Version")
    
    return(plts)
}

plot_survival_analysis = function(program_activity_diff, survival_analysis){
    plts = list()
    
    # distributions of program activity differences in primary tumors
    X = program_activity_diff
    plts[["survival_analysis-cancer_vs_activity_diff-violin"]] = X %>%
        ggviolin(x="cancer_type", y="activity_diff", color=NA, fill="cancer_type", trim=TRUE) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        geom_text(
            aes(y=-2.5, label=label),
            . %>% count(cancer_type) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette(get_palette("Paired", 14)) +
        theme_pubr(x.text.angle=45) +
        guides(fill="none") +
        labs(x="Cancer Type", y="Median Program Act. Diff.")
    
    # associations with survival
    X = survival_analysis
    plts[["survival_analysis-cancer_vs_coxph-bar"]] = X %>%
        mutate(
            surv_metric = factor(surv_metric, levels=c("OS","PFI","DFI"))
        ) %>%
        ggbarplot(x="cancer_type", y="coxph_coef", color=NA, fill="cancer_type") +
        geom_text(
            aes(y=1.5, label=label),
            . %>% mutate(label=sprintf("p=%s", scales::pvalue(coxph_pvalue))),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette(get_palette("Paired", 14)) +
        facet_wrap(~surv_metric, nrow=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Cancer Type", y="Cox PH Coeff.") +
        coord_flip()
    
    # kaplan meier for interesting cases?
    ## KICH
    cancer_type_oi = "KICH"
    surv_metric = "OS"
    surv_time_col = sprintf("%s.time", surv_metric)
    surv_event_col = surv_metric
    vars_formula = "activity_diff"
    result_cut = surv_cutpoint(
        program_activity_diff %>% filter(cancer_type==cancer_type_oi), 
        time=surv_time_col, event=surv_event_col, variables="activity_diff"
    )
    result_cat = surv_categorize(result_cut) %>% data.frame()
    fit = surv_fit(
        as.formula(sprintf("Surv(%s, %s) ~ %s", surv_time_col, surv_event_col, vars_formula)), 
        data=result_cat
    )
    plts[["survival_analysis-KICH-km"]] = ggsurvplot(
        fit, data=result_cat, 
        risk.table=TRUE, conf.int=TRUE, pval=TRUE, pval.size=FONT_SIZE+2,
        risk.table.fontsize=FONT_SIZE+2, risk.table.font.family=FONT_FAMILY,
        palette = get_palette("Dark2", 2)
    )
    
    ## BRCA
    cancer_type_oi = "BRCA"
    result_cut = surv_cutpoint(
        program_activity_diff %>% filter(cancer_type==cancer_type_oi), 
        time=surv_time_col, event=surv_event_col, variables="activity_diff"
    )
    result_cat = surv_categorize(result_cut) %>% data.frame()
    fit = surv_fit(
        as.formula(sprintf("Surv(%s, %s) ~ %s", surv_time_col, surv_event_col, vars_formula)), 
        data=result_cat
    )
    plts[["survival_analysis-BRCA-km"]] = ggsurvplot(
        fit, data=result_cat, 
        risk.table=TRUE, conf.int=TRUE, pval=TRUE, pval.size=FONT_SIZE+2,
        risk.table.fontsize=FONT_SIZE+2, risk.table.font.family=FONT_FAMILY,
        palette = get_palette("Dark2", 2)
    )
    
    return(plts)
}


plot_mutation_analysis = function(program_activity_diff, mutation_analysis){
    plts = list()
    
    # distributions on samples with a mutated splicing factor vs WT
    X = program_activity_diff %>%
        distinct(index, activity_diff, cancer_type, sf_mutated, mutated_sf3b1, mutated_u2af1, mutated_driver_type)
    
    plts[["mutation_analysis-sf_mutated-overall-violin"]] = X %>%
        ggviolin(x="sf_mutated", y="activity_diff", fill="sf_mutated", color=NA, trim=TRUE) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("Dark2") +
        stat_compare_means(
            method="wilcox", size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(sf_mutated) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="SF(s) with Non-silent Mutation", y="Median Program Act. Diff.")
    
    plts[["mutation_analysis-sf_mutated-by_cancer_type-violin"]] = X %>%
        ggplot(aes(x=cancer_type, y=activity_diff, group=interaction(cancer_type, sf_mutated))) +
        geom_violin(aes(fill=sf_mutated), color=NA, trim=TRUE, position=position_dodge(0.9)) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA, position=position_dodge(0.9)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("Dark2") +
        stat_compare_means(
            method="wilcox", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=45
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(cancer_type, sf_mutated) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr(x.text.angle=45) +
        labs(x="Cancer Cohort", y="Median Program Act. Diff.", fill="SF(s) with Non-silent Mutation")
    
    plts[["mutation_analysis-mutated_sf3b1-overall-violin"]] = X %>%
        ggviolin(x="mutated_sf3b1", y="activity_diff", fill="mutated_sf3b1", color=NA, trim=TRUE) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("jco") +
        stat_compare_means(
            method="wilcox", size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(mutated_sf3b1) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="SF3B1 with Non-silent Mutation", y="Median Program Act. Diff.")    
    
    plts[["mutation_analysis-mutated_sf3b1-by_cancer_type-violin"]] = X %>%
        ggplot(aes(x=cancer_type, y=activity_diff, group=interaction(cancer_type, mutated_sf3b1))) +
        geom_violin(aes(fill=mutated_sf3b1), color=NA, trim=TRUE, position=position_dodge(0.9)) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA, position=position_dodge(0.9)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("jco") +
        stat_compare_means(
            method="wilcox", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=45
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(cancer_type, mutated_sf3b1) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr(x.text.angle=45) +
        labs(x="Cancer Cohort", y="Median Program Act. Diff.", fill="SF3B1 with Non-silent Mutation")
    
    plts[["mutation_analysis-mutated_u2af1-overall-violin"]] = X %>%
        ggviolin(x="mutated_u2af1", y="activity_diff", fill="mutated_u2af1", color=NA, trim=TRUE) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("jco") +
        stat_compare_means(
            method="wilcox", size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(mutated_u2af1) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="U2AF1 with Non-silent Mutation", y="Median Program Act. Diff.")    
    
    plts[["mutation_analysis-mutated_u2af1-by_cancer_type-violin"]] = X %>%
        ggplot(aes(x=cancer_type, y=activity_diff, group=interaction(cancer_type, mutated_u2af1))) +
        geom_violin(aes(fill=mutated_u2af1), color=NA, trim=TRUE, position=position_dodge(0.9)) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA, position=position_dodge(0.9)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("jco") +
        stat_compare_means(
            method="wilcox", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=45
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(cancer_type, mutated_u2af1) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr(x.text.angle=45) +
        labs(x="Cancer Cohort", y="Median Program Act. Diff.", fill="U2AF1 with Non-silent Mutation")
    
    plts[["mutation_analysis-mutated_driver_type-overall-violin"]] = X %>%
        mutate(mutated_driver_type = factor(mutated_driver_type, levels=c("WT SFs","Non-driver-like","Onco only","TS only","OncoTS"))) %>%
        ggviolin(x="mutated_driver_type", y="activity_diff", fill="mutated_driver_type", color=NA, trim=TRUE) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("futurama") +
        stat_compare_means(
            ref.group="WT SFs", method="wilcox", label="p.format", size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(mutated_driver_type) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        labs(x="SF Class with Non-silent Mutation", y="Median Program Act. Diff.")
    
    plts[["mutation_analysis-mutated_driver_type-by_cancer_type-violin"]] = X %>%
        mutate(mutated_driver_type = factor(mutated_driver_type, levels=c("WT SFs","Non-driver-like","Onco only","TS only","OncoTS"))) %>%
        ggviolin(x="mutated_driver_type", y="activity_diff", fill="mutated_driver_type", color=NA, trim=TRUE) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        fill_palette("futurama") +
        stat_compare_means(
            ref.group="WT SFs", method="wilcox", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=45
        ) +
        geom_text(
            aes(y=-2.3, label=label),
            . %>% count(mutated_driver_type) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        facet_wrap(~cancer_type, nrow=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="SF Class with Non-silent Mutation", y="Median Program Act. Diff.", fill="SF Class & Mutation Status")
    
    # distributions of program activity differences between types of mutations
    X = mutation_analysis %>%
        distinct(sample, effect, activity_diff, gene, cancer_type) %>%
        drop_na()
    
    effects_oi = X %>% 
        count(effect, gene) %>%
        filter(n>=5) %>%
        pull(effect) %>%
        unique()
    
    cancers_oi = X %>% 
        count(cancer_type, effect, gene) %>%
        filter(n>=10) %>%
        pull(effect) %>%
        unique()
    
    # distributions of program activity differences
    ## overall - by mutation type and gene
    plts[["mutation_analysis-effects-overall-violin"]] = X %>%
        filter(effect%in%effects_oi) %>%
        ggviolin(x="effect", y="activity_diff", fill="gene", color=NA, trim=TRUE) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        stat_compare_means(
            ref.group="Silent", method="wilcox", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=45
        ) +
        geom_text(
            aes(y=-1, label=label),
            . %>% count(effect, gene) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        facet_wrap(~gene, nrow=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        theme_pubr(x.text.angle=45) +
        guides(fill="none") +
        labs(x="Effect Type", y="Median Program Act. Diff.")
    
    ## by cancer type and mutation type and gene
    plts[["mutation_analysis-effects-by_cancer-violin"]] = X %>%
        filter(effect%in%c("Silent","Missense_Mutation")) %>%
        ggplot(aes(x=cancer_type, y=activity_diff, group=interaction(cancer_type, effect))) +
        geom_violin(aes(fill=effect), color=NA, trim=TRUE, position=position_dodge(0.9)) +
        geom_boxplot(outlier.shape=NA, width=0.1, fill=NA, position=position_dodge(0.9)) +
        fill_palette(get_palette("Paired", 14)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=LINE_SIZE) +
        stat_compare_means(
            method="wilcox", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=45
        ) +
        geom_text(
            aes(y=-1, label=label),
            . %>% count(effect, gene, cancer_type) %>% mutate(label=sprintf("n=%s", n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        facet_wrap(~gene, nrow=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        theme_pubr(x.text.angle=45) +
        labs(x="Effect Type", y="Median Program Act. Diff.")
    
    return(plts)
}


make_plots = function(diff_activity, program_activity_diff, survival_analysis, mutation_analysis){
    plts = list(
        plot_driver_selection(diff_activity),
        plot_survival_analysis(program_activity_diff, survival_analysis),
        plot_mutation_analysis(program_activity_diff, mutation_analysis)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(diff_activity){
    figdata = list(
        "cancer_program" = list(
            "enrichments_reactome" = diff_activity
        )
    )
    return(figdata)
}

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
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
    # cancer splicing program definition
    save_plt(plts, "driver_selection-n_signif_vs_driver_type-activity-bar", '.pdf', figs_dir, width=5, height=7)
    
    save_plt(plts, "survival_analysis-cancer_vs_activity_diff-violin", '.pdf', figs_dir, width=10, height=5)
    save_plt(plts, "survival_analysis-cancer_vs_coxph-bar", '.pdf', figs_dir, width=12, height=5)
    save_plt(plts, "survival_analysis-KICH-km", '.pdf', figs_dir, width=10, height=13, format=FALSE)
    save_plt(plts, "survival_analysis-BRCA-km", '.pdf', figs_dir, width=10, height=13, format=FALSE)

    save_plt(plts, "mutation_analysis-sf_mutated-overall-violin", '.pdf', figs_dir, width=3, height=8)
    save_plt(plts, "mutation_analysis-sf_mutated-by_cancer_type-violin", '.pdf', figs_dir, width=13, height=8)
    save_plt(plts, "mutation_analysis-mutated_sf3b1-overall-violin", '.pdf', figs_dir, width=3, height=8)
    save_plt(plts, "mutation_analysis-mutated_sf3b1-by_cancer_type-violin", '.pdf', figs_dir, width=13, height=8)
    save_plt(plts, "mutation_analysis-mutated_u2af1-overall-violin", '.pdf', figs_dir, width=3, height=8)
    save_plt(plts, "mutation_analysis-mutated_u2af1-by_cancer_type-violin", '.pdf', figs_dir, width=13, height=8)

    save_plt(plts, "mutation_analysis-mutated_driver_type-overall-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "mutation_analysis-mutated_driver_type-by_cancer_type-violin", '.pdf', figs_dir, width=13, height=9)
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
        make_option("--diff_activity_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--gene_annotation_file", type="character"),
        make_option("--figs_dir", type="character"),
        make_option("--random_seed", type="integer", default=1234)
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    diff_activity_file = args[["diff_activity_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    gene_annotation_file = args[["gene_annotation_file"]]
    figs_dir = args[["figs_dir"]]
    random_seed = args[["random_seed"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_activity = read_tsv(diff_activity_file)
    splicing_factors = read_tsv(splicing_factors_file)
    gene_annotation = read_tsv(gene_annotation_file) %>%
        dplyr::rename(
            GENE = `Approved symbol`,
            ENSEMBL = `Ensembl gene ID`
        )
    program_activity_diff = read_tsv(program_activity_diff_file)
    metadata = read_tsv(metadata_file)
    mutations = read_tsv(mutations_file)
    driver_types = read_tsv(driver_types_file)
    
    # prep
    only_new_network = splicing_factors %>% filter(!in_viper_alone & in_new_w_viper) %>% pull(ENSEMBL)
    diff_activity = diff_activity %>%
        group_by(cancer_type) %>% # correct p-values for each type of cancer
        mutate(
            padj = p.adjust(pvalue, method="fdr"),
            is_significant = padj < THRESH_FDR,
            only_new_network = regulator %in% only_new_network
        ) %>%
        ungroup() %>%
        dplyr::rename("ENSEMBL"="regulator") %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = "ENSEMBL"
        )
    
    driver_activity = diff_activity %>%
        mutate(
            driver_type = ifelse(
                `condition_a-median`>0, "Oncogenic", "Tumor suppressor"
            )
        ) %>%
        filter(is_significant)
    
    # survival analysis
    program_activity_diff = program_activity_diff %>%
        left_join(
            metadata %>% distinct(sample, OS.time, OS, PFI.time, PFI, DFI.time, DFI), 
            by=c("index"="sample")
        )
    
    surv_metrics = c("OS","PFI","DFI")
    survival_analysis = lapply(surv_metrics, function(surv_metric){
        surv_time_col = sprintf("%s.time", surv_metric)
        surv_event_col = surv_metric
        vars_formula = "activity_diff"
        cancer_types = program_activity_diff %>% pull(cancer_type) %>% unique()
        result = lapply(cancer_types, function(cancer_type_oi){
            
            df = program_activity_diff %>% 
                filter(cancer_type==cancer_type_oi) %>% 
                distinct(.data[[surv_time_col]], .data[[surv_event_col]], activity_diff) %>%
                drop_na()
            
            fit_coxph = survival::coxph(
                as.formula(sprintf("Surv(%s, %s) ~ %s", surv_time_col, surv_event_col, vars_formula)), 
                data=df
            )

            result_surv = data.frame(
                n_obs = df %>% nrow(),
                surv_metric = surv_metric,
                cancer_type = cancer_type_oi,
                coxph_coef = fit_coxph[["coefficients"]][["activity_diff"]],
                coxph_pvalue = summary(fit_coxph)[["waldtest"]][["pvalue"]]
            )

            return(result_surv)
        }) %>% bind_rows()
        
        return(result)
    }) %>% bind_rows()
    
    # mutation analysis
    genes_oi = c("SF3B1","U2AF1")
    mutation_analysis = mutations %>%
        filter(gene %in% genes_oi) %>%
        left_join(program_activity_diff, by=c("sample"="index"))
    
    program_activity_diff = program_activity_diff %>%
        # overall
        left_join(
            mutations %>% 
                filter((effect!="Silent") & (gene%in%splicing_factors[["GENE"]])) %>%
                distinct(sample) %>%
                mutate(sf_mutated = TRUE),
            by = c("index"="sample")
        ) %>%
        mutate(sf_mutated = replace_na(sf_mutated, FALSE)) %>%
        # SF3B1
        left_join(
            mutations %>% 
                filter((effect!="Silent") & (gene%in%c("SF3B1"))) %>%
                distinct(sample) %>%
                mutate(mutated_sf3b1 = TRUE),
            by = c("index"="sample")
        ) %>%
        mutate(mutated_sf3b1 = replace_na(mutated_sf3b1, FALSE)) %>%
        # U2AF1
        left_join(
            mutations %>% 
                filter((effect!="Silent") & (gene%in%c("U2AF1"))) %>%
                distinct(sample) %>%
                mutate(mutated_u2af1 = TRUE),
            by = c("index"="sample")
        ) %>%
        mutate(mutated_u2af1 = replace_na(mutated_u2af1, FALSE)) %>%
        # by cancer program
        ## oncogenic
        left_join(
            mutations %>% 
                filter((effect!="Silent") & (gene%in% (driver_types %>% filter(driver_type=="Oncogenic") %>% pull(GENE)))) %>%
                distinct(sample) %>%
                mutate(mutated_oncogenic = TRUE),
            by = c("index"="sample")
        ) %>%
        mutate(mutated_oncogenic = replace_na(mutated_oncogenic, FALSE)) %>%
        ## tumor suppressor
        left_join(
            mutations %>% 
                filter((effect!="Silent") & (gene%in% (driver_types %>% filter(driver_type=="Tumor suppressor") %>% pull(GENE)))) %>%
                distinct(sample) %>%
                mutate(mutated_suppressor = TRUE),
            by = c("index"="sample")
        ) %>%
        mutate(mutated_suppressor = replace_na(mutated_suppressor, FALSE)) %>%
        ## classify
        mutate(
            mutated_driver_type = case_when(
                mutated_suppressor & mutated_oncogenic ~ "OncoTS",
                mutated_suppressor & !mutated_oncogenic ~ "TS only",
                !mutated_suppressor & mutated_oncogenic ~ "Onco only",
                !mutated_suppressor & !mutated_oncogenic & sf_mutated ~ "Non-driver-like",
                !mutated_suppressor & !mutated_oncogenic & !sf_mutated ~ "WT SFs"
            )
        )
    
    # plot
    plts = make_plots(driver_activity, program_activity_diff, survival_analysis, mutation_analysis, mutation_analysis)
    
    # make figdata
    figdata = make_figdata(driver_activity, program_activity_diff, survival_analysis, mutation_analysis, mutation_analysis)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
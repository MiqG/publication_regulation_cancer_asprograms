require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)

# variables
RANDOM_SEED = 1234

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "darkgreen"

PAL_DRIVER_TYPE = c(
    #"Non-driver"="lightgrey",
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
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

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_programs_in_differentiation")
# metadata_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets","metadata","CardosoMoreira2020.tsv.gz")
# genexpr_file = file.path(RAW_DIR,"viper_splicing_intermediate_files","datasets","genexpr_tpm","CardosoMoreira2020.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","CardosoMoreira2020-EX.tsv.gz")
# driver_types_file = file.path(ROOT,"results","new_empirical_network",'files','PANCAN','cancer_program.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","CardosoMoreira2020")

##### FUNCTIONS #####
plot_cancer_programs = function(protein_activity, genexpr){
    plts = list()
    
    X = protein_activity %>%
        drop_na(driver_type) %>%
        group_by(tissue, time, driver_type, study_accession, GENE) %>%
        summarize(activity = median(activity)) %>%
        ungroup()
    
#     for (tissue_oi in TISSUES){
#         plts[[sprintf("cancer_programs-differentiation_vs_activity-%s-violin",tissue_oi)]] = X %>%
#             filter(tissue==tissue_oi) %>%
#             mutate(time = cut(log10(time+1), breaks=10)) %>%
#             ggplot(aes(x=time, y=activity, group=interaction(time,driver_type))) +
#             geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
#             geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
#             fill_palette(PAL_DRIVER_TYPE) +
#             stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
#             geom_text(
#                 aes(y=-3, label=label, group=driver_type),
#                 . %>% count(time, driver_type) %>% mutate(label=paste0("n=",n)),
#                 size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
#             ) +
#             theme_pubr(x.text.angle=45) +
#             facet_wrap(~tissue, scales="free_x") +
#             theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#             labs(x="log10(Days Post Conception + 1)", y="Protein Activity", fill="Driver Type")
        
#         plts[[sprintf("cancer_programs-differentiation_vs_activity_diff-%s-line",study_oi)]] = X %>%
#             filter(tissue==tissue_oi) %>%
#             mutate(time = cut(log10(time+1), breaks=10)) %>%
#             mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
#             group_by(study_accession, time, driver_type) %>%
#             summarize(activity = median(activity)) %>%
#             ungroup() %>%
#             group_by(study_accession, time) %>%
#             summarize(activity_diff = sum(activity)) %>%
#             ungroup() %>%
#             ggline(
#                 x="time", y="activity_diff", color=PAL_DARK, 
#                 size=LINE_SIZE, point.size=0.05
#             ) +
#             theme_pubr(x.text.angle=45) +
#             facet_wrap(~tissue, scales="free") +
#             theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#             labs(x="log10(Days Post Conception + 1)", y="Protein Activity Diff.", fill="Driver Type")

#         plts[[sprintf("cancer_programs-differentiation_vs_mki67-%s-box",study_oi)]] = genexpr %>%
#             filter(tissue==tissue_oi & ID=="ENSG00000148773") %>%
#             mutate(time = cut(log10(time+1), breaks=10)) %>%
#             ggplot(aes(x=time, y=genexpr_tpm)) +
#             geom_boxplot(color=PAL_DARK, width=0.5, fill=NA, outlier.shape=NA) +
#             geom_point(color=PAL_DARK, position=position_jitter(0.1), size=1) +
#             geom_text(
#                 aes(y=-0.2, label=label),
#                 . %>% count(time) %>% mutate(label=paste0("n=",n)),
#                 size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
#             ) +
#             theme_pubr(x.text.angle=0) +
#             facet_wrap(~tissue, scales="free") +
#             theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#             labs(x="log10(Days Post Conception + 1)", y="MKI67 log2(TPM+1)", fill="Driver Type")
#     }
    
    plts[["cancer_programs-tissue_vs_differentiation-violin"]] = genexpr %>%
        distinct(time, tissue, sampleID) %>%
        mutate(time = log10(time+1)) %>%
        mutate(tissue = factor(tissue, levels=TISSUES)) %>%
        ggviolin(x="tissue", y="time", fill="tissue", color=NA, trim=TRUE) +
        geom_text(
                aes(y=1, label=label),
                . %>% count(tissue) %>% mutate(label=paste0("n=",n)),
                size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
            ) +
        fill_palette(get_palette("Dark2", length(TISSUES))) +
        guides(fill="none") +
        labs(x="Tissue", y="log10(Days Post Conception + 1)")
    
    
    plts[["cancer_programs-differentiation_vs_n_samples-bar"]] = genexpr %>%
        distinct(time, tissue, sampleID) %>%
        group_by(tissue) %>%
        mutate(time = as.numeric(cut(log10(time+1), breaks=10))) %>%
        ungroup() %>%
        count(tissue, time) %>%
        mutate(tissue = factor(tissue, levels=TISSUES)) %>%
        ggbarplot(x="time", y="n", fill="tissue", color=NA) +
        fill_palette(get_palette("Dark2", length(TISSUES))) +
        labs(x="log10(Days Post Conception + 1) Binned", y="Protein Activity Diff.", fill="Tissue")
    
    
    plts[["cancer_programs-differentiation_vs_activity_diff-line"]] = X %>%
        mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
        group_by(tissue) %>%
        mutate(time = as.numeric(cut(log10(time+1), breaks=10))) %>%
        ungroup() %>%
        group_by(tissue, study_accession, time, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(tissue, study_accession, time) %>%
        summarize(activity_diff = sum(activity)) %>%
        ungroup() %>%
        mutate(tissue = factor(tissue, levels=TISSUES)) %>%
        ggline(
            x="time", y="activity_diff", color="tissue", numeric.x.axis=TRUE,
            size=LINE_SIZE, linetype="solid", point.size=0.05
        ) +
        color_palette(get_palette("Dark2", length(TISSUES))) +
        geom_hline(yintercept=0, color="black", linetype="dashed", linewidth=LINE_SIZE) +
        stat_cor(aes(color=tissue), method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY, label.y = -1.5) +
        labs(x="log10(Days Post Conception + 1) Binned", y="Protein Activity Diff.", color="Tissue")
    
    
    
    return(plts)
}

make_plots = function(protein_activity, genexpr){
    plts = list(
        plot_cancer_programs(protein_activity, genexpr)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity, genexpr){
    figdata = list(
        "differentiation" = list(
            "genexpr" = genexpr,
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
    
    # activity
    save_plt(plts, "cancer_programs-tissue_vs_differentiation-violin",'.pdf', figs_dir, width=4, height=5)
    save_plt(plts, "cancer_programs-differentiation_vs_n_samples-bar",'.pdf', figs_dir, width=8, height=4.5)
    save_plt(plts, "cancer_programs-differentiation_vs_activity_diff-line",'.pdf', figs_dir, width=8, height=6)
    
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
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    genexpr = read_tsv(genexpr_file)
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    driver_types = read_tsv(driver_types_file)
    gc()
    
    # prep
    ## presence of potential cancer driver neojunctions in development vs gtex
    ## consider as detected junctions with at least 10 read counts
    metadata = metadata %>% 
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
    
    protein_activity = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        mutate(time = time_norm) %>%
        drop_na(tissue, time, activity) %>%
        mutate(
            condition_lab = sprintf(
                "%s | %s | %s", time, tissue, study_accession
            )
        ) %>%
        
        # summarize replicates
        group_by(condition_lab, time, tissue, study_accession, regulator) %>%
        summarize(
            activity = median(activity, na.rm=TRUE),
            abs_activity = abs(activity),
        ) %>%
        ungroup() %>%
        
        # add activity
        group_by(condition_lab) %>%
        arrange(activity) %>%
        mutate(
            activity_ranking = row_number()
        ) %>%
        arrange(abs_activity) %>%
        mutate(
            abs_activity_ranking = row_number()
        ) %>%
        ungroup() %>%
        left_join(driver_types, by=c("regulator"="ENSEMBL"))
    
    genexpr = genexpr %>%
        filter(ID%in%driver_types[["ENSEMBL"]] | ID=="ENSG00000148773") %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr_tpm") %>%
        left_join(metadata, by="sampleID") %>%
        mutate(time = time_norm) %>%
        drop_na(tissue, time, genexpr_tpm) %>%
        mutate(
            condition_lab = sprintf(
                "%s | %s | %s", sampleID, time, tissue, study_accession
            )
        ) %>%
        left_join(driver_types, by=c("ID"="ENSEMBL"))
    
    # plot
    plts = make_plots(protein_activity, genexpr)

    # make figdata
    figdata = make_figdata(protein_activity, genexpr)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}

##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

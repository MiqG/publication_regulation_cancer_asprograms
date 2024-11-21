#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - Interesting SFs
#    - https://www.nature.com/articles/s41588-021-00851-w :
#        - QKI (opposite brain vs heart)
#    - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6839889/ :
#        - Brain:
#            - PTBP1
#            - PTBP2
#            - SRRM4
#            - RBFOX1
#            - RBFOX3
#            - NOVA2
#            - KHDRBS3
#            - CTCF (brain epigenetics)
#            - TDP43 (neurological disorders)
#        - Muscle:
#            - CELF1
#            - RBFOX1
#            - RBFOX2
#            - RBM24
#            - MBNL1
#            - RBM20
#            - SF3B1
#            - PTBP1
#            - QKI
#        - Pancreas:
#            - NOVA1
#            - RBM4
#            - SRSF3
#            - SRSF10
#            - SLU7
#            - ESRP2
#        - Differentiation
#            - PTBP1 (smooth muscle cells)
#            - HNRNPA
#            - HNRNPB
#            - SNRP70
#            - HNRPLL
#            - MBNL2


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

GENES_OI = c(
    "PTBP1",
    "PTBP2",
    "SRRM4",
    "RBFOX1",
    "RBFOX3",
    "NOVA2",
    "KHDRBS3",
    "CTCF",
    "TDP43",
    "CELF1",
    "RBFOX1",
    "RBFOX2",
    "RBM24",
    "MBNL1",
    "RBM20",
    "SF3B1",
    "PTBP1",
    "QKI",
    "NOVA1",
    "RBM4",
    "SRSF3",
    "SRSF10",
    "SLU7",
    "ESRP2",
    "PTBP1",
    "HNRNPA",
    "HNRNPB",
    "SNRP70",
    "HNRPLL",
    "MBNL2"
)


DEV_STAGES = list(
    "PRJEB1195" = c(
        "iPSC","DEFINITIVE_ENDODERM","PRIMITIVE_GUT_TUBE","POSTERIOR_FOREGUT",
        "PANCREATIC_ENDODERM","LATE_PANCREATIC_ENDODERM","ENDOCRINE_CELLS","ADULT_BETA_CELL"
    ),
    "PRJNA379280"=c(
        "iPSC","ALVEOLAR_EPITHELIAL_PROGENITOR","ALVEOLAR_EPITHELIAL_TYPE_II","PRIMARY_ALVEOLAR_EPITHELIAL_TYPE_II"
    ),
    "PRJNA596331"=c(
        "iPSC","ACC_DORSAL",
        "NPC","ROSETTE","NEURONS_ALONE","NEURONS_PLUS_ASTROS"
    ),
    "PRJNA665705"=c(
        "iPSC","MESODERM","CARDIOMYOCYTE"
    )
)

STUDY_IPSC = data.frame(
    study_accession = c("PRJNA665705","PRJNA596331","PRJNA379280","PRJEB1195"),
    ipsc_type = c("Heart","Brain","Lung","Pancreas")
)


# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_programs_in_differentiation")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","ipsc_differentiation.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","ipsc_differentiation-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","ipsc_differentiation.tsv.gz")
# driver_types_file = file.path(ROOT,"results","new_empirical_network",'files','PANCAN','cancer_program.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","validation_ipsc_differentiation")

##### FUNCTIONS #####
plot_cancer_programs = function(protein_activity, genexpr){
    plts = list()
    
    X = protein_activity %>%
        drop_na(driver_type) %>%
        group_by(condition, driver_type, study_accession, GENE) %>%
        summarize(activity = median(activity)) %>%
        ungroup()
    
    for (study_oi in names(DEV_STAGES)){
        plts[[sprintf("cancer_programs-differentiation_vs_activity-%s-violin",study_oi)]] = X %>%
            filter(condition!="iPSC" & study_accession==study_oi & condition%in%DEV_STAGES[[study_oi]]) %>%
            mutate(condition=factor(condition, levels=DEV_STAGES[[study_oi]])) %>%
            ggplot(aes(x=condition, y=activity, group=interaction(condition,driver_type))) +
            geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
            geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
            fill_palette(PAL_DRIVER_TYPE) +
            stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
            geom_text(
                aes(y=-3, label=label, group=driver_type),
                . %>% count(condition, driver_type) %>% mutate(label=paste0("n=",n)),
                size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
            ) +
            theme_pubr(x.text.angle=45) +
            facet_wrap(~study_accession, scales="free") +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(x="IPSC Differentiation Stage", y="Protein Activity", fill="Driver Type")
        
        plts[[sprintf("cancer_programs-differentiation_vs_activity_diff-%s-line",study_oi)]] = X %>%
            filter(condition!="iPSC" & study_accession==study_oi & condition%in%DEV_STAGES[[study_oi]]) %>%
            mutate(activity = ifelse(driver_type=="Tumor suppressor", -activity, activity)) %>%
            group_by(study_accession, condition, driver_type) %>%
            summarize(activity = median(activity)) %>%
            ungroup() %>%
            group_by(study_accession, condition) %>%
            summarize(activity_diff = sum(activity)) %>%
            ungroup() %>%
            mutate(condition=factor(condition, levels=DEV_STAGES[[study_oi]])) %>%
            ggline(
                x="condition", y="activity_diff", color=PAL_DARK, 
                size=LINE_SIZE, point.size=0.05
            ) +
            theme_pubr(x.text.angle=45) +
            facet_wrap(~study_accession, scales="free") +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(x="IPSC Differentiation Stage", y="Protein Activity Diff.", fill="Driver Type")

        plts[[sprintf("cancer_programs-differentiation_vs_mki67-%s-box",study_oi)]] = genexpr %>%
            filter(ID=="ENSG00000148773" & study_accession==study_oi) %>%
            mutate(condition=factor(condition, levels=DEV_STAGES[[study_oi]])) %>%
            ggplot(aes(x=condition, y=genexpr_tpm)) +
            geom_boxplot(color=PAL_DARK, width=0.5, fill=NA, outlier.shape=NA) +
            geom_point(color=PAL_DARK, position=position_jitter(0.1), size=1) +
            geom_text(
                aes(y=-0.2, label=label),
                . %>% count(condition) %>% mutate(label=paste0("n=",n)),
                size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
            ) +
            theme_pubr(x.text.angle=0) +
            facet_wrap(~study_accession, scales="free") +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(x="IPSC Differentiation Stage", y="MKI67 log2(TPM+1)", fill="Driver Type")
    }
    
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
        "validation_drug_target_activity" = list(
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
    save_plt(plts, "cancer_programs-differentiation_vs_activity-PRJEB1195-violin", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_activity-PRJNA379280-violin", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_activity-PRJNA596331-violin", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_activity-PRJNA665705-violin", '.pdf', figs_dir, width=5, height=8)
    
    # activity diff
    save_plt(plts, "cancer_programs-differentiation_vs_activity_diff-PRJEB1195-line", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_activity_diff-PRJNA379280-line", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_activity_diff-PRJNA596331-line", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_activity_diff-PRJNA665705-line", '.pdf', figs_dir, width=5, height=8)
    
    # MKI67
    save_plt(plts, "cancer_programs-differentiation_vs_mki67-PRJEB1195-box", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_mki67-PRJNA379280-box", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_mki67-PRJNA596331-box", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "cancer_programs-differentiation_vs_mki67-PRJNA665705-box", '.pdf', figs_dir, width=5, height=8)
    
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
    protein_activity = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, activity) %>%
        mutate(
            condition_lab = sprintf(
                "%s | %s | %s", condition, cell_line_name, study_accession
            )
        ) %>%
        
        # summarize replicates
        group_by(condition_lab, condition, cell_line_name, study_accession, regulator) %>%
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
        drop_na(condition, genexpr_tpm) %>%
        mutate(
            condition_lab = sprintf(
                "%s | %s | %s", sampleID, condition, cell_line_name, study_accession
            )
        ) %>%
        
        # summarize replicates
        group_by(condition_lab, condition, cell_line_name, study_accession, ID) %>%
        summarize(
            genexpr_tpm = median(genexpr_tpm, na.rm=TRUE),
        ) %>%
        ungroup() %>%
        left_join(driver_types, by=c("ID"="ENSEMBL"))
    
    # plot
    plts = make_plots(protein_activity, genexpr)
    
    # make figdata
    #figdata = make_figdata(protein_activity, genexpr)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
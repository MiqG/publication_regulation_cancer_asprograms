require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","activity_estimation_w_genexpr")
# losses_file = file.path(RESULTS_DIR,"files","model_sf_activity","losses-merged.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","activity_modeling")

##### FUNCTIONS #####
plot_losses = function(losses){
    plts = list()
    
    X = losses %>%
        separate(loss_type, sep="_", into=c("subset","metric"), remove=FALSE)
    
    plts[["losses-epoch_vs_loss-line"]] = X %>%
        filter(metric=="loss") %>%
        ggscatter(x="epoch", y="loss", color="omic_regulon", size=0.1, alpha=0.5, palette="Paired") +
        geom_smooth(aes(fill=omic_regulon), color="black", linewidth=LINE_SIZE, linetype="dashed") +
        facet_wrap(~model_type + subset, nrow=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Epoch", y="SmoothL1Loss")
    
    plts[["losses-epoch_vs_pearson-line"]] = X %>%
        filter(metric=="pearson") %>%
        ggscatter(x="epoch", y="loss", color="omic_regulon", size=0.1, alpha=0.5, palette="Paired") +
        geom_smooth(aes(fill=omic_regulon), color="black", linewidth=LINE_SIZE, linetype="dashed") +
        facet_wrap(~model_type + subset, nrow=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Epoch", y="Pearson Corr. Coef.")
    
    plts[["losses-last_epoch_vs_pearson_vs_omic_regulon-box"]] = X %>%
        filter(metric=="pearson") %>%
        group_by(subset, omic_regulon, model_type) %>%
        slice_max(epoch) %>%
        ungroup() %>%
        ggplot(aes(x=omic_regulon, y=loss, group=interaction(omic_regulon, subset))) +
        geom_boxplot(aes(color=subset), fill=NA, outlier.shape=NA, position=position_dodge(0.9)) +
        geom_jitter(aes(color=subset), size=0.25, fill=NA, position=position_jitterdodge(jitter.width=0.5, dodge.width=0.9)) +
        color_palette("Paired") +
        theme_pubr() +
        facet_wrap(~model_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Network & Adjustment", y="Pearson Corr. Coef.")
    
    return(plts)
}

make_plots = function(losses){
    plts = list(
        plot_losses(losses)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(losses){
    figdata = list(
        "activity_modeling" = list(
            "losses" = losses
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
    save_plt(plts, "losses-epoch_vs_loss-line", '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, "losses-epoch_vs_pearson-line", '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, "losses-last_epoch_vs_pearson_vs_omic_regulon-box", '.pdf', figs_dir, width=7, height=5)
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
        make_option("--losses_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    losses_file = args[["losses_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    losses = read_tsv(losses_file)
    
    # plot
    plts = make_plots(losses)
    
    # make figdata
    figdata = make_figdata(losses)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
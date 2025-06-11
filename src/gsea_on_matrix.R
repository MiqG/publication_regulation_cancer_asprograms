# Note
# ----
# WE ARE RUNNING METAVIPER!

Sys.setenv(VROOM_CONNECTION_SIZE='100000000')
require(optparse)
require(tidyverse)
require(clusterProfiler)
require(BiocParallel)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# signature_file = file.path(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'log2_fold_change_tpm.tsv.gz')
# signature_file = "~/projects/publication_regulation_cancer_asprograms/results/network_inference/files/signatures/Hodis2022-invitro_eng_melanoc-genexpr.tsv.gz"
# msigdb_dir = file.path(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# gene_info_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# random_seed = 1234
# ontology_oi = "hallmarks"
# n_jobs = 4

##### FUNCTIONS #####
load_ontologies = function(msigdb_dir, cosmic_genes_file){
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "hallmarks_nomyc" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")), # ugly workaround
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "GO_CC" = read.gmt(file.path(msigdb_dir,"c5.go.cc.v7.4.symbols.gmt"))
    )
    return(ontologies)
}

gsea_for_sample = function(sample_name, signature, ontology) {
    cols_oi = c("Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "rank")
    
    # Create a named vector of fold changes for the sample
    genes = signature[[sample_name]]
    names(genes) = rownames(signature)
    genes = sort(genes, decreasing=TRUE)
    
    # run depending on the percentage of ties
    ratio_ties = table(genes)
    ratio_ties = ratio_ties / sum(ratio_ties)
    max_ratio_ties = ratio_ties[which.max(ratio_ties)]
    
    if (max_ratio_ties > 0.5){
        str = sprintf("Skipping %s because more than 0.5 tie ratio...", sample_name)
        print(str)
        
        # placeholder
        result = data.frame(matrix(NA, ncol = length(cols_oi), nrow = 1))
        colnames(result) = cols_oi
        
    } else {
        
        # Perform GSEA using clusterProfiler
        result = GSEA(geneList=genes, TERM2GENE=ontology, pvalueCutoff=1)

        # Convert to a data frame and add the sample name
        result = result@result %>%
                    as.tibble() %>% 
                    select(all_of(cols_oi))
    }
    
    result = result %>% mutate(sampleID = sample_name)

    return(result)
}

parseargs = function(){
    
    option_list = list( 
        make_option("--signature_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--ontology_oi", type="character"),
        make_option("--gene_info_file", type="character", default=NULL),
        make_option("--output_file", type="character"),
        make_option("--random_seed", type="integer", default=1234),
        make_option("--n_jobs", type="integer", default=1)
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    signature_file = args[["signature_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    ontology_oi = args[["ontology_oi"]]
    gene_info_file = args[["gene_info_file"]]
    output_file = args[["output_file"]]
    random_seed = args[["random_seed"]]
    n_jobs = args[["n_jobs"]]
    
    set.seed(args[["random_seed"]])
    
    # load
    signature = read_tsv(signature_file)
    ontologies = load_ontologies(msigdb_dir)
    
    # prep
    ## signature
    signature = signature %>% as.data.frame()
    rownames(signature) = signature[,1]
    signature = signature[,2:ncol(signature)]
    signature = signature %>% 
        dplyr::select(where(is.numeric))
    
    ## ontologies
    ontology = ontologies[[ontology_oi]]
    if (ontology_oi == "hallmarks_nomyc"){
        # ugly workaround
        genes_todrop = ontology %>% filter(str_detect(term,"_MYC_")) %>% pull(gene)
        ontology = ontology %>%
            filter(!(gene %in% genes_todrop))
    }
    
    ## translate ontology
    if (!is.null(gene_info_file)){
        gene_info = read_tsv(gene_info_file)
        ontology = ontology %>%
            left_join(
                gene_info,
                by=c("gene"="Approved symbol")
            ) %>%
            distinct(term, `Ensembl gene ID`) %>%
            drop_na() %>%
            rename(gene = `Ensembl gene ID`)
    }
    
    # run GSEA
    ## set up parallel processing
    register(MulticoreParam(workers=n_jobs))

    ## run
    result = bplapply(
        colnames(signature), gsea_for_sample, signature, ontology
    ) %>% bind_rows()

    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

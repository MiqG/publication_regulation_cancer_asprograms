Sys.setenv(VROOM_CONNECTION_SIZE = 500000)
require(optparse)
require(tidyverse)
require(viper)

##### FUNCTIONS #####
as_regulon_network = function(regulons){

    regulators = regulons[['regulator']] %>% unique()
    regulons = sapply(regulators, function(regulator_oi){
            X = regulons %>%
                filter(regulator %in% regulator_oi)
            X = list(
                tfmode = setNames(X[['tfmode']], X[['target']]),
                likelihood = X[['likelihood']]
            )
            return(X)
        }, simplify=FALSE)
    
    return(regulons)
}


load_networks = function(network_path, n_tails="two", patt=NULL){
    if (file.exists(network_path) && !dir.exists(network_path)){
        # network_path is a file, we load only that network (we'll run regular VIPER)
        network_files = list(network_path)
    }else if (dir.exists(network_path)){
        # network_path is a directory, we load all networks contained (we'll run metaVIPER)
        network_files = list.files(network_path, pattern=patt, full.names=TRUE)
    }else {
        stop("Invalid network_path.")
    }
    
    networks = sapply(network_files, function(network_file){
        print(network_file)
        network = read_tsv(network_file)
        if (nrow(network)>1 & n_tails=="one"){
            network = network %>% mutate(tfmode=abs(tfmode))
        }
        network = network %>%
            mutate(
                network_file = basename(network_file),
                study_accession = gsub("-","_",gsub(".tsv.gz","",network_file)),
                PERT_ID = sprintf("%s__%s", study_accession, ENSEMBL)
            )
        return(network)
    }, simplify=FALSE) %>% bind_rows()
    
    return(networks)
}


prep_regulons = function(networks){
    
    # from dataframe to regulon lists
    network_files = networks[["network_file"]] %>% unique()
    regulons = sapply(network_files, function(network_file){
        network = networks %>% filter(network_file==network_file)
        network = as_regulon_network(network)
        return(network)
    }, simplify=FALSE)
    
    # drop regulons that cannot be used
    regulons = sapply(regulons, function(regulon){
        to_keep = sapply(regulon, function(x){ length(x[[1]]) }) >= 25 # viper's default minsize
        regulon = regulon[to_keep]
        return(regulon)
    }, simplify=FALSE)
    
    # drop networks that cannot be used
    to_keep = sapply(regulons, length)>1
    regulons = regulons[to_keep]
    
    return(regulons)
}

run_viper = function(signature, regulons, shadow_correction="no"){
    # runs VIPER or metaVIPER depending on whether there are multiple regulons
    # in `regulons`
    
    pleiotropy = (shadow_correction=="yes") # TRUE/FALSE
    protein_activities = viper(signature, regulons, verbose=FALSE, pleiotropy=pleiotropy)
    
    return(protein_activities)
}

parseargs = function(){
    
    option_list = list( 
        make_option("--signature_file", type="character"),
        make_option("--regulons_path", type="character"),
        make_option("--output_file", type="character"),
        make_option("--random_seed", type="integer", default=1234),
        make_option("--shadow_correction", type="character", default="no"),
        make_option("--n_tails", type="character", default="two")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    signature_file = args[["signature_file"]]
    regulons_path = args[["regulons_path"]]
    random_seed = args[["random_seed"]]
    shadow_correction = args[["shadow_correction"]]
    n_tails = args[["n_tails"]]
    output_file = args[["output_file"]]
    
    set.seed(args[["random_seed"]])
    
    # load
    signature = read_tsv(signature_file)
    networks = load_networks(regulons_path, n_tails)
    
    # prep
    ## signature
    signature = signature %>% as.data.frame()
    rownames(signature) = signature[,1]
    signature = signature[,2:ncol(signature)]
    signature = signature %>% 
        dplyr::select(where(is.numeric))
    
    if (n_tails=="one"){
        signature = abs(signature)
    }
    
    # run regular viper
    regulons = prep_regulons(networks)
    result = run_viper(signature, regulons, shadow_correction)
    result = result %>% as.data.frame() %>% rownames_to_column('regulator')
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
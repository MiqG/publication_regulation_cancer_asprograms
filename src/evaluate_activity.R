require(optparse)
require(tidyverse)
require(viper)
require(pROC)
require(clusterProfiler)

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# signature_file="~/projects/publication_regulation_cancer_asprograms/data/raw/viper_splicing_intermediary_files/benchmark/signatures_tpm_encorekd_hepg2.tsv.gz"
# regulons_path="~/projects/publication_regulation_cancer_asprograms/results/network_inference/files/experimentally_derived_regulons_pruned-scgenexpr"
# eval_labels_file="~/projects/publication_regulation_cancer_asprograms/results/network_inference/files/regulon_evaluation_labels/ENCOREKD_HepG2.tsv.gz"
# shadow_correction = "no"
# n_tails = "two"
# method_activity = "viper"
# precomputed_activity_file = NULL

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
                study_accession = gsub("-","_",gsub(".tsv.gz","",network_file))
            )
        gc()
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

compute_correlations = function(signature, networks, cor_method, batch_size=NULL){
    
    avail_regulators = unique(networks[["regulator"]])
    
    if (!is.null(batch_size)){
        regulators_batches = split(avail_regulators, ceiling(seq_along(avail_regulators) / batch_size))
    }else{
        regulators_batches = list(avail_regulators)
    }
    
    protein_activities = sapply(names(signature), function(sample_oi){
            gc()
        
            x = signature[, sample_oi, drop=FALSE] %>%
                rownames_to_column("target") %>%
                drop_na()
            
            activities = lapply(regulators_batches, function(regulators_oi){
                gc()
                
                correlation = networks %>%
                    filter(regulator%in%regulators_oi) %>%
                    left_join(x, by="target") %>%
                    group_by(regulator) %>%
                    summarize(
                        activity = cor(
                            tfmode*likelihood, get(sample_oi), 
                            method=cor_method, use="pairwise.complete.obs"
                        )
                    ) %>%
                    ungroup()
                    
                gc()
                return(correlation)
            }) %>% 
            bind_rows() %>%
            deframe()
            
            gc()
            return(activities)
        }, simplify=FALSE) %>% 
        do.call(cbind, .)

    
    return(protein_activities)
}


compute_enrichment = function(signature, networks_eval){
    # adapt networks
    ontology = networks_eval %>%
        distinct(regulator, target) %>%
        rename(term=regulator, gene=target)
    
    # compute enrichment for each signature
    protein_activities = sapply(names(signature), function(sample_oi){
            x = signature[, sample_oi, drop=FALSE] %>%
                rownames_to_column("target") %>%
                drop_na() %>%
                deframe() %>%
                sort(decreasing=TRUE)
            x = x[x!=0] # remove 0s to avoid errors
            
            # are detected splicing changes in network?
            in_network = sum(names(x) %in% ontology[["gene"]]) > 1
            if (in_network){
                result = GSEA(geneList=x, TERM2GENE=ontology, pvalueCutoff=1.1)

                activities = result@result %>%
                    as.tibble() %>% 
                    distinct(Description, NES) %>%
                    deframe()

                s = sprintf("sample:%s and data length:%s", sample_oi, length(activities))
                print(s)
            } else {
                print("No exon found in network, making activities placeholder...")
                
                regulators = unique(ontology[["term"]])
                activities = setNames(rep(NA, length(regulators)), regulators)
            }
        
            gc()
            return(activities)
        }, simplify=FALSE) %>%
        do.call(cbind, .)
        
    return(protein_activities)
}

compute_activity = function(signature, networks_eval, method_activity, shadow_correction="no", batch_size=NULL){
    
    if (!is.null(batch_size)){
        sample_batches = split(names(signature), ceiling(seq_along(names(signature)) / batch_size))
    }else{
        sample_batches = list(names(signature))
    }
    
    protein_activities = lapply(sample_batches, function(batch){
        
        if (method_activity=="viper"){

            # run metaVIPER with all possible regulons
            regulons = prep_regulons(networks_eval)
            pleiotropy = (shadow_correction=="yes")
            protein_activities = viper(signature[, batch, drop=FALSE], regulons, verbose=FALSE, pleiotropy=pleiotropy)

        } else if (method_activity=="correlation_pearson"){

            # correlate networks with each sample
            protein_activities = compute_correlations(signature[, batch, drop=FALSE], networks_eval, "pearson")

        } else if (method_activity=="correlation_spearman"){

            # correlate networks with each sample
            protein_activities = compute_correlations(signature[, batch, drop=FALSE], networks_eval, "spearman")

        } else if (method_activity=="gsea"){

            protein_activities = compute_enrichment(signature[, batch, drop=FALSE], networks_eval)

        }
    
    })
    protein_activities = do.call(cbind, protein_activities)
    
    return(protein_activities)
}


compute_average_precision = function(precision, recall) {
    
    precision[is.na(precision)] = 0
    recall[is.na(recall)] = 0

    if(length(precision) != length(recall)) {
    stop("Precision and recall vectors must be of the same length")
    }

    # Apply trapezoidal rule
    ap = sum((recall[-1] - recall[-length(recall)]) * precision[-1])
  
  return(ap)
}


compute_auc = function(x, y){
    x[is.na(x)] = 0
    y[is.na(y)] = 0
    
    x = sort(x)
    y = sort(y)
    
    diffs.x = x[-1] - x[-length(x)]
    means.vert = (y[-1] + y[-length(y)])/2
    auc = sum(means.vert * diffs.x)

    return(auc)
}

compute_curves = function(x, grouping_var){
    # precision recall curves
    ## for each perturbation experiment (PERT_ID), across splicing factors
    ## we only consider those experiments in the dataset that perturbed a splicing
    ## factor whose activity was measured
    
    x = x %>%
        group_by(get(grouping_var)) %>%

        # drop those with exactly the same activity
        filter(n_distinct(activity) > 1) %>%
        ungroup %>%
        
        group_by(get(grouping_var)) %>%
        # we changed the sign of overexpression activities, so we expect
        # perturbed splicing factors to have low activity; to evaluate 
        # we normalize such that the lower the closer to one and
        # the higher the closer to 0
        mutate(
            activity_norm = -activity,
            activity_norm = (activity_norm - min(activity_norm)) / (max(activity_norm) - min(activity_norm))
        ) %>%
        ungroup()
    
    curves_by_group = x %>% 
        group_by(get(grouping_var)) %>%
        do({
            # summary
            pert_sf_clean = .$pert_sf[is.finite(.$activity_norm)] # we replaced NAs
            
            # rank percentile
            mean_rp = mean(.$activity_norm[.$pert_sf=="Regulator"], na.rm=TRUE)
            median_rp = median(.$activity_norm[.$pert_sf=="Regulator"], na.rm=TRUE)
            
            # ROC curves
            roc_curve = roc(
                response=.$pert_sf, predictor=.$activity_norm, 
                direction="<", levels=c("NotRegulator", "Regulator")
            )
            coords_curve = coords(roc_curve, transpose = FALSE, ret = "all")
            
            # Add PERT_ID to the results by getting it from the grouped data
            coords_curve %>%
                mutate(precision = ifelse(is.na(precision) & threshold==Inf, 1, precision)) %>%
                mutate(
                    auc_roc = compute_auc(fpr, sensitivity),
                    auc_pr = compute_auc(recall, precision),
                    avg_pr = compute_average_precision(precision, recall),
                    mean_rank_percentile = mean_rp,
                    median_rank_percentile = median_rp,
                    grouping_var = grouping_var,
                    n_pos_class = sum(pert_sf_clean=="Regulator"),
                    n_total = length(pert_sf_clean),
                    curves_type = "by_group"
                )
        }) %>%
        ungroup()

    curves_combined = x %>% 
        do({
            # summary
            pert_sf_clean = .$pert_sf[is.finite(.$activity_norm)]

            # rank percentile
            mean_rp = mean(.$activity_norm[.$pert_sf=="Regulator"], na.rm=TRUE)
            median_rp = median(.$activity_norm[.$pert_sf=="Regulator"], na.rm=TRUE)
            
            # ROC curves
            roc_curve = roc(
                response=.$pert_sf, predictor=.$activity_norm, 
                direction="<", levels=c("NotRegulator", "Regulator")
            )
            coords_curve = coords(roc_curve, transpose = FALSE, ret = "all")
            
            # Add PERT_ID to the results by getting it from the grouped data
            coords_curve %>%
                mutate(precision = ifelse(is.na(precision) & threshold==Inf, 1, precision)) %>%
                mutate(
                    auc_roc = compute_auc(fpr, sensitivity),
                    auc_pr = compute_auc(recall, precision),
                    avg_pr = compute_average_precision(precision, recall),
                    mean_rank_percentile = mean_rp,
                    median_rank_percentile = median_rp,
                    grouping_var = "all",
                    n_pos_class = sum(pert_sf_clean=="Regulator"),
                    n_total = length(pert_sf_clean),
                    curves_type = "combined"
                )
        })
    
    curves = bind_rows(curves_by_group, curves_combined)
    
    return(curves)
}


evaluate_activities = function(protein_activities, eval_labels){
    # prep
    X = protein_activities %>%
        as.data.frame() %>%
        rownames_to_column("regulator") %>%
        pivot_longer(-regulator, names_to="PERT_ID", values_to="activity") %>%
        left_join(eval_labels, by="PERT_ID") %>%
        mutate(
            # change sign of protein activities from OVEREXPRESSION signatures
            activity = ifelse(PERT_TYPE=="OVEREXPRESSION",-activity, activity),
            # label perturbed regulator (splicing factor)
            pert_sf = ifelse(regulator==PERT_ENSEMBL, "Regulator", "NotRegulator")
        ) %>%
        drop_na()
    
    # we consider only perturbations where the regulator was measured, keeping all regulators
    # but regulators have different activities
    perts_oi = X %>%
        filter(PERT_ENSEMBL==regulator) %>%
        distinct(PERT_ID) %>%
        pull()
    curves_real_pert = X %>% 
        filter(PERT_ID %in% perts_oi) %>%
        compute_curves(., grouping_var="PERT_ID") %>%
        rename(
            PERT_ID = `get(grouping_var)`
        )  %>%
        mutate(
            eval_direction = "perturbations",
            eval_type = "real"
        )
    curves_random_pert = X %>% 
        filter(PERT_ID %in% perts_oi) %>%
        group_by(PERT_ID) %>%
        mutate(activity = sample(activity)) %>%
        ungroup() %>%
        compute_curves(., grouping_var="PERT_ID") %>%
        rename(
            PERT_ID = `get(grouping_var)`
        )  %>%
        mutate(
            eval_direction = "perturbations",
            eval_type = "random"
        )
    
    # we consider only regulators that were perturbed, keeping all perturbations
    sfs_oi = X %>%
        filter(PERT_ENSEMBL == regulator) %>%
        distinct(regulator) %>%
        pull()
    curves_real_sf = X %>% 
        filter(regulator %in% sfs_oi) %>%
        compute_curves(., grouping_var="regulator") %>%
        rename(
            regulator = `get(grouping_var)`
        )  %>%
        mutate(
            eval_direction = "regulators",
            eval_type = "real"
        )
    curves_random_sf = X %>% 
        filter(regulator %in% sfs_oi) %>%
        group_by(regulator) %>%
        mutate(activity = sample(activity)) %>%
        ungroup() %>%
        compute_curves(., grouping_var="regulator") %>%
        rename(
            regulator = `get(grouping_var)`
        )  %>%
        mutate(
            eval_direction = "regulators",
            eval_type = "random"
        )
    
    # merge
    evaluation = curves_random_pert %>%
        bind_rows(curves_real_pert) %>%
        bind_rows(curves_random_sf) %>%
        bind_rows(curves_real_sf)
    
    return(evaluation)
}

parseargs = function(){
    
    option_list = list( 
        make_option("--signature_file", type="character", default=NULL),
        make_option("--regulons_path", type="character", default=NULL),
        make_option("--eval_labels_file", type="character", default=NULL),
        make_option("--precomputed_activity_file", type="character", default=NULL),
        make_option("--output_file", type="character"),
        make_option("--random_seed", type="integer", default=1234),
        make_option("--shadow_correction", type="character", default="no"),
        make_option("--method_activity", type="character", default="no"),
        make_option("--n_tails", type="character", default="two")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    signature_file = args[["signature_file"]]
    regulons_path = args[["regulons_path"]]
    eval_labels_file = args[["eval_labels_file"]]
    precomputed_activity_file = args[["precomputed_activity_file"]]
    random_seed = args[["random_seed"]]
    shadow_correction = args[["shadow_correction"]]
    n_tails = args[["n_tails"]]
    method_activity = args[["method_activity"]]
    output_file = args[["output_file"]]
    
    set.seed(args[["random_seed"]])
    
    # load
    eval_labels = read_tsv(eval_labels_file)
    if (is.null(precomputed_activity_file)){
    
        # we must compute protein activity here
        signature = read_tsv(signature_file)
        networks = load_networks(regulons_path, n_tails)
    
    } else {
        
        precomputed_activity = read_tsv(precomputed_activity_file)
    
    } 
    
    # prep
    if (is.null(precomputed_activity_file)){

        ## signature
        signature = signature %>% as.data.frame()
        rownames(signature) = signature[,1]
        signature = signature[,2:ncol(signature)]
        signature = signature %>% 
            dplyr::select(where(is.numeric))

        if (n_tails=="one"){
            signature = abs(signature)
        }
        
        # filter out networks found in the evaluation labels
        eval_dataset = gsub(".tsv.gz","",basename(eval_labels_file))
        networks_eval = networks %>% filter(!str_detect(study_accession, eval_dataset))
        networks_kept = networks_eval %>% distinct(study_accession) %>% pull(study_accession)
        sprintf("Will evaluate %s on networks from %s .", signature_file, paste0(networks_kept, collapse=", "))

        networks_per_regulator = networks_eval %>% 
            distinct(regulator, network_file) %>%
            count(regulator, name="n_networks_per_regulator")

        # compute activity
        sprintf("Computing activities with %s ...", method_activity)
        protein_activity = compute_activity(signature, networks_eval, method_activity, shadow_correction)
        
    } else {
        
        # no need to compute protein activity
        protein_activity = precomputed_activity %>% as.data.frame()
        
        # make compatible
        rownames(protein_activity) = protein_activity[,1]
        protein_activity = protein_activity[,2:ncol(protein_activity)]
        protein_activity = protein_activity %>% dplyr::select(where(is.numeric))
        
    }
    
    # run viper and evaluate predicted protein activities
    print("Evaluating...")
    result = evaluate_activities(protein_activity, eval_labels)
    
    if (is.null(precomputed_activity_file)){
        
        result = result %>%
                left_join(networks_per_regulator, by="regulator") %>%
                mutate(
                    regulon_set_id = basename(regulons_path),
                    signature_id = basename(eval_labels_file) %>% gsub(".tsv.gz","",.),
                    shadow_correction = shadow_correction,
                    n_tails = n_tails,
                    method_activity = method_activity
                )
        
    } else {
        
        result = result %>%
            mutate(
                regulon_set_id = "precomputed",
                signature_id = "precomputed",
                shadow_correction = NA,
                n_tails = NA,
                method_activity = "precomputed"
            )
        
    }
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
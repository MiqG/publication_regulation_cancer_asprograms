import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
NETWORKS_DIR = os.path.join(ROOT,"results","network_inference")
RESULTS_DIR = os.path.join(ROOT,"results","carcinogenic_switch_regulation")

REGULON_DIR = os.path.join(ROOT,"results","activity_estimation_w_genexpr")
REGULON_DIRS = {
    "bulkgenexpr": os.path.join(REGULON_DIR,"files","experimentally_derived_regulons_pruned-bulkgenexpr"),
    "scgenexpr": os.path.join(REGULON_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr"),
    "bulkscgenexpr": os.path.join(REGULON_DIR,"files","experimentally_derived_regulons_pruned-bulkscgenexpr")
}

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

PERT_GENEXPR_FILES = {
    "ReplogleWeissman2022_rpe1": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    #"ReplogleWeissman2022_K562_essential": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz")
}

MODEL_TYPES = ["fclayer"]
OMIC_REGULONS = ["bulkgenexpr"]
GENE_SETS = ["pert_splicing_factors","pert_splicing_factors_random"]
K_CROSS_VALIDATION = 5

##### RULES #####
rule all:
    input:
        # estimate splicing factor activity from gene expression
        ## raw
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"), dataset=PERT_GENEXPR_FILES.keys(), omic_regulon=OMIC_REGULONS),
        ## adjusted
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}-adjusted_{model_type}.tsv.gz"), dataset=PERT_GENEXPR_FILES.keys(), model_type=MODEL_TYPES, omic_regulon=OMIC_REGULONS),
        
        # shortest paths among splicing factors
        os.path.join(SUPPORT_DIR,'pert_splicing_factors_random-symbol.txt'),
        expand(os.path.join(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_{set_oi}.tsv.gz'), set_oi=GENE_SETS),
        
        # make figures
        os.path.join(RESULTS_DIR,"figures","upstream_regulators"),
        

        
rule compute_protein_activity:
    input:
        signature = lambda wildcards: PERT_GENEXPR_FILES[wildcards.dataset],
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_regulon]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
       """

rule adjust_activity:
    input:
        activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"),
        weights = [os.path.join(REGULON_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth").format(k=k, omic_regulon="{omic_regulon}", model_type="{model_type}") for k in range(K_CROSS_VALIDATION)],
        input_regulators = os.path.join(REGULON_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","input_regulators.tsv.gz"),
        output_regulators = os.path.join(REGULON_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","output_regulators.tsv.gz")
    params:
        weights = ",".join([os.path.join(REGULON_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth").format(k=k, omic_regulon="{omic_regulon}", model_type="{model_type}") for k in range(K_CROSS_VALIDATION)]),
        model_type = "{model_type}",
        script_dir = SRC_DIR
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}-adjusted_{model_type}.tsv.gz")
    shell:
        """
        python {params.script_dir}/vipersp/scripts/adjust_genexpr_sf_activity.py \
                    --activity_file={input.activity} \
                    --weights_files="{params.weights}" \
                    --input_regulators_file={input.input_regulators} \
                    --output_regulators_file={input.output_regulators} \
                    --model_type={params.model_type} \
                    --output_file={output}
        """

rule random_gene_set:
    input:
        gt_gene_set = os.path.join(SUPPORT_DIR,"pert_splicing_factors-symbol.txt"),
        avail_genes = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        random_gene_set = os.path.join(SUPPORT_DIR,'pert_splicing_factors_random-symbol.txt')
    params:
        random_seed = 1234
    run:
        import pandas as pd
        
        # load
        gt_gene_set = pd.read_table(input.gt_gene_set, header=None)[0].to_list()
        avail_genes = pd.read_table(input.avail_genes)
        random_seed = params.random_seed
        
        # get random genes
        n_genes = len(gt_gene_set)
        subset = avail_genes.sample(n=n_genes, replace=False, random_state=random_seed)
        
        # save
        subset[["Approved symbol"]].drop_duplicates().to_csv(output.random_gene_set, sep="\t", index=False, header=False)
        
        print("Done!")
    
rule shortest_paths_stringdb:
    input:
        ppi = os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz'),
        sources = os.path.join(SUPPORT_DIR,"{set_oi}-symbol.txt"),
        targets = os.path.join(SUPPORT_DIR,"pert_splicing_factors-symbol.txt")
    output:
        os.path.join(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_{set_oi}.tsv.gz')
    threads: 10
    shell:
        """
        nice python scripts/ppi_path_lengths.py \
                    --ppi_file={input.ppi} \
                    --sources_file={input.sources} \
                    --targets_file={input.targets} \
                    --output_file={output} \
                    --n_jobs={threads}
        """ 
        
rule figures_upstream_regulators:
    input:
        protein_activity_rpe1 = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_rpe1-bulkgenexpr-adjusted_fclayer.tsv.gz"),
        metadata_rpe1 = os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-conditions.tsv.gz"),
        gene_info = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz"),
        cancer_program = os.path.join(ROOT,"results","new_empirical_network",'files','PANCAN','cancer_program.tsv.gz'),
        msigdb_dir = os.path.join(RAW_DIR,"MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs"),
        cosmic_genes = os.path.join(RAW_DIR,"COSMIC","cancer_gene_census.tsv"),
        regulons_dir = os.path.join(ROOT,"results","new_empirical_network","files","experimentally_derived_regulons_pruned_w_viper_networks-EX"),
        event_info = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv'),
        splicing_factors = os.path.join(SUPPORT_DIR,"supplementary_tables","splicing_factors.tsv"),
        shortest_paths_pert_sfs = os.path.join(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_pert_splicing_factors.tsv.gz'),
        shortest_paths_random = os.path.join(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_pert_splicing_factors_random.tsv.gz')
    output:
        directory(os.path.join(RESULTS_DIR,"figures","upstream_regulators"))
    shell:
        """
        Rscript scripts/figures_upstream_regulators.R \
                    --protein_activity_rpe1_file={input.protein_activity_rpe1} \
                    --metadata_rpe1_file={input.metadata_rpe1} \
                    --gene_info_file={input.gene_info} \
                    --cancer_program_file={input.cancer_program} \
                    --msigdb_dir={input.msigdb_dir} \
                    --cosmic_genes_file={input.cosmic_genes} \
                    --regulons_dir={input.regulons_dir} \
                    --event_info_file={input.event_info} \
                    --splicing_factors_file={input.splicing_factors} \
                    --shortest_paths_pert_sfs_file={input.shortest_paths_pert_sfs} \
                    --shortest_paths_random_file={input.shortest_paths_random} \
                    --figs_dir={output}
        """

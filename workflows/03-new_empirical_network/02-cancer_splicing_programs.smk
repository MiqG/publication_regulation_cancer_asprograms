import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
SRC_DIR = os.path.join(ROOT,"src")
RESULTS_DIR = os.path.join(ROOT,"results","new_empirical_network")
TCGA_DIR = os.path.join(RAW_DIR,"viper_splicing_intermediate_files","tcga")
REGULONS_PATH = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned_w_viper_networks-EX")

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}
PADJ_METHOD = 'fdr_bh'
OMICS = ["protein_activity"]

CANCER_TYPES_PTSTN = [
    'BLCA',
    'BRCA',
    'COAD',
    'HNSC',
    'KICH',
    'KIRC',
    'KIRP',
    'LIHC',
    'LUAD',
    'LUSC',
    'PRAD',
    'STAD',
    'THCA',
    'UCEC'
]

DIFF_CONDITIONS = {
    "PrimaryTumor_vs_SolidTissueNormal":{
        "a":"PrimaryTumor",
        "b":"SolidTissueNormal"
    }
}

DIFF_CANCER_TYPES = {
    "PrimaryTumor_vs_SolidTissueNormal": CANCER_TYPES_PTSTN
}

rule all:
    input:
        # calculate signatures
        expand(os.path.join(RESULTS_DIR,"files","signatures","{cancer}-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz"), cancer=CANCER_TYPES_PTSTN),

        # compute viper SF activities
        ## PT vs STN
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample}-EX.tsv.gz"), cancer=CANCER_TYPES_PTSTN, sample=["PrimaryTumor_vs_SolidTissueNormal"]),
        
        # differential analyses
        ## SF activities
        expand(os.path.join(RESULTS_DIR,'files',"diff_protein_activity",'{cancer}-{comparison}.tsv.gz'), comparison=["PrimaryTumor_vs_SolidTissueNormal"], cancer=CANCER_TYPES_PTSTN),
        expand(os.path.join(RESULTS_DIR,'files','PANCAN','{omic}-mannwhitneyu-{comparison}.tsv.gz'), comparison=DIFF_CANCER_TYPES.keys(), omic=OMICS),
        
        ## define cancer program
        os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz'),
        
        ## compute program activity differences across cancers
        os.path.join(RESULTS_DIR,"files","protein_activity",'PANCAN-PrimaryTumor_vs_SolidTissueNormal-program_activity_diff.tsv.gz'),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","cancer_splicing_program")
        
        
rule compute_signature_pt_vs_stn:
    input:
        splicing_pt = os.path.join(PREP_DIR,"event_psi","{cancer}-PrimaryTumor-EX.tsv.gz"),
        splicing_stn = os.path.join(PREP_DIR,"event_psi","{cancer}-SolidTissueNormal-EX.tsv.gz")
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{cancer}-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz")
    run:
        import pandas as pd
        
        splicing_pt = pd.read_table(input.splicing_pt, index_col=0)
        splicing_stn = pd.read_table(input.splicing_stn, index_col=0)
        
        # subtract median from the other group
        signature = pd.concat([
            splicing_pt - splicing_stn.median(axis=1).values.reshape(-1,1), 
            splicing_stn - splicing_pt.median(axis=1).values.reshape(-1,1)
        ], axis=1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{cancer}-{sample}-EX.tsv.gz"),
        regulons_path = REGULONS_PATH
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample}-EX.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
        
rule compute_differential_protein_activity:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{comparison}-EX.tsv.gz"),
        metadata = os.path.join(TCGA_DIR,'metadata','{cancer}.tsv.gz')
    output:
        os.path.join(RESULTS_DIR,'files',"diff_protein_activity",'{cancer}-{comparison}.tsv.gz')
    params:
        script_dir = SRC_DIR,
        padj_method = PADJ_METHOD,
        condition_a = lambda wildcards: DIFF_CONDITIONS[wildcards.comparison]["a"],
        condition_b = lambda wildcards: DIFF_CONDITIONS[wildcards.comparison]["b"],
        comparison_col = "sample_type_clean",
        sample_col = "sampleID"
    shell:
        """
        python {params.script_dir}/MannWhitneyU.py \
                    --data_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --sample_col={params.sample_col} \
                    --comparison_col={params.comparison_col} \
                    --condition_a={params.condition_a} \
                    --condition_b={params.condition_b} \
                    --output_file={output} \
                    --padj_method={params.padj_method} 
        """
        
        
rule combine_differential_results:
    input:
        diff_files = lambda wildcards: [os.path.join(RESULTS_DIR,'files',"diff_{omic}",'{cancer}-{comparison}.tsv.gz').format(cancer=cancer, comparison=wildcards.comparison, omic="{omic}") for cancer in DIFF_CANCER_TYPES[wildcards.comparison]],
    params:
        comparison = '{comparison}'
    output:
        os.path.join(RESULTS_DIR,'files','PANCAN','{omic}-mannwhitneyu-{comparison}.tsv.gz')
    run:
        import os
        import pandas as pd
        
        dfs = []
        for diff_file in input.diff_files:
            # combine
            df = pd.read_table(diff_file)
            
            # add cancer type
            cancer_type = os.path.basename(diff_file).replace(".tsv.gz","").split("-")[0]
            print(cancer_type)
            df['cancer_type'] = cancer_type
            
            dfs.append(df)
            
            del df
            
        dfs = pd.concat(dfs)  
        dfs.to_csv(output[0], sep='\t', index=False, compression='gzip')
        
        print("Done!")
        
        
rule define_cancer_program:
    input:
        diff_activity = os.path.join(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz'),
        gene_annotation = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        cancer_program = os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
    params:
        thresh_fdr = 0.05,
        thresh_n_sum = 5
    run:
        import pandas as pd
        import numpy as np
        
        # load data
        diff_activity = pd.read_table(input.diff_activity)
        gene_annotation = pd.read_table(input.gene_annotation)
        
        # oncogenic and tumor suppressor splicing factors are those that preferentially
        # activate or inactivate recurrently across cancer types
        
        ## keep significant
        diff_activity = diff_activity.loc[diff_activity["padj"] < params.thresh_fdr].copy()
        
        ## add driver type
        diff_activity["driver_type"] = diff_activity["condition_a-median"].apply(
            lambda x: "Oncogenic" if x>0 else "Tumor suppressor"
        )
        
        ## count recurrence
        cancer_program = diff_activity.groupby(
            ["regulator","driver_type"]
        ).size().reset_index().rename(columns={0:"n"})
        cancer_program["n_sign"] = [
            n if driver_type=="Oncogenic" else -n 
            for n, driver_type in cancer_program[["n","driver_type"]].values
        ]
        n_sum = cancer_program.groupby("regulator")["n_sign"].sum().reset_index().rename(
            columns={"n_sign":"n_sum"}
        )
        cancer_program = pd.merge(cancer_program, n_sum, on="regulator", how="left")
        
        ## classify
        cancer_program = cancer_program.loc[
            np.abs(cancer_program["n_sum"]) > params.thresh_n_sum
        ].copy()
        cancer_program = cancer_program.loc[
            cancer_program.groupby("regulator")["n"].idxmax()
        ].copy()
        
        # add gene annotations s
        gene_annotation = gene_annotation.rename(
            columns={"Approved symbol":"GENE", "Ensembl gene ID":"ENSEMBL"}
        )
        cancer_program = cancer_program.rename(
            columns={"regulator":"ENSEMBL"}
        )
        cancer_program = pd.merge(
            cancer_program, gene_annotation[["GENE","ENSEMBL"]], on="ENSEMBL", how="left"
        )
        
        # save
        cancer_program.to_csv(output.cancer_program, **SAVE_PARAMS)
        
        print("Done!")
        

rule compute_program_activity_differences:
    input:
        activity = [os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample}-EX.tsv.gz").format(cancer=canc, sample=samp) for canc in CANCER_TYPES_PTSTN for samp in ["PrimaryTumor_vs_SolidTissueNormal"]],
        cancer_program = os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
    output:
        program_activity_diff = os.path.join(RESULTS_DIR,"files","protein_activity",'PANCAN-PrimaryTumor_vs_SolidTissueNormal-program_activity_diff.tsv.gz')
    run:
        import os
        import pandas as pd
        
        activity_files = input.activity
        cancer_program = pd.read_table(input.cancer_program)
        sf_onco = cancer_program.loc[cancer_program["driver_type"]=="Oncogenic","ENSEMBL"]
        sf_ts = cancer_program.loc[cancer_program["driver_type"]=="Tumor suppressor","ENSEMBL"]
        
        program_activity_diff = []
        for f in activity_files:
            activity = pd.read_table(f, index_col=0)
            
            cancer_type = os.path.basename(f).split("-")[0]
            comparison = os.path.basename(f).split("-")[1]
            feature = os.path.basename(f).split("-")[2].replace(".tsv.gz","")
            
            # compute median activity across each driver type for each sample
            median_onco = activity.loc[activity.index.isin(sf_onco)].median(axis=0) 
            median_ts = activity.loc[activity.index.isin(sf_ts)].median(axis=0)
            diff = (median_onco - median_ts).reset_index(name="activity_diff")
            diff["cancer_type"] = cancer_type
            diff["comparison"] = comparison
            diff["feature"] = feature
            
            program_activity_diff.append(diff)
        program_activity_diff = pd.concat(program_activity_diff)
            
        program_activity_diff.to_csv(output.program_activity_diff, **SAVE_PARAMS)
        
        print("Done!")
        
rule figures_cancer_splicing_program:
    input:
        diff_activity = os.path.join(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz'),
        splicing_factors = os.path.join(RESULTS_DIR,"figures","empirical_psi_networks","figdata","eda","splicing_factors.tsv.gz"),
        gene_annotation = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz"),
        program_activity_diff = os.path.join(RESULTS_DIR,"files","protein_activity",'PANCAN-PrimaryTumor_vs_SolidTissueNormal-program_activity_diff.tsv.gz'),
        metadata = os.path.join(RAW_DIR,'UCSCXena','TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.tsv'),
        mutations = os.path.join(RAW_DIR,'UCSCXena','TCGA','snv','mc3.v0.2.8.PUBLIC.xena.gz'),
        driver_types = os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz'),
    output:
        directory(os.path.join(RESULTS_DIR,"figures","cancer_splicing_program"))
    shell:
        """
        Rscript scripts/figures_cancer_splicing_program.R \
                    --diff_activity_file={input.diff_activity} \
                    --splicing_factors_file={input.splicing_factors} \
                    --gene_annotation_file={input.gene_annotation} \
                    --program_activity_diff_file={input.program_activity_diff} \
                    --metadata_file={input.metadata} \
                    --mutations_file={input.mutations} \
                    --driver_types_file={input.driver_types} \
                    --figs_dir={output}
        """
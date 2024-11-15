import os
import pandas as pd
import numpy as np

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
ARTICLES_DIR = os.path.join(RAW_DIR,'articles')
ARTICLE_DIR = os.path.join(ARTICLES_DIR,'Rogalska2024')
VASTDB_DIR = os.path.join(RAW_DIR,'VastDB')

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

# load metadata
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-PRJEB49033-Rogalska2024.tsv'))
metadata = metadata.loc[metadata["library_source"]=="TRANSCRIPTOMIC"]

## URLS to download
URLS = metadata['fastq_ftp'].str.split(';').str[0].apply(os.path.dirname).to_list()
URLS = {os.path.basename(url): url for url in URLS}
SAMPLES = list(URLS.keys())
N_SAMPLES = len(SAMPLES)

## fastq sizes
SIZES = metadata.set_index("run_accession")["fastq_bytes"].str.split(";").apply(lambda x: max(np.array(x, dtype=int))).to_dict()
SIZE_THRESH = 5e9

ENDS = ["1","2"]

##### RULES #####
rule all:
    input:
        # download metadata
        os.path.join(ARTICLE_DIR,"metadata.tsv"),
        
        # download .fastq
        expand(os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}'), end=ENDS, sample=SAMPLES),
        
        # quantify event PSI and gene expression TPM
        expand(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES),
        
        # combine
        os.path.join(ARTICLE_DIR,'vast_out','.done/vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES),
        
        # tidy PSI
        os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        '.done/Rogalska2024.done'


rule download_metadata:
    params:
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB49033&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created,bam_ftp,bam_bytes,bam_md5&format=tsv&download=true&limit=0"
    output:
        os.path.join(RAW_DIR,"articles","Rogalska2024","metadata.tsv")
    threads: 1
    resources:
        # runtime = int(3600*0.5), # 30 min in seconds
        runtime = int(60*0.5), # 30 min in minutes
        memory = 2
    shell:
        """
        set -eo pipefail
        
        wget --user-agent="Chrome" \
             --no-check-certificate \
             "{params.url}" \
             -O {output}
             
        echo "Done!"
        """
        
        
rule download:
    params:
        sample = '{sample}',
        end = "{end}",
        url = lambda wildcards: URLS[wildcards.sample],
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir="~/repositories/vast-tools/"
    output:
        download_done = os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}')
    threads: 1
    resources:
        # runtime = 24*3600, # 24h in seconds
        runtime = 24*60, # 24h in minutes
        memory = 10
    shell:
        """
        # download
        echo "Downloading {params.sample}..."
        
        wget --user-agent="Chrome" \
             --no-check-certificate \
             {params.url}/{params.sample}_{params.end}.fastq.gz \
             -O {params.fastqs_dir}/{params.sample}_{params.end}.fastq.gz
        
        touch {output.download_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule align:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir="~/repositories/vast-tools/",
        vast_out = directory(os.path.join(ARTICLE_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
        download_done = [os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS]
    output:
        align_done = touch(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'))
    threads: 12
    resources:
        runtime = lambda wildcards: 24*60 if SIZES[wildcards.sample]>SIZE_THRESH else 6*60, # 6h or 24h in minutes
        memory = 20
    shell:
        """
        # align paired reads
        echo "Aligning {params.sample}..."
        {params.bin_dir}/vast-tools align \
                    {params.fastqs_dir}/{params.sample}_1.fastq.gz \
                    {params.fastqs_dir}/{params.sample}_2.fastq.gz \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --expr \
                    --EEJ_counts \
                    --cores {threads} \
                    --output {params.vast_out}
        echo "Finished aligning {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule vasttools_combine:
    input:
        done = [os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}').format(sample=sample) for sample in SAMPLES],
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(ARTICLE_DIR,'vast_out','.done','vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES)),
        tpm = os.path.join(ARTICLE_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES),
        psi = os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES)
    params:
        bin_dir="~/repositories/vast-tools/",
        folder = os.path.join(ARTICLE_DIR,'vast_out')
    threads: 12
    resources:
        # runtime = 86400, # 24h in seconds
        runtime = 24*60, # 24h in minutes
        memory = 60
    shell:
        """
        # group results
        echo "Grouping results..."
        mkdir -p {params.folder}/to_combine
        ln -s {params.folder}/*/to_combine/* {params.folder}/to_combine/
        mkdir -p {params.folder}/expr_out
        ln -s {params.folder}/*/expr_out/* {params.folder}/expr_out/
        
        # combine runs
        echo "Combining runs..."
        {params.bin_dir}/vast-tools combine \
                    --cores {threads} \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --keep_raw_reads \
                    --keep_raw_incl \
                    --output {params.folder} \
                    --TPM \
                    -C
        
        # compress outputs
        echo "Compressing outputs..."
        gzip -f {params.folder}/raw_incl/*
        gzip -f {params.folder}/raw_reads/*
        gzip -f {params.folder}/*.tab
        
        # remove grouped results
        echo "Removing grouped results..."
        rm -r {params.folder}/to_combine
        rm -r {params.folder}/expr_out
        
        echo "Done!"
        """
    
    
rule vasttools_tidy:
    input:
        os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES) ## combined table
    output:
        touch('.done/Rogalska2024.done'),
        tidy = os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir="~/repositories/vast-tools/"
    threads: 1
    resources:
        # runtime = 3600*12, # 12h in seconds
        runtime = 60*12, # 12h in minutes
        memory = 60
    shell:
        """
        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        echo "Done!"
        """
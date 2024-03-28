"""
Download and align
------------------
- metadata KD RNAseq: https://www.encodeproject.org/metadata/?assay_title=shRNA+RNA-seq&files.file_type=fastq&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&searchTerm=shrna&type=Experiment
- metadata eCLIP: https://www.encodeproject.org/metadata/?files.biological_replicates_formatted=Rep+1%2C+Rep+2&files.output_type=peaks&files.file_type=bed+narrowPeak&type=Experiment&cart=%2Fcarts%2Faf6b3405-3695-4658-8748-eb393dde1ede%2F&files.assembly=GRCh38
"""
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
DATA_DIR = os.path.join(ROOT,'data')
SUPPORT_DIR = os.path.join(ROOT,'support')
ENCODE_DIR = os.path.join(DATA_DIR,'ENCODE') # parent
ENCORE_DIR = os.path.join(ENCODE_DIR,'ENCORE') # child
VASTDB_DIR = os.path.join(DATA_DIR,'VastDB')

INDEXES_DIR = os.path.join(DATA_DIR,"genome_indexes")
ADAPTERS_FILE = os.path.join(SUPPORT_DIR,"trimmomatic_adapters.fasta")
STAR_DIR = os.path.join("~/repositories","STAR-2.7.10a","bin","Linux_x86_64_static")
FASTP_DIR = os.path.join("~/repositories")
TF_DIR = os.path.join(DATA_DIR,"ThermoFisher")

ENDS = ["1","2"]

# load metadata (previously downloaded)
metadata = pd.read_table(os.path.join(ENCORE_DIR,'metadata','ENCORE.tsv'))
metadata["sample"] = metadata["dbxrefs"].str.replace("SRA:","")
metadata["end"] = metadata["Paired end"].astype('str')

# dictionary of paired end fastq files
SAMPLES_ENCORE = {}
for sample_oi in metadata["sample"].unique():
    SAMPLES_ENCORE[sample_oi] = metadata.loc[
        metadata["sample"]==sample_oi,
        ["end","File accession"]]\
        .set_index("end")\
        .to_dict()["File accession"]
    
URLS = metadata.set_index("File accession")["File download URL"].to_dict()
SIZES = metadata.set_index("File accession")["Size"].to_dict()
SIZE_THRESH = 5e9

# load metadata eCLIP
metadata_eclip = pd.read_table(os.path.join(ENCORE_DIR,'metadata','ENCORE-eCLIP.tsv'))
SAMPLES_ENCORE_ECLIP = metadata_eclip["File accession"]
URLS_ECLIP = metadata_eclip.set_index("File accession")["File download URL"].to_dict()

rule all:
    input:
        # RNA-seq
        ## download metadata
        os.path.join(ENCORE_DIR,'metadata','ENCORE.tsv'),
        
        ## Download .fastq files and Quantify splicing and expression
        expand(os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}'), end=ENDS, sample=SAMPLES_ENCORE.keys()),
        
        ## trim adapters
        expand(os.path.join(ENCORE_DIR,"fastqs_trimmed",".done","{sample}"), sample=list(SAMPLES_ENCORE.keys())),
        
        ## quantify mRNA levels
        ### load genome in memory
        # ".done/genome_loaded",
        ### align
        expand(os.path.join(ENCORE_DIR,"STAR",".done","{sample}_trimmed"), sample=list(SAMPLES_ENCORE.keys())),

        ### remove genome from memory
        ".done/genome_unloaded",
        
        ### merge read gene counts
        os.path.join(ENCORE_DIR,"STAR","merged_counts.tab.gz"),
        
        ## quantify splicing
        # expand(os.path.join(ENCORE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES_ENCORE.keys()),
        # expand(os.path.join(ENCORE_DIR,'fastqs','.done_rm','{sample}'), sample=SAMPLES_ENCORE.keys()),
        ## Combine into single tables
        #os.path.join(ENCORE_DIR,'vast_out','.done/vasttools_combine'),
        ## Tidy PSI
        #'.done/ENCORE.done',
        #os.path.join(ENCORE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
        
        # eCLIP
        ## download metadata
        os.path.join(ENCORE_DIR,'metadata','ENCORE-eCLIP.tsv'),
        
        ## download narrow peak bed files
        expand(os.path.join(ENCORE_DIR,"eclip_peaks","{sample}.bed.gz"), sample=SAMPLES_ENCORE_ECLIP),

        ## map peaks to exons (considering neighboring introns)
        expand(os.path.join(ENCORE_DIR,"eclip_peaks_mapped","{sample}.tsv.gz"), sample=SAMPLES_ENCORE_ECLIP),
        
        ## merge peaks
        os.path.join(ENCORE_DIR,"eclip_peaks_mapped","merged.tsv.gz")
        
        
rule download_metadata:
    params:
        metadata = "'https://www.encodeproject.org/metadata/?assay_title=shRNA+RNA-seq&files.file_type=fastq&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&searchTerm=shrna&type=Experiment'"
    output:
        metadata = os.path.join(ENCORE_DIR,'metadata','ENCORE.tsv'),
        readme = os.path.join(ENCORE_DIR,'metadata','README.md')
    shell:
        """
        # metadata
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.metadata} -O {output.metadata}
        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """
        
rule download:
    params:
        sample = '{sample}',
        end = '{end}',
        url = lambda wildcards: URLS[SAMPLES_ENCORE[wildcards.sample][wildcards.end]],
        fastqs_dir = os.path.join(ENCORE_DIR,'fastqs'),
    output:
        download_done = os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}')
    threads: 1
    resources:
        runtime = 5400,
        memory = 2
    group: "ENCORE"
    shell:
        """
        # download
        echo "Downloading {params.sample}..."
        
        wget --user-agent="Chrome" --no-check-certificate {params.url} -O {params.fastqs_dir}/{params.sample}_{params.end}.fastq.gz
        
        touch {output.download_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        

rule trim_sequencing_adapters:
    input:
        [os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}').format(sample="{sample}", end=end) for end in ENDS]
    output:
        [os.path.join(ENCORE_DIR,'fastqs_trimmed','{sample}_{end}.fastq.gz').format(sample="{sample}", end=end) for end in ENDS],
        os.path.join(ENCORE_DIR,"fastqs_trimmed","reports","html","{sample}.html"),
        os.path.join(ENCORE_DIR,"fastqs_trimmed","reports","json","{sample}.json"),
        trimming_done = touch(os.path.join(ENCORE_DIR,"fastqs_trimmed",".done","{sample}"))
    params:
        sample = "{sample}",
        fastqs_dir = os.path.join(ENCORE_DIR,"fastqs"),
        trimmer_dir = FASTP_DIR,
        adapters = ADAPTERS_FILE
    threads: 6
    resources:
        runtime = int(3600*1), # 1 h
        memory = 2 # GB
    shell:
        """
        # trim reads
        nice {params.trimmer_dir}/fastp.0.23.2 \
                --thread {threads} \
                --in1 {params.fastqs_dir}/{params.sample}_1.fastq.gz \
                --in2 {params.fastqs_dir}/{params.sample}_2.fastq.gz \
                --adapter_fasta {params.adapters} \
                --detect_adapter_for_pe \
                --out1 {params.fastqs_dir}_trimmed/{params.sample}_1.fastq.gz \
                --out2 {params.fastqs_dir}_trimmed/{params.sample}_2.fastq.gz \
                --html {params.fastqs_dir}_trimmed/reports/html/{params.sample}.html \
                --json {params.fastqs_dir}_trimmed/reports/json/{params.sample}.json
        
        echo "Done!"
        """
        
        
rule load_star_genome:
    input:
        index_dir = os.path.join(INDEXES_DIR,"STAR","hg38_ERCC92","99")
    output:
        touch(".done/genome_loaded")
    params:
        aligner_dir = STAR_DIR
    threads: 6
    resources:
        runtime = int(3600*1), # 1 h
        memory = 35 # GB
    shell:
        """
        set -euo pipefail

        nice {params.aligner_dir}/STAR \
                    --genomeDir {input.index_dir} \
                    --genomeLoad LoadAndExit \
                    --runThreadN {threads}
        echo "Done!"
        """
        
        
rule quantify_genexpr_star:
    input:
        #".done/genome_loaded",
        [os.path.join(ENCORE_DIR,'fastqs_trimmed','{sample}_{end}.fastq.gz').format(sample="{sample}", end=end) for end in ENDS],
        index_dir = os.path.join(INDEXES_DIR,"STAR","hg38_ERCC92","99")
    output:
        align_done = touch(os.path.join(ENCORE_DIR,"STAR",".done","{sample}_trimmed"))
    params:
        sample = "{sample}",
        aligner_dir = STAR_DIR,
        # tmp_dir = os.path.join("/solid","manglada","{sample}"),
        tmp_dir = os.path.join("/tmp","{sample}"),
        fastqs_dir = os.path.join(ENCORE_DIR,"fastqs_trimmed"),
        star_out = directory(os.path.join(ENCORE_DIR,"STAR","{sample}_trimmed"))
    threads: 6
    resources:
        runtime = lambda wildcards: int(3600*1) if any([SIZES[SAMPLES_ENCORE[wildcards.sample]["1"]]>SIZE_THRESH,SIZES[SAMPLES_ENCORE[wildcards.sample]["2"]]>SIZE_THRESH]) else int(3600*2.5), # most 1h is enough; some needed 2.5h (more reads).
        memory = 2 # GB
    shell:
        """
        set -euo pipefail
        
        # align reads
        mkdir -p {params.star_out}
        nice {params.aligner_dir}/STAR \
                    --genomeDir {input.index_dir} \
                    --genomeLoad LoadAndKeep \
                    --readFilesIn \
                    {params.fastqs_dir}/{params.sample}_1.fastq.gz \
                    {params.fastqs_dir}/{params.sample}_2.fastq.gz \
                    --readFilesCommand zcat \
                    --quantMode GeneCounts \
                    --outTmpDir {params.tmp_dir} \
                    --outFileNamePrefix "{params.star_out}/" \
                    --outSAMtype BAM Unsorted \
                    --runThreadN {threads}
        
        echo "Done!"
        """
        
        
rule remove_loaded_genome:
    input:
        # ".done/genome_loaded",
        [os.path.join(ENCORE_DIR,"STAR",".done","{sample}_trimmed").format(sample=sample) for sample in SAMPLES_ENCORE.keys()],
        index_dir = os.path.join(INDEXES_DIR,"STAR","hg38_ERCC92","99")
    output:
        touch(".done/genome_unloaded")
    params:
        aligner_dir = STAR_DIR
    resources:
        runtime = int(3600*0.25), # 1 h
        memory = 1 # GB
    shell:
        """
        set -euo pipefail

        nice {params.aligner_dir}/STAR \
                    --genomeDir {input.index_dir} \
                    --genomeLoad Remove
        echo "Done!"
        """        
        
        
rule combine_star:
    input:
        [os.path.join(ENCORE_DIR,"STAR",".done","{sample}_trimmed").format(sample=sample) for sample in SAMPLES_ENCORE.keys()]
    output:
        counts = os.path.join(ENCORE_DIR,"STAR","merged_counts.tab.gz")
    params:
        counts_col = 2,
        star_out = os.path.join(ENCORE_DIR,"STAR"),
        done_dir = os.path.join(ENCORE_DIR,"STAR",".done")
    shell:
        """
        set -euo pipefail

        SAMPLES=$(basename -a {params.done_dir}/*_trimmed)
        COUNTER=0
        for SAMPLE in $SAMPLES; do
            if [ "$COUNTER" -eq "0" ]; then
                # for first iteration, add gene column
                cat {params.star_out}/$SAMPLE/ReadsPerGene.out.tab | sed 1,4d | cut -f1,{params.counts_col} | sed "1i ENSEMBL\t$SAMPLE" > {params.star_out}/tmp
                
            else
                # take only mRNA counts for the rest
                paste {params.star_out}/tmprev <(cat {params.star_out}/$SAMPLE/ReadsPerGene.out.tab | sed 1,4d | cut -f{params.counts_col} | sed "1i $SAMPLE") > {params.star_out}/tmp
            fi
            
            mv {params.star_out}/tmp {params.star_out}/tmprev
            
            COUNTER=$[$COUNTER +1]
        done
        
        mv {params.star_out}/tmprev {output.counts}
        gzip -f {output.counts}
        mv {output.counts}.gz {output.counts}
        
        echo "Done!"
        """        
        
        
rule quantify_splicing:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(ENCORE_DIR,'fastqs'),
        bin_dir="~/repositories/vast-tools/",
        vast_out = directory(os.path.join(ENCORE_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
        download_done = [os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS]
    output:
        align_done = touch(os.path.join(ENCORE_DIR,'vast_out','.done','{sample}'))
    threads: 16
    resources:
        runtime = lambda wildcards: 86400 if any([SIZES[SAMPLES_ENCORE[wildcards.sample]["1"]]>SIZE_THRESH,SIZES[SAMPLES_ENCORE[wildcards.sample]["2"]]>SIZE_THRESH]) else 21600, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    group: "ENCORE"
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

        
rule delete_fastqs:
    input:
        download_done = [os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS],
        align_done = os.path.join(ENCORE_DIR,'vast_out','.done','{sample}')
    output:
        touch(os.path.join(ENCORE_DIR,'fastqs','.done_rm','{sample}'))
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(ENCORE_DIR,'fastqs')
    threads: 1
    resources:
        runtime = 300,
        memory = 1
    group: "ENCORE"
    shell:
        """
        rm {params.fastqs_dir}/{params.sample}*
        
        echo "Done!"
        """
        
        
rule vasttools_combine:
    input:
        [os.path.join(ENCORE_DIR,'vast_out','{sample}').format(sample=sample) for sample in SAMPLES_ENCORE.keys()],
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(ENCORE_DIR,'vast_out','.done','vasttools_combine')),
        tpm = os.path.join(ENCORE_DIR,'vast_out','TPM-hg38-1097.tab.gz'),
        psi = os.path.join(ENCORE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-1097.tab.gz')
    params:
        bin_dir="~/repositories/vast-tools/",
        folder = os.path.join(ENCORE_DIR,'vast_out')
    threads: 16
    resources:
        runtime = 172800, # 48h
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
                    --output {params.folder} \
                    --TPM
        
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
        os.path.join(ENCORE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-1097.tab.gz') ## combined table
    output:
        touch('.done/ENCORE.done'),
        tidy = os.path.join(ENCORE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir="~/repositories/vast-tools/"
    threads: 1
    resources:
        runtime = 43200, # 12h
        memory = 60
    shell:
        """
        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        echo "Done!"
        """

rule download_metadata_eclip:
    params:
        metadata = "'https://www.encodeproject.org/metadata/?files.biological_replicates_formatted=Rep+1%2C+Rep+2&files.output_type=peaks&files.file_type=bed+narrowPeak&type=Experiment&cart=%2Fcarts%2Faf6b3405-3695-4658-8748-eb393dde1ede%2F&files.assembly=GRCh38'"
    output:
        metadata = os.path.join(ENCORE_DIR,'metadata','ENCORE-eCLIP.tsv')
    shell:
        """
        # metadata
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.metadata} -O {output.metadata}
        
        echo Done!
        """
        
rule download_eclip:
    params:
        sample = '{sample}',
        url = lambda wildcards: URLS_ECLIP[wildcards.sample],
    output:
        peaks = os.path.join(ENCORE_DIR,'eclip_peaks','{sample}.bed.gz')
    threads: 1
    resources:
        runtime = 5400,
        memory = 2
    group: "ENCORE_eCLIP"
    shell:
        """
        # download
        echo "Downloading {params.sample}..."
        
        nice wget --user-agent="Chrome" --no-check-certificate {params.url} -O {output.peaks}
        
        echo "Finished downloading {params.sample}."
        echo "Done!"
        """
        
rule map_peaks_to_exons:
    input:
        peaks = os.path.join(ENCORE_DIR,'eclip_peaks','{sample}.bed.gz'),
        metadata = os.path.join(ENCORE_DIR,"metadata","ENCORE-eCLIP.tsv"),
        event_info = os.path.join(DATA_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
    output:
        os.path.join(ENCORE_DIR,"eclip_peaks_mapped","{sample}.tsv.gz")
    run:
        import pandas as pd
        import pyranges as pr
        
        # load
        ## read peaks file
        col_names = [
            "Chromosome",
            "Start",
            "End",
            "sample_id",
            "is_significant",
            "Strand",
            "log10_pvalue",
            "log2_fc",
            "unknown1",
            "unknown2"
        ]
        peaks = pd.read_table(
            input.peaks,
            header=None,
            names=col_names
        )
        peaks = peaks.reset_index()
        peaks = peaks.rename(columns={"index":"peak_id"})
        ## read peaks metadata
        metadata = pd.read_table(input.metadata)
        ## read vastdb events
        events = pd.read_table(input.event_info)
        
        # assign exons to peaks if a peak falls on the exon or one of its neighboring introns
        ## only event type of interest
        events = events.loc[events["EVENT"].str.contains("EX")].copy()
        ## process event coordinates
        events["Chromosome"] = (
            events["COORD_o"].str.split(":").str[0]
        )
        events["EVENT_start"] = (
            events["COORD_o"].str.split(":").str[1].str.split("-").str[0]
        ).astype("int")
        events["EVENT_end"] = (
            events["COORD_o"].str.split(":").str[1].str.split("-").str[1].astype("int")
        )
        events["Strand"] = (
            events["REF_CO"].str.split(":").str[2]
        )
        ## start and end will be the next splice site junctions
        events["Start"] = events["CO_C1"].str.split(":").str[1].str.split("-").str[1]
        events["End"] = events["CO_C2"].str.split(":").str[1].str.split("-").str[0]
        
        # intersect
        ## prepare
        X = pr.PyRanges(
            events[["EVENT", "EVENT_start", "EVENT_end", "Chromosome", "Start", "End"]],
            int64=True,
        )
        Y = pr.PyRanges(peaks, int64=True)
        ## join
        margin = 0
        mapping = X.join(Y, slack=margin, how="left", report_overlap=True)
               
        # prepare outputs
        mapping = mapping.df.loc[~mapping.df["is_significant"].isin([-1])]
        mapping = mapping.rename(columns={
            "Start_b":"peak_start",
            "End_b":"peak_end"
        })
        
        # save
        mapping.to_csv(output[0], sep="\t", compression="gzip", index=False)
        
        print("Done!")
        
        
rule merge_mapped:
    input:
        mapped_peaks = [os.path.join(ENCORE_DIR,"eclip_peaks_mapped","{sample}.tsv.gz").format(sample=sample) for sample in SAMPLES_ENCORE_ECLIP]
    output:
        os.path.join(ENCORE_DIR,"eclip_peaks_mapped","merged.tsv.gz")
    run:
        import pandas as pd
        
        df = pd.concat([
            pd.read_table(file) for file in input.mapped_peaks
        ])
        
        df.to_csv(output[0], sep="\t", compression="gzip", index=False)
        
        print("Done!")

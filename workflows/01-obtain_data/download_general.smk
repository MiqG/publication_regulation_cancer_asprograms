import os
import pandas as pd
import numpy as np

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### RULES #####
rule all:
    input:
        # HGNC - gene ids
        os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz"),
        
        # harmonizome
        os.path.join(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz"),
        
        # MSigDB
        os.path.join(RAW_DIR,"MSigDB",'msigdb_v7.4'),

        # STRINGDB
        os.path.join(RAW_DIR,"STRINGDB",'9606.protein.links.full.v11.5.txt.gz'),
        os.path.join(RAW_DIR,"STRINGDB",'9606.protein.aliases.v11.5.txt.gz'),
        
        # preprocessed files from Anglada-Girotto2024
        os.path.join(RAW_DIR,"viper_intermediate_files")
        
        
rule download_hgnc:
    params:
        annot = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_pub_ensembl_id&col=gd_aliases&col=gd_enz_ids&col=family.id&col=gd_app_name&col=gd_status&col=gd_locus_type&col=gd_locus_group&col=gd_prev_sym&col=gd_prev_name&col=gd_name_aliases&col=gd_pub_chrom_map&col=gd_date2app_or_res&col=gd_date_mod&col=gd_date_sym_change&col=gd_date_name_change&col=gd_pub_acc_ids&col=gd_pub_eg_id&col=gd_mgd_id&col=gd_other_ids&col=gd_other_ids_list&col=gd_pubmed_ids&col=gd_pub_refseq_ids&col=family.name&col=gd_ccds_ids&col=gd_vega_ids&col=gd_lsdb_links&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
    output:
        annot = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    shell:
        """
        wget --user-agent="Chrome" --no-check-certificate "{params.annot}" -O {output.annot}
        
        gzip --force {output.annot}
        
        mv {output.annot}.gz {output.annot}
        """
        
rule download_stringdb:
    params:
        links = 'https://stringdb-static.org/download/protein.links.full.v11.5/9606.protein.links.full.v11.5.txt.gz',
        aliases = 'https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz'
    output:
        links = os.path.join(RAW_DIR,"STRINGDB",'9606.protein.links.full.v11.5.txt.gz'),
        aliases = os.path.join(RAW_DIR,"STRINGDB",'9606.protein.aliases.v11.5.txt.gz'),
        readme = os.path.join(RAW_DIR,"STRINGDB",'README.md')
    shell:
        """
        # download
        wget --no-check-certificate {params.links} -O {output.links}
        wget --no-check-certificate {params.aliases} -O {output.aliases}
        
        # add readme
        echo "Downloaded on $(date)." > {output.readme}
        echo Done!
        """
        
rule download_harmonizome_gene_sets:
    params:
        chea = "https://maayanlab.cloud/static/hdfs/harmonizome/data/cheappi/gene_set_library_crisp.gmt.gz"
    output:
        chea = os.path.join(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz"),
        readme = os.path.join(RAW_DIR,"Harmonizome",'README.md')
    shell:
        """
        # aneuploidy
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.chea} -O {output.chea}        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """
        
rule download_msigdb:
    params:
        db = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/msigdb_v7.4_files_to_download_locally.zip'
    output:
        db = directory(os.path.join(RAW_DIR,"MSigDB",'msigdb_v7.4')),
        readme = os.path.join(RAW_DIR,"MSigDB",'README.md')
    shell:
        """
        # download
        wget --user-agent="Chrome" --no-check-certificate {params.db} -O {output.db}.zip
        
        # unzip
        unzip {output.db}.zip -d {output.db}
        
        echo "Downloaded on $(date)." > {output.readme}
        echo Done!
        """
        
rule download_viper_intermediate_files:
    params:
        files = "https://figshare.com/ndownloader/files/50598684"
    output:
        files = directory(os.path.join(RAW_DIR,'viper_splicing_intermediate_files'))
    shell:
        """
        # download
        wget --user-agent="Chrome" --no-check-certificate {params.files} -O {output.files}.zip
        
        # unzip
        unzip {output.files}.zip -d {output.files}
        
        echo Done!
        """
        
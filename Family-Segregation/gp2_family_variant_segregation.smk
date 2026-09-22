#!/usr/bin/env python
import glob
import os
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
from collections import defaultdict

shell.prefix("set -o pipefail; umask 002; ")  # set g+w
configfile: "configs/config.yaml"

#Parse wildcards for the main Snakemake
chr_list = list(range(1,23))+["X", "Y"]
chr_id= ["chr" + str(i) for i in chr_list]
ped_path = "meta/segregation_report_ped.tsv"
FAMS = pd.read_csv(ped_path, sep='\t', usecols=['family_id'])['family_id'].dropna().unique().tolist()

rule all:
    input:
        expand("families/reports/{fam}_segregation_variants.csv", fam=FAMS)

rule subset_plink_vcf:
    input:
        sample = "meta/samples_to_keep_family.txt",
        pgen = "input/pfiles/{chrom}.pgen"
    output:
        vcf = "vcfs/qced/{chrom}.vcf.gz"
    params:
       in_prefix = "input/pfiles/{chrom}",
       out_prefix = "vcfs/qced/{chrom}"
    threads:
        10
    resources:
        nodes = 1,
        mem_mb_per_cpu = 10000
    shell:
        r"""
        # convert plink pfiles to vcf 

        bins/plink2 --pfile {params.in_prefix} \
                     --keep {input.sample} \
                     --threads {threads} \
                     --mac 1 \
                     --export id-paste=iid bgz vcf-4.2 \
                     --output-chr chrM \
                     --out {params.out_prefix}
        """

rule slivar_clean:
    input:
        vcf = "vcfs/qced/{chrom}.vcf.gz"
    output:
        "vcfs/cleaned/{chrom}.vcf.gz"
    params:
        slivar_function = "script/slivar-functions.v0.2.8.js",
    container:
        "docker://zihhuafang/slivar_modified:0.2.8"
    threads:
        1
    resources:
        nodes = 1,
        mem_mb_per_cpu = 8000
    shell:
        r"""
        # annotate with gnomAD and TOPMed freq using slivar
        
        slivar expr --js /mnt/slivar/{params.slivar_function} \
                    -g /mnt/slivar/TOPMed_freeze8_PASS.zip \
                    -g /mnt/slivar/gnomADv4.0_wfilter.zip \
                    --info 'variant.ALT[0] != "*" && variant.call_rate >= 0.95' \
                    --vcf {input.vcf} \
        | bgzip -c > {output}

        tabix -p vcf --force {output}
        """

rule vep:
    input:
        vcf = rules.slivar_clean.output,
        CADD = "CADD/{chrom}.tsv.gz",
    output:
        vcf = "vcfs/vep/{chrom}.vcf.gz"
    container:
        "docker://zihhuafang/ensembl_vep_loftee:v111"
    threads: 4
    resources:
        nodes = 1,
        mem_mb_per_cpu = 8000
    params:
        ref_genome= config["docker_ref_genome"]
    shell:
        r"""
        vep \
        -i {input.vcf} \
        -o {output.vcf} \
        --compress_output bgzip \
        --vcf \
        --no_stats \
        --ASSEMBLY GRCh38 \
        --buffer_size 3000 \
        --fasta {params.ref_genome} \
        --sift b \
        --polyphen b \
        --ccds \
        --hgvs \
        --symbol \
        --numbers \
        --domains \
        --regulatory \
        --canonical \
        --protein \
        --biotype \
        --af_gnomade \
        --af_gnomadg \
        --max_af \
        --pubmed \
        --uniprot \
        --mane \
        --tsl \
        --appris \
        --variant_class \
        --gene_phenotype \
        --mirna \
        --var_synonyms \
        --check_existing \
        --nearest symbol \
        --terms SO \
        --check_existing \
        --clin_sig_allele 1 \
        --force_overwrite \
        --cache \
        --pick \
        --pick_order mane_select,canonical,biotype,rank,appris,tsl,ccds,length \
        --offline \
        --dir_cache /opt/vep/.vep  \
        --dir_plugins /opt/vep/.vep/Plugins/ \
        --plugin AlphaMissense,file=/opt/vep/.vep/vep_annotation/AlphaMissense_hg38.tsv.gz,transcript_match=1 \
        --plugin CADD,{input.CADD} \
        --plugin LoF,loftee_path:/opt/vep/.vep/Plugins/loftee-1.0.4_GRCh38,\
human_ancestor_fa:/opt/vep/.vep/Plugins/loftee-1.0.4_GRCh38/dat/human_ancestor.fa.gz,\
conservation_file:/opt/vep/.vep/Plugins/loftee-1.0.4_GRCh38/dat/loftee.sql,\
gerp_bigwig:/opt/vep/.vep/Plugins/loftee-1.0.4_GRCh38/dat/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
        --custom file={config[clinvar]},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%%CLNDN \
        --plugin SpliceAI,snv=/opt/vep/.vep/vep_annotation/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/opt/vep/.vep/vep_annotation/spliceai_scores.masked.indel.hg38.vcf.gz
        """

rule fam_slivar_filter:
    input:
        vcf = rules.vep.output.vcf,
        ped = "meta/segregation_report_fam.ped"
    output:
        vcf = "families/seg/{chrom}.vcf.gz"
    params:
        slivar_function = "script/slivar-functions.v0.2.8.js"
    threads: 4
    resources:
        nodes=1,
        mem_mb=8000
    container:
        "docker://zihhuafang/slivar_modified:0.2.8"
    shell:
        r"""
        slivar expr \
            --vcf {input.vcf} \
            --ped {input.ped} \
            --js /mnt/slivar/{params.slivar_function} \
            --pass-only \
            --info "variant.ALT[0] != '*' && variant.call_rate >= 0.95" \
            --family-expr "HOM:(variant.CHROM != 'X' && variant.CHROM != 'chrX') && fam.every(segregating_recessive) && INFO.gnomADg_nhomalt <= 1" \
            --family-expr "X_HOM:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_recessive_x) && INFO.gnomADg_nhomalt <= 1" \
            --family-expr "denovo:(variant.CHROM != 'X' && variant.CHROM != 'chrX') && fam.every(segregating_denovo) && INFO.gnomADg_Popmax_AF <= 0.005" \
            --family-expr "x_denovo:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_denovo_x) && INFO.gnomADg_Popmax_AF <= 0.005" \
            --family-expr "dominant:(variant.CHROM != 'X' && variant.CHROM != 'chrX') && fam.every(segregating_dominant) && !fam.every(segregating_denovo) && INFO.gnomADg_Popmax_AF <= 0.005" \
            --family-expr "x_dominant:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_dominant_x) && !fam.every(segregating_denovo_x) && INFO.gnomADg_Popmax_AF <= 0.005" \
            --family-expr "fam_comphet_side:fam.every(function(s) {{return (s.het || s.hom_ref)}}) && fam.some(function(s) {{return s.het && s.affected}}) && INFO.gnomADg_nhomalt <= 1" \
            --trio "comphet_side:comphet_side(kid, mom, dad) && INFO.gnomADg_nhomalt <= 1" \
        | bgzip -c > {output.vcf}

        tabix -p vcf {output.vcf}
        """


skip_list = [
    'intergenic',
    '5_prime_UTR',
    '3_prime_UTR',
    'TFBS_ablation',
    'TF_binding_site',
    'regulatory_region',
    'non_coding_transcript',
    'non_coding',
    'upstream_gene',
    'downstream_gene',
    'non_coding_transcript_exon',
    'NMD_transcript']

rule fam_slivar_comphet:
    input:
        vcf = rules.fam_slivar_filter.output.vcf,
        ped = "meta/segregation_report_fam.ped"
    output:
        vcf = "families/comphet/{chrom}.vcf.gz"
    params:
        skip = ",".join(skip_list)
    container:
        "docker://zihhuafang/slivar_modified:0.2.8"
    threads: 1
    resources:
        nodes = 1,
        mem_mb=10000
    shell:
        r"""
        slivar compound-hets \
               -v {input.vcf} \
               --allow-non-trios --sample-field comphet_side --sample-field fam_comphet_side --sample-field denovo --sample-field x_denovo -p {input.ped} | bgzip -c > {output.vcf}
        """

info_fields = [
    'gnomADg_AF',
    'gnomADg_Popmax_AF',
    'gnomADg_nhomalt',
    'gnomADg_AC',
    'gnomADg_Popmax_AF_filter',
    'TOPMed8_AF']

csq_columns= [
    'Existing_variation',
    'gnomADe_AF',
    'MAX_AF_POPS',
    'LoF',
    'LoF_filter',
    'LoF_flags',
    'LoF_info',
    'CADD_PHRED',
    'SpliceAI_pred_DS_AG',
    'SpliceAI_pred_DS_AL',
    'SpliceAI_pred_DS_DG',
    'SpliceAI_pred_DS_DL',
    'SpliceAI_pred_SYMBOL',
    'Gene',
    'IMPACT',
    'BIOTYPE',
    'STRAND',
    'CANONICAL',
    'MANE_SELECT',
    'NEAREST',
    'EXON',
    'Codons',
    'Amino_acids',
    'HGVSc',
    'HGVSp',
    'ClinVar_CLNSIG',
    'ClinVar_CLNDN',
    'VAR_SYNONYMS']

gene_descs= [
    '/mnt/slivar/pLI_lookup.txt',
    '/mnt/slivar/oe_lof_upper_lookup.txt',
    '/mnt/slivar/oe_mis_upper_lookup.txt',
    '/mnt/slivar/oe_syn_upper_lookup.txt',
    '/mnt/slivar/clinvar_gene_desc.txt',
    '/mnt/slivar/hgncSymbol.inheritance.tsv',
    '/mnt/slivar/gene_fullname_lookup.txt',
    '/mnt/slivar/OMIM_entry_genesymbol.txt']

rule fam_slivar_tsv:
    input:
        ped = "meta/segregation_report_fam.ped",
        vcf = rules.fam_slivar_filter.output.vcf,
    output:
        "families/tmp_report/filtered_{chrom}.tsv"
    params:
        info = "".join([f"--info-field {x} " for x in info_fields]),
        csq = "".join([f"--csq-column {x} " for x in csq_columns]),
        gene_desc = "".join([f"-g {x} " for x in gene_descs]),
        sample_expr = lambda wildcards: (
            "-s x_denovo -s X_HOM -s x_dominant"
            if wildcards.chrom in ["chrX", "X"]
            else "-s denovo -s HOM -s dominant"
        )
    container:
        "docker://zihhuafang/slivar_modified:0.2.8"
    threads: 1
    resources:
        nodes = 1
    envmodules:
        "singularity/3.7.1"
    shell:
        r"""
        slivar tsv \
            {params.info} \
            {params.sample_expr} \
            -c CSQ \
            {params.csq} \
            {params.gene_desc} \
            -p {input.ped} \
            {input.vcf} > {output}
        """


rule fam_comphet_tsv:
    input:
        vcf = rules.fam_slivar_comphet.output.vcf,
        ped = "meta/segregation_report_fam.ped"
    output:
        "families/tmp_report/comphet_{chrom}.tsv"
    params:
        info = "".join([f"--info-field {x} " for x in info_fields]),
        csq = "".join([f"--csq-column {x} " for x in csq_columns]),
        gene_desc = "".join([f"-g {x} " for x in gene_descs])
    container:
        "docker://zihhuafang/slivar_modified:0.2.8"
    threads: 1
    resources:
        nodes = 1
    envmodules:
        "singularity/3.7.1"
    shell:
        r"""
        slivar tsv \
            -s slivar_comphet \
               {params.info} \
               -c CSQ \
               {params.csq} \
               {params.gene_desc} \
             -p {input.ped} \
             {input.vcf} \
             | {{ grep -v ^# || true; }} >> {output}
        """

rule chrom_merge_tsv:
    input:
        filtered = rules.fam_slivar_tsv.output,
        comphet =  rules.fam_comphet_tsv.output
    output:
        "families/tmp_report/merged_{chrom}.tsv"
    params:
        prefix = "families/tmp_report/merged_{chrom}"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        # get header from first file and drop it from other files
        # and make sure slivar_comphet id is unique
        #awk 'NR == FNR || FNR > 1 {{ sub(/^slivar_comphet/, "comphet", $0); print; }}' {input.filtered} {input.comphet}
        ## change slivar_comphet_* to comphet
        cat {input.comphet} | sed -r -e 's/slivar_comphet_[0-9]+/comphet/g' > {params.prefix}_comphet.tsv
        awk 'FNR==1 && NR!=1 {{ next; }} {{ print }}' {input.filtered} > {params.prefix}_filtered.tsv

        cat {params.prefix}_filtered.tsv {params.prefix}_comphet.tsv \
                | sed '1 s/gene_description_1/gnomAD_pLI/;s/gene_description_2/gnomAD_oe_lof_upper/;s/gene_description_3/gnomAD_oe_mis_upper/;s/gene_description_4/gnomAD_oe_syn_upper/;s/gene_description_5/ClinVar_gene_description/;s/gene_description_6/MOI/;s/gene_description_7/Gene_Fullname/;s/gene_description_8/OMIM/;' > {output}
        """

rule report_per_fam:
    input:
        ped = "meta/segregation_report_ped.tsv",
        tsv = "families/tmp_report/merged_{chrom}.tsv"
    output:
        outdir = directory("families/reports_per_chrom/{chrom}")
    threads: 10
    resources:
        mem_mb = 8000
    params:
        prefix = "families/reports_per_chrom/{chrom}"
    script:
        "script/report_format_center.py"

rule gather_fam_report:
    input:
        dirs = lambda wildcards: expand("families/reports_per_chrom/{chrom}", chrom=chr_id)
    output:
        csv = "families/reports/{fam}_segregation_variants.csv"
    run:
        valid_files = []
        for d in input.dirs:
            f_path = os.path.join(d, f"{wildcards.fam}_segregation_variants.csv")
            # Safely check if the file was actually generated by the script
            if os.path.exists(f_path) and os.path.getsize(f_path) > 0:
                valid_files.append(f_path)
        if not valid_files:
            pd.DataFrame().to_csv(output.csv, index=False)
            return
        df = pd.concat(map(pd.read_csv, valid_files), axis=0, ignore_index=True)
        if df.empty:
            df.to_csv(output.csv, index=False)
            return
        comphet = df.loc[df['mode'] == 'comphet']
        other = df.loc[df['mode'] != 'comphet']
        other = other.sort_values(['mode', 'highest_impact'])
        comphet = comphet.sort_values(['gene', 'highest_impact'])
        out = pd.concat([other, comphet], axis=0)
        out.to_csv(output.csv, index=False)

import os
from os.path import join
import yaml

# Global variables
BASE_DIR = workflow.basedir
INPUT_FASTQ_DIR = join(BASE_DIR, "fastq/")
LOG_DIR = join(BASE_DIR, "logs/")
FASTQC_DIR = join(BASE_DIR, "fastqc_reports/")
TRIMMED_FASTQ_DIR = join(INPUT_FASTQ_DIR, "trimmed_polished/")
MAPPED_DIR = join(BASE_DIR, "mapped/")
TRACK_DIR = join(BASE_DIR, "tracks/")

# Configs
configfile: join(BASE_DIR, "config.yaml")
CLUSTER = yaml.load(open(join(BASE_DIR, "cluster.yaml")))

# Get sample names
SAMPLES, = sorted(glob_wildcards(join(INPUT_FASTQ_DIR, '{sample,[^/]+}_R1.fastq.gz')))
DAM = ["Dam"]
FUSION = [sample.split("_")[0] for sample in SAMPLES if not sample.lower().startswith("dam")]
EXP = [sample.split("_")[1] for sample in SAMPLES if not sample.lower().startswith("dam")]

# Create log folder for cluster
os.makedirs(join(LOG_DIR, "cluster/"), exist_ok=True)

# Define targets
FASTQC_RAW_REPORTS = expand(
    join(FASTQC_DIR, "raw/", "{sample}_{read}_fastqc.{ext}"), 
    sample=SAMPLES, 
    read=["R1", "R2"], 
    ext=["zip", "html"])
FASTQC_TRIMMED_POLISHED_REPORTS = expand(
    join(FASTQC_DIR, "trimmed_polished/", "{sample}_{read}_paired.trimmed_polished_fastqc.{ext}"),
    sample=SAMPLES,
    read=["R1", "R2"], 
    ext=["zip", "html"])
BW_FILES = expand(
    join(TRACK_DIR, "{fusion}-vs-Dam.{exp}.bw"), 
    zip, 
    fusion=FUSION, 
    exp=EXP)
BED_GZ_FILES = expand(
    join(TRACK_DIR, "{fusion}-vs-Dam.{exp}.bedgraph.gz"), 
    zip, 
    fusion=FUSION, 
    exp=EXP)
BAI_FILES = expand(
    join(MAPPED_DIR, "{exp}/", "{sample}.{type}.bam.bai"), 
    zip,
    sample=DAM,
    exp=EXP,
    type="filtered") + expand(
    join(MAPPED_DIR, "{exp}/", "{sample}.{type}.bam.bai"), 
    zip,
    sample=FUSION,
    exp=EXP,
    type="filtered")
BAI_FILES = BAI_FILES + [i.replace(".filtered.", ".unfiltered.") for i in BAI_FILES]


# Rules
rule targets:
    input: FASTQC_RAW_REPORTS + FASTQC_TRIMMED_POLISHED_REPORTS + BW_FILES + BAI_FILES + BED_GZ_FILES


rule fastqc_raw:
    input:  
        join(INPUT_FASTQ_DIR, "{sample}_{read}.fastq.gz")
    output:
        expand(join(FASTQC_DIR, "raw/", "{sample}_{read}_fastqc.{ext}"), sample="{sample}", read="{read}", ext=["zip", "html"])
    params:
        out_dir = join(FASTQC_DIR, "raw/")
    log: join(LOG_DIR, "fastqc/raw/{sample}_{read}_fastqc.log")
    # message: """--- QC {sample}_{read} ---"""
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        fastqc \
        --outdir={params.out_dir} \
        {input} 2> {log}
        """


rule adaptor_trimming_polishing_pe:
    input:
        R1 = join(INPUT_FASTQ_DIR, "{sample}_R1.fastq.gz"),
        R2 = join(INPUT_FASTQ_DIR, "{sample}_R2.fastq.gz"),
    output:
        R1_paired = join(TRIMMED_FASTQ_DIR, "{sample}_R1_paired.trimmed_polished.fastq.gz"),
        R2_paired = join(TRIMMED_FASTQ_DIR, "{sample}_R2_paired.trimmed_polished.fastq.gz"),
        R1_unpaired = join(TRIMMED_FASTQ_DIR, "{sample}_R1_unpaired.trimmed_polished.fastq.gz"),
        R2_unpaired = join(TRIMMED_FASTQ_DIR, "{sample}_R2_unpaired.trimmed_polished.fastq.gz")
    params:
        adapter_fa = config["trimming_adaptors_fa_path"]
    log: join(LOG_DIR, "adaptor_trimming_polishing/{sample}_trimmomatic.log")
    #message: """--- Trimming ---""" #"""--- Trimming {sample} ---"""
    threads: CLUSTER["adaptor_trimming_polishing_pe"]["cpu"]
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        trimmomatic PE \
        -threads {threads} \
        -phred33 \
        {input.R1} \
        {input.R2} \
        {output.R1_paired} {output.R1_unpaired} \
        {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:{params.adapter_fa}.fa:2:30:10 \
        TRAILING:28 \
        MINLEN:36 \
        2> {log}
        """


rule fastqc_trimmed:
    input:  
        join(TRIMMED_FASTQ_DIR, "{sample}_{read}_paired.trimmed_polished.fastq.gz")
    output:
        expand(join(FASTQC_DIR, "trimmed_polished/", "{sample}_{read}_paired.trimmed_polished_fastqc.{ext}"), sample="{sample}", read="{read}", ext=["zip", "html"])
    params:
        out_dir = join(FASTQC_DIR, "trimmed_polished/")
    log: join(LOG_DIR, "fastqc/trimmed_polished/{sample}_{read}.trimmed_polished_fastqc.log")
    # message: """--- QC trimmed and polished {sample}_{read} ---"""
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        fastqc \
        --outdir={params.out_dir} \
        {input} 2> {log}
        """


rule mapping_bowtie2_pe:
    input:
        R1 = join(TRIMMED_FASTQ_DIR, "{sample}_{exp}_R1_paired.trimmed_polished.fastq.gz"),
        R2 = join(TRIMMED_FASTQ_DIR, "{sample}_{exp}_R2_paired.trimmed_polished.fastq.gz")
    output:
        join(MAPPED_DIR, "{exp}/", "{sample}.unfiltered.bam")
    params:
        genome_index = config["bowtie2_genome_index_bashpath"]
    log: join(LOG_DIR, "mapping/{sample}_{exp}.mapping.log")
    # message: """--- Mapping {sample}_{exp} ---"""
    threads: CLUSTER["mapping_bowtie2_pe"]["cpu"]
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        bowtie2 \
        -p {threads} \
        -x {params.genome_index} \
        -1 {input.R1} \
        -2 {input.R2} | \
        samtools view -@ {threads} -bhSu - | \
        samtools sort -@ {threads} -T {wildcards.sample}_{wildcards.exp} -O bam - > {output} 2> {log}
        """


rule filter_bam_mapq:
    input:
        join(MAPPED_DIR, "{exp}/", "{sample}.unfiltered.bam")
    output:
        join(MAPPED_DIR, "{exp}/", "{sample}.filtered.bam")
    threads: CLUSTER["filter_bam_mapq"]["cpu"]
    log: join(LOG_DIR, "bam_filtering/{sample}_{exp}.bam_filtering.log")
    # message: """--- Filtering {sample}_{exp} ---"""
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        samtools view -@ {threads} -bhu -q 30 {input} | samtools sort -@ {threads} -T {wildcards.sample}_{wildcards.exp} -O bam -> {output} 2> {log}
        """


rule bam_index:
    input:
        join(MAPPED_DIR, "{exp}/" , "{sample}.{type}.bam")
    output:
        join(MAPPED_DIR, "{exp}/", "{sample}.{type}.bam.bai")
    log: join(LOG_DIR, "bam_filtering/{sample}_{exp}_{type}.bam_filtering.log")
    # message: """--- Filtering {sample}_{exp}_{type} ---"""
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell: 
        """
        samtools index {input} 2> {log}
        """


rule damidseq_analysis:
    input:
        fusion = join(MAPPED_DIR, "{exp}/", "{fusion}.filtered.bam"),
        dam = join(MAPPED_DIR, "{exp}/", "{dam}.filtered.bam")
    output:
        join(TRACK_DIR, "{fusion}-vs-{dam}.{exp}.bedgraph")
    params:
        genome_index = config["bowtie2_genome_index_bashpath"],
        gatc_file = config["gatc_sites_file"],
        out_dir = join(MAPPED_DIR, "{exp}/"),
        temp_bedgraph = join(MAPPED_DIR, "{exp}/", "{fusion}-vs-{dam}.gatc.bedgraph")
    log: join(LOG_DIR, "damidseq/{fusion}_{dam}_{exp}.damidseq.log")
    # message: """--- Damming {fusion} vs. {dam} for {exp} ---"""
    threads: CLUSTER["damidseq_analysis"]["cpu"]
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        cd {params.out_dir}
        damidseq_pipeline \
        --threads={threads} \
        --bowtie2_genome_dir={params.genome_index}  \
        --gatc_frag_file={params.gatc_file} \
        --dam="{input.dam}" "{input.fusion}"
        mv pipeline*.log {log}
        mv {params.temp_bedgraph} {output}
        """


rule bedgraph_to_bigwig:
    input:
        join(TRACK_DIR, "{fusion}-vs-{dam}.{exp}.bedgraph")
    output:
        join(TRACK_DIR, "{fusion}-vs-{dam}.{exp}.bw"),
    params:
        chrom_sizes = config["chrom_sizes"]
    log: join(LOG_DIR, "bedgraph_conversion/{fusion}_{dam}_{exp}.bedgraph_conversion.log")
    # message: """--- Damming {fusion} vs. {dam} for {exp} ---"""
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        bedGraphToBigWig {input} {params.chrom_sizes} {output} 2> {log}
        """

rule call_peaks:
    input:
        join(TRACK_DIR, "{fusion}-vs-{dam}.{exp}.bedgraph")
    output:
        join(TRACK_DIR, "{fusion}-vs-{dam}.{exp}.bedgraph.gz")
    log: join(LOG_DIR, "peak_calling/{fusion}_{dam}_{exp}.peak_calling.log")
    params:
        out_dir = TRACK_DIR
    # message: """--- Calling peaks {fusion} vs. {dam} for {exp} ---"""
    # conda:
    #     "envs/damidseq-dpnii.yaml"
    shell:
        """
        cd {params.out_dir}
        find_peaks {input} 2> {log}
        gzip {input}
        """
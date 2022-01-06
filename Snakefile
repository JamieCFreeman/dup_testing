#### 2021-8-1 JCF ####
# Quick mapping w/ bwa mem to look at library complexity
# Need to make sure to make input files in directory:
# 	"config.yaml"
# 	"units.tsv" (sample table)
# 	as well as directory "envs" with includes specs for conda envs snakemake will use to execute
#  Maps w/ bwa mem, marks duplicates, and calculates some summary stats on different size partitions of the libraries
#  	(b/c I have noticed that in our Tn5 library preps small fragments are the source of more than their share of duplicates

configfile: "config.yaml"

import pandas as pd

from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

# Specify name of results output directory in config file (Snakemake will create) 
OUTDIR = "dup_test"

# Read in sample sheet
samples_table = pd.read_table("units.tsv", dtype=str).set_index("sample", drop=False)

# Get sample wildcards as a list
SAMPLES= samples_table['sample'].values.tolist()

# Functions to return fq file path from sample name
def fq1_from_sample(wildcards):
	return samples_table.loc[wildcards.sample, "fq1"]

def fq2_from_sample(wildcards):
        return samples_table.loc[wildcards.sample, "fq2"]

# Function to fetch read group info from sample table (future- write function to format RG info from sample table
def RG_from_sample(wildcards):
	return samples_table.loc[wildcards.sample, "RG"]

# The final desired files from the pipeline
rule all:
	input:
		expand(f"{OUTDIR}/logs/bwa_mem/{{sample}}_dups.txt", sample=SAMPLES),
		expand(f"{OUTDIR}/logs/bwa_mem/fragment_distrib/{{sample}}_frag_dist.png", sample=SAMPLES),
		expand(f"{OUTDIR}/logs/dupl_partition/{{sample}}_ge151.stats", sample=SAMPLES),
		expand(f"{OUTDIR}/logs/dupl_partition/{{sample}}_le150.stats", sample=SAMPLES),
		expand(f"{OUTDIR}/logs/summary_stats/dupl_partition_summary.txt", sample=SAMPLES)

# Takes only first 12000 reads from each sample through the pipeline (for testing)
#rule test:
#	input:
#		fq1 = fq1_from_sample,
#		fq2 = fq2_from_sample
#	output:
#		cut1 = "test/{sample}_1.fq",
#		cut2 = "test/{sample}_2.fq"
#	shell:
#		"zcat {input.fq1} | head -n 12000 > {output.cut1}; zcat {input.fq2} | head -n 12000 > {output.cut2}"

#rule fastqc:
#	input:
#		fq1 = fq1_from_sample,
#		fq2 = fq2_from_sample
#	output:
#		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.html",
#		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.zip"
#	threads: 4
#	params: 
#		outdir = f"{OUTDIR}/results/fastqc"
#	log: f"{OUTDIR}/logs/fastqc/{{sample}}.log"
#	conda: "envs/fastqc.yaml"	
#	shell:
#		"fastqc --outdir {params.outdir} --format fastq --threads {threads} {input}"

#rule multiqc:
#	input:
#		expand(f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.html")
#	output:
#              expand(f"{OUTDIR}/results/multiqc/{{pre}}_multiqc_report.html",
#                        pre="EF_set1_lib_complexity")
#	params:
#		indir = f"{OUTDIR}/results/fastqc",
#		name = "EF_set1_lib_complexity",
#		outdir = f"{OUTDIR}/results/multiqc"
#	conda:
#		"envs/multiqc.yaml"
#	shell:
#		"multiqc {params.indir} --outdir {params.outdir} --title {params.name}"

#rule multiqc:
#	input:
#              [f"{OUTDIR}/results/qualimap/{sample}/qualimapReport.html" for sample in config["samples"]]
#	output:
#              expand(f"{OUTDIR}/results/multiqc/{{pre}}_multiqc_report.html", 
#			pre=config["prefix"])
#	params:
#		indir = f"{OUTDIR}/results",
#		name = config["prefix"],
#		outdir =  f"{OUTDIR}/results/multiqc"
#	conda: "envs/multiqc.yaml"
#	shell:
#               "multiqc {params.indir} --outdir {params.outdir} --title {params.name}"

rule bwa_mem:
	input:
		fq1 = fq1_from_sample,
		fq2 = fq2_from_sample
	params: 
		REF = config["genome"]
	log:  f"{OUTDIR}/logs/bwa_mem/{{sample}}_map1.log"
	output:
		temp(f"{OUTDIR}/bwa_mem/{{sample}}.bam")
	threads: 8
	shell:
		"bwa mem -M -t {threads} {params.REF} {input.fq1} {input.fq2} | samtools view -bS - > {output}"

rule flagstat:
	input:
		f"{OUTDIR}/bwa_mem/{{sample}}.bam"
	output:
		f"{OUTDIR}/logs/bwa_mem/{{sample}}.stats"
	shell:
		"samtools flagstat {input} > {output}"

rule sort_bam:
	input:
		f"{OUTDIR}/bwa_mem/{{sample}}.bam"
	output:
		temp(f"{OUTDIR}/bwa_mem/{{sample}}_sort.bam")
	shell:
		"samtools sort {input} > {output}"

rule mark_dups:
	input:
		f"{OUTDIR}/bwa_mem/{{sample}}_sort.bam"
	output:
		bam = f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam",
		metrics = f"{OUTDIR}/logs/bwa_mem/{{sample}}_dups.txt"
	conda:  "envs/picard.yaml"
	shell:
		"picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"	
#		"java -Xmx4g -jar picard.jar INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"

rule index_dup:
	input:
		f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam"
	output:
		f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam.bai"
	shell:
		"samtools index {input}"

rule frag_plot:
	input:
		bai = f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam.bai",
		bam = f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam"
	output:	
		plot = f"{OUTDIR}/logs/bwa_mem/fragment_distrib/{{sample}}_frag_dist.png",
		table = f"{OUTDIR}/logs/bwa_mem/fragment_distrib/{{sample}}_frag_dist.txt"
	conda:  "envs/deeptools.yaml"
	shell:
		"bamPEFragmentSize -b {input.bam} -o {output.plot} --table {output.table} --maxFragmentLength 2000"

rule dup_large_frag:
        input:
                bai = f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam.bai",
                bam = f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam"
        output:
                temp(f"{OUTDIR}/bwa_mem/mark_dup/partition/{{sample}}_ge151.bam")
        log:
                f"{OUTDIR}/logs/dupl_partition/{{sample}}_ge151_log.txt"
        conda: "envs/deeptools.yaml"
        shell:
                "alignmentSieve -b {input.bam} -o {output} --minFragmentLength 151 --filterMetrics {log}"

rule dup_small_frag:
        input:
                bai = f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam.bai",
                bam = f"{OUTDIR}/bwa_mem/mark_dup/{{sample}}.bam"
        output: 
                temp(f"{OUTDIR}/bwa_mem/mark_dup/partition/{{sample}}_le150.bam")
        log:
                f"{OUTDIR}/logs/dupl_partition/{{sample}}_le150_log.txt"
        conda: "envs/deeptools.yaml"
        shell:  
                "alignmentSieve -b {input.bam} -o {output} --maxFragmentLength 150 --filterMetrics {log}"

rule flagstat_large:
        input:  
                f"{OUTDIR}/bwa_mem/mark_dup/partition/{{sample}}_ge151.bam"
        output: 
                f"{OUTDIR}/logs/dupl_partition/{{sample}}_ge151.stats"
        shell:
                "samtools flagstat {input} > {output}"

rule flagstat_small:
        input:  
                f"{OUTDIR}/bwa_mem/mark_dup/partition/{{sample}}_le150.bam"
        output: 
                f"{OUTDIR}/logs/dupl_partition/{{sample}}_le150.stats"
        shell:  
                "samtools flagstat {input} > {output}"

rule summary_partion:
	input:
		expand(f"{OUTDIR}/logs/dupl_partition/{{sample}}_ge151.stats", sample=SAMPLES),
                expand(f"{OUTDIR}/logs/dupl_partition/{{sample}}_le150.stats", sample=SAMPLES)
	output:
		f"{OUTDIR}/logs/summary_stats/dupl_partition_summary.txt"
	shell:
		"./helper_scripts/partition_dup_summary.sh"


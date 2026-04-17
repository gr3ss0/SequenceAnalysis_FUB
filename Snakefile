import os

READ_DIR = "tiny"
IDS, = glob_wildcards(READ_DIR + "/{id}.sam")

# files at the end of analysis
rule all:
	input:
		expand("results/stats/{sample}.txt", sample=IDS),
		expand("results/filtered/{sample}_mapping.bam", sample=IDS)

rule convert:
	input:
		READ_DIR+"/{sample}.sam"
	output:
		"results/bam/{sample}.bam"
	shell:
		"samtools view -bS {input} > {output}"

rule sort:
	input:
		"results/bam/{sample}.bam"
	output:
		"results/bam_sorted/{sample}.sorted.bam"
	shell:
		"samtools sort -o {output} {input}"

rule index:
	input:
		"results/bam_sorted/{sample}.sorted.bam"
	output:
		"results/bam_sorted/{sample}.sorted.bam.bai"
	shell:
		"samtools index {input}"

rule calculate_stats:
	input:
		bam = "results/bam_sorted/{sample}.sorted.bam",
		bai = "results/bam_sorted/{sample}.sorted.bam.bai"
	output:
		"results/stats/{sample}.txt"
	shell:
		"samtools idxstats {input.bam} > {output}"
rule extract_mapping:
	input:
		bam = "results/bam_sorted/{sample}.sorted.bam",
		bai = "results/bam_sorted/{sample}.sorted.bam.bai"
	output:
		"results/filtered/{sample}_mapping.bam"
	shell:
		"samtools view -b {input.bam} NZ_AMKI01000040.1 NZ_AMKI01000041.1 > {output}"

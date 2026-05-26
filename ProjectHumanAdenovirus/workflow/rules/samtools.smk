# We don't produce any SAM files.
# rule convert:
# 	input:
# 		"results/sam/{sample}.sam"
# 	output:
# 		"results/bam/{sample}.bam"
# 	threads: 4
# 	log:
# 		"logs/convert/{sample}.log"
# 	conda:
# 		"../envs/mapping.yaml"
# 	shell:
# 		"samtools view -@ {threads} -bS {input} > {output} 2>{log}"

# Sorting is done by minimap2 directly.
# rule sort:
# 	input:
# 		"results/bam/{sample}.bam"
# 	output:
# 		"results/bam_sorted/{sample}_sorted.bam"
# 	threads: 4
# 	log:
# 		"logs/sort/{sample}.log"
# 	conda:
# 		"../envs/mapping.yaml"
# 	shell:
# 		"samtools sort -@ {threads} -o {output} {input} 2>{log}"

rule index:
	input:
		"results/bam_sorted/{sample}_sorted.bam"
	output:
		"results/bam_sorted/{sample}_sorted.bam.bai"
	threads: 4
	log:
		"logs/index/{sample}.log"
	conda:
		"../envs/mapping.yaml"
	shell:
		"samtools index -@ {threads} {input} 2>{log}"

rule calculate_stats:
	input:
		bam = "results/bam_sorted/{sample}_sorted.bam",
		bai = "results/bam_sorted/{sample}_sorted.bam.bai"
	output:
		"results/stats/{sample}.stats"
	threads: 4
	log:
		"logs/stats/{sample}.log"
	conda:
		"../envs/mapping.yaml"
	shell:
		"samtools idxstats -@ {threads} {input.bam} > {output} 2>{log}"

# rule extract_mapping:
# 	input:
# 		bam = "results/bam_sorted/{sample}_sorted.bam",
# 		bai = "results/bam_sorted/{sample}_sorted.bam.bai"
# 	output:
# 		"results/filtered/{sample}_mapping.bam"
# 	threads: 4
# 	log:
# 		"logs/filter/{sample}.log"
# 	conda:
# 		"../envs/mapping.yaml"
# 	shell:
# 		"samtools view -@ {threads} -b {input.bam} NZ_AMKI01000040.1 NZ_AMKI01000041.1 > {output} 2>{log}"


rule samtools_pileup:
	input:
		bamlist = expand("results/mapped/{sample}_aln.sorted.bam", sample=SAMPLES.index),


	
	
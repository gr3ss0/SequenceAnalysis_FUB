rule convert:
	input:
		"results/sam/{sample}.sam"
	output:
		"results/bam/{sample}.bam"
	log:
		"logs/convert/{sample}.log"
	shell:
		"samtools view -bS {input} > {output} 2>{log}"

rule sort:
	input:
		"results/bam/{sample}.bam"
	output:
		"results/bam_sorted/{sample}_sorted.bam"
	log:
		"logs/sort/{sample}.log"
	shell:
		"samtools sort -o {output} {input} 2>{log}"

rule index:
	input:
		"results/bam_sorted/{sample}_sorted.bam"
	output:
		"results/bam_sorted/{sample}_sorted.bam.bai"
	log:
		"logs/index/{sample}.log"
	shell:
		"samtools index {input} 2>{log}"

rule calculate_stats:
	input:
		bam = "results/bam_sorted/{sample}_sorted.bam",
		bai = "results/bam_sorted/{sample}_sorted.bam.bai"
	output:
		"results/stats/{sample}.stats"
	log:
		"logs/stats/{sample}.log"
	shell:
		"samtools idxstats {input.bam} > {output} 2>{log}"
rule extract_mapping:
	input:
		bam = "results/bam_sorted/{sample}_sorted.bam",
		bai = "results/bam_sorted/{sample}_sorted.bam.bai"
	output:
		"results/filtered/{sample}_mapping.bam"
	log:
		"logs/filter/{sample}.log"
	shell:
		"samtools view -b {input.bam} NZ_AMKI01000040.1 NZ_AMKI01000041.1 > {output} 2>{log}"

rule convert:
	input:
		"results/assembly/polish/{sample}_aligned_{read}.sam"
	output:
		"results/assembly/polish/{sample}_aligned_{read}.bam"
	threads: 4
	log:
		"logs/polish/convert/{sample}.log"
	conda:
		"../envs/assembly.yaml"
	shell:
		"samtools view -@ {threads} -bS {input} > {output} 2>{log}"



rule merge_map:
    input:
        sam1 = "results/assembly/polish/{sample}_aligned_1.sam",
        sam2 = "results/assembly/polish/{sample}_aligned_2.sam"
    output:
        paired_bam = temp("results/assembly/polish/{sample}/recombined_paired.bam")
    threads: 4
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        # 1. Sort both SAMs by read name (-n) and convert to BAM
        samtools sort -n -@ {threads} -o tmp1_map.bam {input.sam1}
        samtools sort -n -@ {threads} -o tmp2_map.bam {input.sam2}
        
        # 2. Merge them back into a single paired-end BAM file
        samtools merge -@ {threads} -n {output.paired_bam} tmp1_map.bam tmp2_map.bam
        
        # 3. Clean up temporary files
        rm tmp1_map.bam tmp2_map.bam
        """

rule sort_map:
	input:
		rules.merge_map.output.paired_bam,
	output:
		bam="results/assembly/polish/{sample}/recombined_paired_sorted.bam",
	threads: 4
	log:
		"logs/polish/sort/{sample}.log"
	conda:
		"../envs/assembly.yaml"
	shell:
		"samtools sort -@ {threads} -o {output.bam} {input} 2>{log}"


rule index_map:
	input:
		rules.sort_map.output.bam
	output:
		bai="results/assembly/polish/{sample}/recombined_paired_sorted.bam.bai"
	threads: 4
	log:
		"logs/polish/index/{sample}.log"
	conda:
		"../envs/assembly.yaml"
	shell:
		"samtools index -@ {threads} {input} 2>{log}"

rule merge_filtered:
    input:
        sam1 = "results/assembly/polish/{sample}_filtered_1.sam",
        sam2 = "results/assembly/polish/{sample}_filtered_2.sam"
    output:
        paired_bam = temp("results/assembly/polish/{sample}/recombined_paired_filtered.bam")
    threads: 4
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        # 1. Sort both SAMs by read name (-n) and convert to BAM
        samtools sort -n -@ {threads} -o tmp1_filter.bam {input.sam1}
        samtools sort -n -@ {threads} -o tmp2_filter.bam {input.sam2}
        
        # 2. Merge them back into a single paired-end BAM file
        samtools merge -@ {threads} -n {output.paired_bam} tmp1_filter.bam tmp2_filter.bam
        
        # 3. Clean up temporary files
        rm tmp1_filter.bam tmp2_filter.bam
        """
rule sort_filtered:
    input:
        rules.merge_filtered.output.paired_bam,
    output:
        bam="results/assembly/polish/{sample}/recombined_paired_filtered_sorted.bam",
    threads: 4
    log:
        "logs/polish/sort_filtered/{sample}.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2>{log}"

rule index_filtered:
    input:
        rules.sort_filtered.output.bam
    output:
        bai="results/assembly/polish/{sample}/recombined_paired_filtered_sorted.bam.bai"
    threads: 4
    log:
        "logs/polish/index_filtered/{sample}.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        "samtools index -@ {threads} {input} 2>{log}"

rule calculate_stats:
    input:
        bam = "results/bam_sorted/{sample}_sorted.bam",
        bai = "results/bam_sorted/{sample}_sorted.bam.bai"
    output:
        idxstats = "results/stats/{sample}.idxstats",
        flagstat = "results/stats/{sample}.flagstat",
        stats = "results/stats/{sample}.stats"
    threads: 2
    log:
        "logs/stats/{sample}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats} 2>{log}
        samtools flagstat {input.bam} > {output.flagstat} 2>>{log}
        samtools stats {input.bam} > {output.stats} 2>>{log}
        """


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



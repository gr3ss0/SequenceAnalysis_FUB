rule convert:
	input:
		"results/assembly/polish/{sample}_aligned_{read}.sam"
	output:
		"results/assembly/polish/{sample}_aligned_{read}.bam"
	threads: 4
	log:
		"logs/polish/convert/{sample}_{read}.log"
	conda:
		"../envs/assembly.yaml"
	shell:
		"samtools view -@ {threads} -bS {input} > {output} 2>{log}"


rule recombine_and_pair:
    input:
        sam1 = "results/assembly/polish/{sample}_{step}_1.sam",
        sam2 = "results/assembly/polish/{sample}_{step}_2.sam"
    output:
        paired_bam = temp("results/assembly/polish/{sample}/recombined_paired_{step}.bam")
    threads: 4
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        # Create unique temporary names based on the output file path
        TMP1="{output.paired_bam}.tmp1.bam"
        TMP2="{output.paired_bam}.tmp2.bam"

        # 1. Sort both SAMs by read name (-n) and convert to BAM
        samtools sort -n -@ {threads} -o "$TMP1" {input.sam1}
        samtools sort -n -@ {threads} -o "$TMP2" {input.sam2}
        
        # 2. Merge them back into a single paired-end BAM file
        samtools merge -@ {threads} -n {output.paired_bam} "$TMP1" "$TMP2"
        
        # 3. Clean up temporary files safely
        rm "$TMP1" "$TMP2"
        """
# rule recombine_and_pair:
#     input:
#         sam1 = "results/assembly/polish/{sample}_{step}_1.sam",
#         sam2 = "results/assembly/polish/{sample}_{step}_2.sam"
#     output:
#         paired_bam = temp("results/assembly/polish/{sample}/recombined_paired_{step}.bam")
#     threads: 4
#     conda:
#         "../envs/assembly.yaml"
#     shell:
#         """
#         # 1. Sort both SAMs by read name (-n) and convert to BAM
#         samtools sort -n -@ {threads} -o tmp1.bam {input.sam1}
#         samtools sort -n -@ {threads} -o tmp2.bam {input.sam2}
        
#         # 2. Merge them back into a single paired-end BAM file
#         samtools merge -@ {threads} -n {output.paired_bam} tmp1.bam tmp2.bam
        
#         # 3. Clean up temporary files
#         rm tmp1.bam tmp2.bam
#         """

rule sort:
	input:
		rules.recombine_and_pair.output.paired_bam,
	output:
		bam = temp(rules.recombine_and_pair.output.paired_bam + ".sorted")
	threads: 4
	log:
		"logs/polish/sort/{sample}_{step}.log"
	conda:
		"../envs/assembly.yaml"
	shell:
		"samtools sort -@ {threads} -o {output.bam} {input} 2>{log}"


rule index:
	input:
		rules.sort.output.bam
	output:
		temp(rules.sort.output.bam +".bai")
	threads: 4
	log:
		"logs/polish/index/{sample}_{step}.log"
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



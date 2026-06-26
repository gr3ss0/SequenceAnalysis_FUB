
rule sort_by_name:
    input:
        "results/star/pe/{sample}/pe_aligned.bam",
    output:
        "results/star/pe/{sample}/{sample}_aligned.sorted.bam"
    threads:8
    log:
        "logs/sort_by_name/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -n -@ {threads} -o {output} {input} 2>{log}"

# rule index:
# 	input:
# 		"results/bam_sorted/{sample}_sorted.bam"
# 	output:
# 		"results/bam_sorted/{sample}_sorted.bam.bai"
# 	threads: 4
# 	log:
# 		"logs/index/{sample}.log"
# 	conda:
# 		"../envs/mapping.yaml"
# 	shell:
# 		"samtools index -@ {threads} {input} 2>{log}"

# rule calculate_stats:
#     input:
#         bam = "results/bam_sorted/{sample}_sorted.bam",
#         bai = "results/bam_sorted/{sample}_sorted.bam.bai"
#     output:
#         idxstats = "results/stats/{sample}.idxstats",
#         flagstat = "results/stats/{sample}.flagstat",
#         stats = "results/stats/{sample}.stats"
#     threads: 2
#     log:
#         "logs/stats/{sample}.log"
#     conda:
#         "../envs/mapping.yaml"
#     shell:
#         """
#         samtools idxstats {input.bam} > {output.idxstats} 2>{log}
#         samtools flagstat {input.bam} > {output.flagstat} 2>>{log}
#         samtools stats {input.bam} > {output.stats} 2>>{log}
#         """

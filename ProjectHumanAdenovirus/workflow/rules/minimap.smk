configfile: "config/config.yaml"


# https://snakemake-wrappers.readthedocs.io/en/v7.8.0/wrappers/bio/minimap2/index.html
rule minimap2_index:
    input:
        target=config["ref"]
    output:
        "results/index/reference.mmi"
    log:
        "logs/minimap2_index/ref.log"
    params:
        extra=""  # optional additional args
    threads: 4
    wrapper:
        "v7.1.0/bio/minimap2/index"


# # https://snakemake-wrappers.readthedocs.io/en/v7.8.0/wrappers/bio/minimap2/aligner.html
# rule minimap2_bam_sorted:
#     input:
#         target = rules.minimap2_index.output,  # REF, can be either genome index or genome fasta
#         query = ["query/reads1.fasta", "query/reads2.fasta"], 
#     output:
#         "results/mapped/{sample}_aln.sorted.bam",
#     log:
#         "logs/minimap2/{sample}.log",
#     params:
#         extra="-x map-pb",  # optional
#         sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
#         sort_extra="",  # optional: extra arguments for samtools/picard
#     threads: 3
#     wrapper:
#         "v7.6.0/bio/minimap2/aligner"


def get_map_input(wildcards):
    # Case 1: Trimming was skipped, use raw data from SAMPLES dataframe
    if config["analysis_options"].get("skip_trimming", False):
        return {
            "r1": SAMPLES.at[wildcards.sample, 'fq1'],
            "r2": SAMPLES.at[wildcards.sample, 'fq2']
        }
    
    # Case 2: Trimming was performed, use results from your trimming rule
    # Note: Use the actual filenames produced by your Trimmomatic rule
    else:
        return {
            "r1": f"results/trimmed/{wildcards.sample}_R1_paired.fq.gz",
            "r2": f"results/trimmed/{wildcards.sample}_R2_paired.fq.gz"
        }


rule minimap2:
    input:
        unpack(get_map_input), 
        ref_index = rules.minimap2_index.output,  # REF, can be either genome index or genome fasta
        
    output:
        "results/bam_sorted/{sample}_sorted.bam",
    log:
        "logs/minimap2/{sample}.log",
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        "minimap2 -t {threads} -ax sr {input.ref_index} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output} - "

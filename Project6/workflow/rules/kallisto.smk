# !!! Watch out: Index and Quant has to be run by same version of Kallisto, otherwise the quantification will fail.

# https://snakemake-wrappers.readthedocs.io/en/v3.13.0/wrappers/kallisto/index.html
rule kallisto_index:
    input:
        fasta="results/trinity_out_dir.Trinity.fasta",
    output:
        index="results/kallisto/assembly.idx",
    params:
        extra="",  # optional parameters
    log:
        "logs/kallisto_index.log",
    threads: 30
    wrapper:
        "v9.13.0/bio/kallisto/index"

# https://snakemake-wrappers.readthedocs.io/en/v3.13.0/wrappers/kallisto/quant.html
rule kallisto_quant:
    input:
        fastq=lambda wildcards: [
            get_map_input(wildcards)["fq1"], 
            get_map_input(wildcards)["fq2"]
        ],
        index=rules.kallisto_index.output.index,
    output:
        directory("results/kallisto/quant_{sample}"),
    params:
        extra="",
    log:
        "logs/kallisto/quant_{sample}.log",
    threads: 1
    wrapper:
        "v9.13.0/bio/kallisto/quant"
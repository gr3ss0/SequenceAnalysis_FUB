# https://snakemake-wrappers.readthedocs.io/en/v3.13.0/wrappers/kallisto/index.html
rule kallisto_index:
    input:
        fasta="{transcriptome}.fasta",
    output:
        index="{transcriptome}.idx",
    params:
        extra="",  # optional parameters
    log:
        "logs/kallisto_index_{transcriptome}.log",
    threads: 1
    wrapper:
        "v3.13.0/bio/kallisto/index"
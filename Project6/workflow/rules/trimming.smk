# https://snakemake-wrappers.readthedocs.io/en/v7.8.0/wrappers/bio/fastp.html
rule fastp_pe:
    input:
        sample = lambda wildcards: [
            SAMPLES.loc[wildcards.sample, 'fq1'],
            SAMPLES.loc[wildcards.sample, 'fq2']
        ]
    output:
        trimmed=[
            "results/trimmed/{sample}.1.fastq.gz", 
            "results/trimmed/{sample}.2.fastq.gz"
        ],
        unpaired1="results/trimmed/{sample}.u1.fastq.gz",
        unpaired2="results/trimmed/{sample}.u2.fastq.gz",
        failed="results/trimmed/{sample}.failed.fastq.gz",
        html="results/report/{sample}.html",
        json="results/report/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
        
    threads: 4
    params:
        extra = lambda wildcards, threads: f"--thread {threads}"
    wrapper:
        "v7.1.0/bio/fastp"

rule trim_all:
    input:
        expand("results/trimmed/{sample}.1.fastq.gz", sample=SAMPLES.index),
        expand("results/trimmed/{sample}.2.fastq.gz", sample=SAMPLES.index),
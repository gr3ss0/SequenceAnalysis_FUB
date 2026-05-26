configfile: "config/config.yaml"

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastqc.html
rule run_raw_qc_per_file:
    input:
        lambda wildcards: SAMPLES.at[wildcards.sample, wildcards.read]
    output:
        html="results/qc/fastqc/raw/{sample}_{read}.html",
        zip="results/qc/fastqc/raw/{sample}_{read}_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/raw/{sample}_{read}.log",
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.0/bio/fastqc"

# https://snakemake-wrappers.readthedocs.io/en/v7.8.0/wrappers/bio/fastp.html
rule fastp_pe:
    input:
        sample=["reads/pe/{sample}.1.fastq", "reads/pe/{sample}.2.fastq"]
    output:
        trimmed=["results/trimmed/{sample}.1.fastq", "results/trimmed/{sample}.2.fastq"],
        # Unpaired reads separately
        unpaired1="results/trimmed/pe/{sample}.u1.fastq",
        unpaired2="results/trimmed/{sample}.u2.fastq",
        merged="results/trimmed/{sample}.merged.fastq",
        failed="results/trimmed/{sample}.failed.fastq",
        html="results/report/{sample}.html",
        json="results/report/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 4
    wrapper:
        "v7.1.0/bio/fastp"


# rule trimmomatic:
#     input:
#         r1 = lambda wildcards:SAMPLES.at[wildcards.sample,'fq1'],
#         r2 = lambda wildcards:SAMPLES.at[wildcards.sample,'fq2']
#     output:
#         r1 = "results/trimmed/{sample}_1.fastq.gz", #trimmed R1 fastq(gz) file 
#         r2 = "results/trimmed/{sample}_2.fastq.gz", #trimmed R2 fastq(gz) file (if PE)
#         r1_unpaired = "results/trimmed/{sample}_1_unpaired.fastq.gz", #unpaired R1 fastq(gz) file (if PE)
#         r2_unpaired = "results/trimmed/{sample}_2_unpaired.fastq.gz" #unpaired R2 fastq(gz) file (if PE)
#     log:
#         "logs/trimmomatic/{sample}.log"
#     conda:
#         "../envs/mapping.yaml"
#     threads: 4
#     params:
#         trim_mode=config["trimmomatic_params"]["mode"],
#         adapter=config["trimmomatic_params"]["adapter_file"],
#         seed_mismatch=config["trimmomatic_params"]["seed_mismatch"],
#         palindrome=config["trimmomatic_params"]["palindrome_treshold"],
#         simple=config["trimmomatic_params"]["simple_treshold"],
#         min_adapter=config["trimmomatic_params"]["min_adapter_length"],
#         keep_reads=config["trimmomatic_params"]["keep_both_reads"]

#     shell:
#         """
#         trimmomatic {params.trim_mode} \
#         -threads {threads} \
#         {input} {output} \
#         ILLUMINACLIP:{params.adapter}:{params.seed_mismatch}:{params.palindrome}:{params.simple}:{params.min_adapter}:{params.keep_reads} \
#         > {log} 2>&1
#         """


# # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastqc.html
# rule run_coocked_qc:
#     input:
#         lambda wildcards: (
#             [f"results/trimmed/{wildcards.sample}_1.fastq.gz"] if wildcards.read == "fq1" 
#             else [f"results/trimmed/{wildcards.sample}_2.fastq.gz"]
#         ) if not config["analysis_options"]["skip_trimming"]==True else []
#     output:
#         html="results/qc/fastqc/processed/{sample}_{read}.html",
#         zip="results/qc/fastqc/processed/{sample}_{read}_fastqc.zip"
#     params:
#         extra="--quiet",
#         mem_overhead_factor=0.1,
#     log:
#         "logs/fastqc/processed/{sample}_{read}.log",
#     threads: 1
#     resources:
#         mem_mb = 1024,
#     wrapper:
#         "v7.6.0/bio/fastqc"


rule qualimap:
    input:
        # fails if not sorted
        bam="results/bam_sorted/{sample}_sorted.bam",
        bai="results/bam_sorted/{sample}_sorted.bam.bai"
    output:
        directory("results/qc/qualimap/{sample}")
    log:
        "logs/qualimap/bamqc/{sample}.log",
    conda:
        "../envs/mapping.yaml"
    threads: 4
    shell:
        """
        qualimap bamqc -nt {threads} \
        -bam {input.bam} \
        -outdir {output} \
        > {log} 2>&1
        """


rule multiqc_all:
    input:
        expand("results/qc/qualimap/{sample}", sample=SAMPLES.index) if not config["analysis_options"]["skip_qualimap"]==True else [],
        expand("results/qc/fastqc/processed/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=['fq1', 'fq2'])
    output:
        report_file="results/qc/multiqc_all.html",
        out_dir=directory("results/qc/multiqc_all_data")
    log:
        "logs/multiqc/all.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        multiqc {input} \
            --filename multiqc_all.html \
            --outdir results/qc \
            > {log} 2>&1
        """


rule run_raw_qc:
    input:
        expand("results/qc/fastqc/raw/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=['fq1', 'fq2'])
    output:
        # Definierte Pfade relativ zum Projektverzeichnis
        report_file="results/qc/multiqc_raw.html",
        out_dir=directory("results/qc/multiqc_raw_data") 
    log:
        "logs/multiqc/raw.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        multiqc {input} \
            --filename multiqc_raw.html \
            --outdir results/qc \
            > {log} 2>&1
        """
    
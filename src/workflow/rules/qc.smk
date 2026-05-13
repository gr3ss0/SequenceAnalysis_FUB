

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



# https://snakemake-wrappers.readthedocs.io/en/v3.3.1/wrappers/trimmomatic/se.html
rule trimmomatic:
    input:
        r1 = lambda wildcards:SAMPLES.at[wildcards.sample,'fq1'],
        r2 = lambda wildcards:SAMPLES.at[wildcards.sample,'fq2']
    output:
        r1 = "results/trimmed/{sample}_1.fastq.gz", #trimmed R1 fastq(gz) file 
        r2 = "results/trimmed/{sample}_2.fastq.gz", #trimmed R2 fastq(gz) file (if PE)
        r1_unpaired = "results/trimmed/{sample}_1_unpaired.fastq.gz", #unpaired R1 fastq(gz) file (if PE)
        r2_unpaired = "results/trimmed/{sample}_2_unpaired.fastq.gz" #unpaired R2 fastq(gz) file (if PE)
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        # optional compression levels from -0 to -9 and -11
        compression_level="-9"
    threads: 4
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "v9.8.0/bio/trimmomatic"

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastqc.html
rule run_coocked_qc:
    input:
        lambda wildcards: "results/trimmed/{sample}_1.fastq.gz" if wildcards.read == "fq1" else "results/trimmed/{sample}_2.fastq.gz"
    output:
        html="results/qc/fastqc/processed/{sample}_{read}.html",
        zip="results/qc/fastqc/processed/{sample}_{read}_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/processed/{sample}_{read}.log",
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.0/bio/fastqc"


# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/qualimap/bamqc.html
rule qualimap:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam="results/bam_sorted/{sample}_sorted.bam",
        bai="results/bam_sorted/{sample}_sorted.bam.bai"
    output:
        directory("results/qc/qualimap/{sample}")
    log:
        "logs/qualimap/bamqc/{sample}.log",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=4096,
    wrapper:
        "v7.6.0/bio/qualimap/bamqc"


rule multiqc_all:
    input:
        expand("results/qc/qualimap/{sample}", sample=SAMPLES.index),
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
    
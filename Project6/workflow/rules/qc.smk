configfile: "config/config.yaml"

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastqc.html
rule raw_qc_per_read:
    input:
        lambda wildcards: SAMPLES.loc[wildcards.sample, f"fq{wildcards.read}"]
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

rule run_raw_qc:
    input:
        expand("results/qc/fastqc/raw/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=['1', '2'])
    output:
        # Definierte Pfade relativ zum Projektverzeichnis
        report_file="results/qc/multiqc_raw.html",
        out_dir=directory("results/qc/multiqc_raw_data") 
    log:
        "logs/multiqc/raw.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {input} \
            --filename multiqc_raw.html \
            --outdir results/qc \
            > {log} 2>&1
        """
    

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastqc.html
rule run_coocked_qc:
    input:
        lambda wildcards: (
            [] if config["analysis_options"].get("skip_trimming", False)
            else [f"results/trimmed/{wildcards.sample}.{wildcards.read}.fastq.gz"]
        )
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

rule qualimap:
    input:
        # sorted by name (see samtools::sort_by_name)
        bam="results/star/pe/{sample}/{sample}_aligned.sorted.bam", 
        annotation=config["annotation"],
    output:
        directory("results/qc/qualimap/{sample}"),
    log:
        "logs/qualimap/bamqc/{sample}.log",
    conda:
        "../envs/qc.yaml"
    threads: 1
    resources:
        mem_mb=14000 
    shell:
        """
        qualimap rnaseq \
        -bam {input.bam} \
        -gtf {input.annotation} \
        -outdir {output} \
        --paired \
        --sorted \
        --java-mem-size=12G \
        > {log} 2>&1
        """



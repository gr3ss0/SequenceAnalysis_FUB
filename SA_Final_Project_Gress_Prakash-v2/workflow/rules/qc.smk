rule raw_qc:
    input:
        lambda wildcards: SAMPLES.at[wildcards.sample, f"short_fq{wildcards.read}"]
    output:
        html="results/qc/fastqc/raw_short/{sample}_{read}.html",
        zip="results/qc/fastqc/raw_short/{sample}_{read}_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/raw_short/{sample}_{read}.log",
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.0/bio/fastqc"

rule raw_qc_long:
    input:
        lambda wildcards: SAMPLES.at[wildcards.sample, "long_fq"]
    output:
        html="results/qc/fastqc/raw_long/{sample}.html",
        zip="results/qc/fastqc/raw_long/{sample}_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/raw_long/{sample}.log",
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.0/bio/fastqc"



# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastqc.html
rule processed_qc_short:
    input:
        lambda wildcards: (
            [] if config["analysis_options"].get("skip_trimming", True)
            else [f"results/trimmed/short/{wildcards.sample}.{wildcards.read}.fastq.gz"]
        )
    output:
        html="results/qc/fastqc/processed_short/{sample}_{read}.html",
        zip="results/qc/fastqc/processed_short/{sample}_{read}_fastqc.zip"
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

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastqc.html
rule processed_qc_long:
    input:
        "results/trimmed/long/{sample}.fastq.gz"
    output:
        html="results/qc/fastqc/processed_long/{sample}.html",
        zip="results/qc/fastqc/processed_long/{sample}_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/processed/{sample}_long.log",
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.0/bio/fastqc"

rule qualimap_polish_map: #this is before the polypolish filter
    input:
        # fails if not sorted
        bam="results/assembly/polish/{sample}/recombined_paired_sorted.bam",
        bai="results/assembly/polish/{sample}/recombined_paired_sorted.bam.bai"
    output:
        directory("results/qc/qualimap/polish/{sample}")
    log:
        "logs/qualimap/bamqc/{sample}_polish.log",
    conda:
        "../envs/qc.yaml"
    threads: 4
    shell:
        """
        qualimap bamqc -nt {threads} \
        -bam {input.bam} \
        -outdir {output} \
        > {log} 2>&1
        """

rule qualimap_polish_filter: #this is after the polypolish filter
    input:
        # fails if not sorted
        bam="results/assembly/polish/{sample}/recombined_paired_filtered_sorted.bam",
        bai="results/assembly/polish/{sample}/recombined_paired_filtered_sorted.bam.bai"
    output:
        directory("results/qc/qualimap/filter/{sample}")
    log:
        "logs/qualimap/bamqc/{sample}_filter.log",
    conda:
        "../envs/qc.yaml"
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
        *(
            expand("results/qc/fastqc/processed_short/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=["1", "2"])
            if TRIMMING_ENABLED
            else expand("results/qc/fastqc/raw_short/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=["1", "2"])
        ),
        *(
            expand("results/qc/fastqc/processed_long/{sample}_fastqc.zip", sample=SAMPLES.index)
            if TRIMMING_ENABLED
            else expand("results/qc/fastqc/raw_long/{sample}_fastqc.zip", sample=SAMPLES.index)
        ),
        *(
            expand("results/qc/qualimap/polish/{sample}", sample=SAMPLES.index)
            if (QUALIMAP_ENABLED and INTEGRATE_SHORT_READS)
            else []
        ),
        *(
            expand("results/qc/qualimap/filter/{sample}", sample=SAMPLES.index)
            if (QUALIMAP_ENABLED and INTEGRATE_SHORT_READS)
            else []
        )
        
    output:
        report_file="results/qc/multiqc_all.html",
        out_dir=directory("results/qc/multiqc_all_data")

    log:
        "logs/multiqc/all.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {input} \
            --filename multiqc_all.html \
            --outdir results/qc \
            > {log} 2>&1
        """


rule run_raw_qc:
    input:
        expand("results/qc/fastqc/raw_short/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=['1', '2']),
        expand("results/qc/fastqc/raw_long/{sample}_fastqc.zip", sample=SAMPLES.index) if len(SAMPLES.index)>0 else []
    output:
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
    
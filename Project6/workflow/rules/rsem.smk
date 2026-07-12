rule unzip_reads:
    input:
        unpack(get_map_input)
    output:
        fq1=temp("results/unzipped/{sample}.1.fastq"),
        fq2=temp("results/unzipped/{sample}.2.fastq")
    threads: 1
    shell:
        """
        gzip -d -c {input.fq1} > {output.fq1}
        gzip -d -c {input.fq2} > {output.fq2}
        """


rule prepare_reference:
    input:
        reference_genome="results/trinity_out_dir.Trinity.fasta",
    output:
        multiext("results/rsem/index/reference", ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"),
        # Bowtie2 index files
        bowtie2=multiext("results/rsem/index/reference", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    threads: 16
    log:
        "logs/rsem/prepare-reference.log"
    conda:
        "../envs/rsem.yaml"
    shell:
        # We target the output directory prefix dynamically
        "rsem-prepare-reference --bowtie2 -p {threads} {input.reference_genome} results/rsem/index/reference > {log} 2>&1"


rule calculate_expression:
    input:
        fq1="results/unzipped/{sample}.1.fastq",
        fq2="results/unzipped/{sample}.2.fastq",
        bowtie2=multiext("results/rsem/index/reference", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    output:
        genes_results="results/rsem/quant/{sample}.genes.results",
        isoforms_results="results/rsem/quant/{sample}.isoforms.results"
    threads: 16
    log:
        "logs/rsem/calculate_expression/{sample}.log"
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-calculate-expression \
            --bowtie2 \
            --paired-end \
            -p {threads} \
            --seed 42 \
            {input.fq1} \
            {input.fq2} \
            results/rsem/index/reference \
            results/rsem/quant/{wildcards.sample} > {log} 2>&1
        """


rule rsem_generate_data_matrix:
    input:
        expand("results/rsem/quant/{sample}.genes.results", sample=SAMPLES.index)
    output:
        "results/rsem/all.results"
    log:
        "logs/rsem/generate_data_matrix.log"
    conda:
        "../envs/rsem.yaml"
    shell:
        "rsem-generate-data-matrix {input} > {output} 2> {log}"
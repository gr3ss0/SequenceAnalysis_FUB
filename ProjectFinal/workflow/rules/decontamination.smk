CONTAMINATION_FASTA = config.get("decontamination", {}).get("contamination_fasta", "")


rule kraken2_short:
    input:
        unpack(decide_trimming_short), 
        db = config["kraken2_db"]
    output:
        report = "results/kraken2/short/{sample}.kraken2.report.txt"
    log:
        "logs/kraken2/{sample}_short.log"
    threads: 8
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {input.db} \
            --threads {threads} \
            --paired \
            --report {output.report} \
            --output - \
            {input.r1} {input.r2} \
            > /dev/null 2> {log}
        """

rule kraken2_long:
    input:
        get_map_input_long, 
        db = config["kraken2_db"]
    output:
        report = "results/kraken2/long/{sample}.kraken2.report.txt"
    log:
        "logs/kraken2/{sample}_long.log"
    threads: 8
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {input.db} \
            --threads {threads} \
            --report {output.report} \
            --output - \
            {input}\
            > /dev/null 2> {log}
        """

rule screen_long:
    input:
        expand("results/kraken2/long/{sample}.kraken2.report.txt", sample=SAMPLES.index)

rule multiqc_screen:
    input:
        expand("results/kraken2/long/{sample}.kraken2.report.txt", sample=SAMPLES.index),
        expand("results/kraken2/short/{sample}.kraken2.report.txt", sample=SAMPLES.index),
    output:
        report_file = "results/qc/multiqc_screen.html",
        out_dir = directory("results/qc/multiqc_screen_data")
    log:
        "logs/multiqc/screen.log"
    conda:
        "../envs/multiqc.yaml" 
    shell:
        """
        multiqc {input} \
            --filename multiqc_screen.html \
            --outdir results/qc \
            > {log} 2>&1
        """

rule screen:
    # this rule just takes as input the html produced by multi_qc_screen
    #   snakemake --use-conda --cores 10 screen can be used to stop here and inspect the multiqc report
    input:
        "results/qc/multiqc_screen.html"


rule decon_index: #building the bowtie2 index for the mapping to the contamination sequences
    input:
        fasta = CONTAMINATION_FASTA
    output:
        multiext(
            "results/decon_index/contamination",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
            ".rev.1.bt2", ".rev.2.bt2"
        )
    params:
        prefix = "results/decon_index/contamination"
    log:
        "logs/decon_index/build.log"
    threads: 64
    conda:
        "../envs/mapping.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1"


rule decon_map_short:
    input:
        unpack(decide_trimming_short), #the input to minimap in 4B is input to decon_map in 5A
        index = rules.decon_index.output,
    output:
        bam = "results/decontamination/short/{sample}_contamination_mapped.bam"
    params:
        prefix = "results/decon_index/contamination"
    log:
        "logs/decon_map/short/{sample}.log"
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        """

        bowtie2 -p {threads} -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} 2> {log} | \
        samtools view -@ {threads} -bS - > {output.bam} 2>> {log}
        """


rule decon_filter_short: # this is the filtering, the output of this rule is the input to the minimap rule if decontamination is enabled, see function get_minimap_input
    input:
        bam = rules.decon_map_short.output.bam,
    output:
        r1 = "results/decontaminated/short/{sample}.1.fastq",
        r2 = "results/decontaminated/short/{sample}.2.fastq",
        name_sorted = temp("results/decontamination/short/{sample}_unmapped_namesorted.bam")
    log:
        "logs/decon_filter/{sample}.log"
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        set -e -o pipefail
        
        # 12 is that both the reads in a pair should be unmapped for it to be kept as a read pair for downstream assembly
        samtools view -@ {threads} -b -f 12 {input.bam} 2> {log} | \
        samtools sort -@ {threads} -n -o {output.name_sorted} - >> {log} 2>&1

        samtools fastq -@ {threads} \
            -1 {output.r1} -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            {output.name_sorted} >> {log} 2>&1
        """

    # Note: there is a folder called decontamination and another called decontaminated
    # They differ in that decontaminated contains the filtered output and decontamination contains everything else

rule unzip_long:
    input:
        reads = decide_trimming_long,
    output:
        reads = temp("results/trimmed/long/{sample}.fastq")
    log:
        "logs/unzip_long/{sample}.log"
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gunzip -c {input.reads} > {output.reads} 2> {log}
        """

rule decon_map_long:
    input:
        reads = rules.unzip_long.output.reads, 
        ref = rules.decon_index.output  # Reference FASTA or .mmi index file
    output:
        bam = "results/decontamination/long/{sample}_contamination_mapped.bam"
    log:
        "logs/decon_map/long/{sample}.log"
    threads: 4
    params:
        preset = "map-ont"  # Use "map-hifi" for PacBio HiFi
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        minimap2 -ax {params.preset} -t {threads} {input.ref} {input.reads} 2> {log} | \
        samtools view -@ {threads} -bS - > {output.bam} 2>> {log}
        """


rule decon_filter_long:
    input:
        bam = rules.decon_map_long.output.bam,
    output:
        fastq = "results/decontaminated/long/{sample}.fastq"
    log:
        "logs/decon_filter/long/{sample}.log"
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        set -e -o pipefail

        # -f 4 keeps unmapped single-end reads
        samtools view -@ {threads} -b -f 4 {input.bam} 2> {log} | \
        samtools fastq -@ {threads} - > {output.fastq} 2>> {log}
        """

rule decon_stats_short: 
    #this is included in rule multiqc_all when decontamination is enabled (can be seen in qc.smk)
    input:
        bam = rules.decon_map_short.output.bam, #only doing flagstats not idxtstats so no need of bai, apart from that it is similar to the rule calculate_stats in samtools.smk
    output:
        flagstat = "results/decontamination/short/{sample}_contamination.flagstat"
    log:
        "logs/decon_stats/short/{sample}.log"
    threads: 1
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} 2>{log}"

rule decon_stats_long: 
    #this is included in rule multiqc_all when decontamination is enabled (can be seen in qc.smk)
    input:
        bam = rules.decon_map_long.output.bam, #only doing flagstats not idxtstats so no need of bai, apart from that it is similar to the rule calculate_stats in samtools.smk
    output:
        flagstat = "results/decontamination/long/{sample}_contamination.flagstat"
    log:
        "logs/decon_stats/long/{sample}.log"
    threads: 1
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} 2>{log}"

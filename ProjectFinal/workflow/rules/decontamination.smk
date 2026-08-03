rule kraken2_short:
    input:
        unpack(get_decon_input_short), #it should not get the decontaminated reads as input
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
        reads = get_decon_input_long, #it should not get the decontaminated reads as input
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
        kraken2 \
            --db {input.db} \
            --threads {threads} \
            --report {output.report} \
            --output /dev/null \
            {input.reads} \
            2> {log}
        """


rule multiqc_screen:
    input:
        expand("results/kraken2/long/{sample}.kraken2.report.txt", sample=SAMPLES_LONG.index) if len(SAMPLES_LONG.index) > 0 else [],
        expand("results/kraken2/short/{sample}.kraken2.report.txt", sample=SAMPLES_SHORT.index) if len(SAMPLES_SHORT.index) > 0 else [],
    output:
        report_file = "results/qc/multiqc_screen.html",
        out_dir = directory("results/qc/multiqc_screen_data")
    log:
        "logs/multiqc/screen.log"
    conda:
        "../envs/qc.yaml" 
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



#rule contaminants_index:
#     input:
#        target=config["contamination_fasta"]
#     output:
#         index="results/index/reference.mmi"
#     log:
#         "logs/minimap2_index/ref.log"
#     threads: 4
#     conda:
#         "../envs/mapping.yaml"
#     shell:
#         "minimap2 -t {threads} -d {output.index} {input.target} > {log} 2>&1"


rule decon_index: #building the bowtie2 index for the mapping to the contamination sequences
    input:
        fasta = config['contamination_fasta']
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
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1"


rule decon_map:
    input:
        unpack(get_decon_input_short), 
        index = rules.decon_index.output
    output:
        bam = "results/decontamination/{sample}_contamination_mapped.bam"
    params:
        prefix = "results/decon_index/contamination"
    log:
        "logs/decon_map/{sample}.log"
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        set -e -o pipefail

        bowtie2 \
            -p {threads} \
            -x {params.prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            2>{log} |

        samtools view -@ {threads} -bS - > {output.bam} 2>>{log}
        """


rule decon_filter: 
    input:
        bam = "results/decontamination/{sample}_contamination_mapped.bam"
    output:
        r1 = "results/decontaminated/{sample}.1.fastq",
        r2 = "results/decontaminated/{sample}.2.fastq",
        name_sorted = temp("results/decontamination/{sample}_unmapped_namesorted.bam")
    log:
        "logs/decon_filter/{sample}.log"
    threads: 8
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

rule decon_index_long: #building the minimap2 index for mapping long reads to the contamination sequences
    input:
        fasta = config['contamination_fasta']
    output:
        mmi = "results/decon_index_long/contamination.mmi"
    log:
        "logs/decon_index_long/build.log"
    threads: 8
    params:
        preset = config.get("minimap2_long_preset", "map-pb")
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        minimap2 \
            -t {threads} \
            -x {params.preset} \
            -d {output.mmi} \
            {input.fasta} \
            > {log} 2>&1
        """


rule decon_map_long:
    input:
        reads = get_decon_input_long,
        index = rules.decon_index_long.output.mmi #this is mmi file
    output:
        bam = "results/decontamination/long/{sample}_contamination_mapped.bam"
    log:
        "logs/decon_map_long/{sample}.log"
    threads: 8
    params:
        preset = config.get("minimap2_long_preset", "map-pb")
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        set -e -o pipefail

        minimap2 \
            -t {threads} \
            -ax {params.preset} \
            {input.index} \
            {input.reads} \
            2>{log} |

        samtools view \
            -@ {threads} \
            -b \
            -o {output.bam} \
            - \
            2>> {log}
        """

rule decon_filter_long:
    input:
        bam = rules.decon_map_long.output.bam
    output:
        fq = "results/decontaminated/long/{sample}.fastq", 
        
    log:
        "logs/decon_filter_long/{sample}.log"
    threads: 8
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools fastq \
            -@ {threads} \
            -f 4 \
            -n \
            {input.bam} \
            > {output.fq} \
            2> {log}
        """
        

#     rule decon_stats: 
#         #this is included in rule multiqc_all when decontamination is enabled (can be seen in qc.smk)
#         input:
#             bam = "results/decontamination/{sample}_contamination_mapped.bam" #only doing flagstats not idxtstats so no need of bai, apart from that it is similar to the rule calculate_stats in samtools.smk
#         output:
#             flagstat = "results/decontamination/{sample}_contamination.flagstat"
#         log:
#             "logs/decon_stats/{sample}.log"
#         threads: 1
#         conda:
#             "../envs/mapping.yaml"
#         shell:
#             "samtools flagstat {input.bam} > {output.flagstat} 2>{log}"



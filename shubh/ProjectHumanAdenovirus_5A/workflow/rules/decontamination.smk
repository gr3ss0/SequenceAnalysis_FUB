# Path to the (optional) contamination multi-FASTA, "" if not enabled
CONTAMINATION_FASTA = config.get("decontamination", {}).get("contamination_fasta", "")


def contamination_enabled():
    """True if a contamination fasta was configured -> decontamination on."""
    return CONTAMINATION_FASTA not in (None, "", "none", "None")


def get_minimap_input(wildcards):
    """
    Reads that are used as input for the mapping based assembly
    (rule minimap2 in rules/minimap.smk).
    
    """
    if contamination_enabled():
        return {
            "r1": f"results/decontaminated/{wildcards.sample}.1.fastq",
            "r2": f"results/decontaminated/{wildcards.sample}.2.fastq",
        }
    return get_map_input(wildcards)


# -------------------------------------------------------------------
# 1) Mandatory taxonomic screening with Kraken2
# -------------------------------------------------------------------


rule kraken2:
    input:
        unpack(get_map_input),
        db = config["kraken2_db"]
    output:
        report = "results/kraken2/{sample}.kraken2.report.txt"
    log:
        "logs/kraken2/{sample}.log"
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


rule multiqc_screen:
    input:
        expand("results/kraken2/{sample}.kraken2.report.txt", sample=SAMPLES.index)
    output:
        report_file = "results/qc/multiqc_screen.html",
        out_dir = directory("results/qc/multiqc_screen_data")
    log:
        "logs/multiqc/screen.log"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        multiqc {input} \
            --filename multiqc_screen.html \
            --outdir results/qc \
            > {log} 2>&1
        """


rule screen:
    # Checkpoint target:
    #   snakemake --use-conda --cores 10 screen
    #
    # Runs QC + Kraken2 for every sample and collects all reports in a
    # single MultiQC report, so the user can inspect it and decide
    # whether (and which) sequences belong in a contamination fasta
    # (config["decontamination"]["contamination_fasta"]).
    input:
        "results/qc/multiqc_screen.html"


# -------------------------------------------------------------------
# 2) Optional decontamination (mapping + filtering)
#
# Only enters the DAG if contamination_enabled() == True, because the
# only consumer of these outputs is get_minimap_input().
# -------------------------------------------------------------------

# https://snakemake-wrappers.readthedocs.io/ (plain bowtie2, same tool
# as used elsewhere in this workflow, see envs/mapping.yaml)
rule decon_index:
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
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1"


rule decon_map:
    input:
        unpack(get_map_input),
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
        bowtie2 -p {threads} -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} 2> {log} | \
        samtools view -@ {threads} -bS - > {output.bam} 2>> {log}
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
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        set -e -o pipefail
        
        # Keep only read pairs where BOTH mates are unmapped against the
        # contamination reference, i.e. reads that do NOT belong to a
        # contaminant (samtools flag 4 = read unmapped, 8 = mate unmapped).
        samtools view -@ {threads} -b -f 12 {input.bam} 2> {log} | \
        samtools sort -@ {threads} -n -o {output.name_sorted} - >> {log} 2>&1

        samtools fastq -@ {threads} \
            -1 {output.r1} -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            {output.name_sorted} >> {log} 2>&1
        """


rule decon_stats:
    # quick flagstat of the decontamination mapping, e.g. to see how many
    # read pairs were classified as contamination; included in
    # multiqc_all when decontamination is enabled (see rules/qc.smk)
    input:
        bam = "results/decontamination/{sample}_contamination_mapped.bam"
    output:
        flagstat = "results/decontamination/{sample}_contamination.flagstat"
    log:
        "logs/decon_stats/{sample}.log"
    threads: 1
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} 2> {log}"

rule star_index:
    input:
        fasta=config["reference"],
        gtf=config["annotation"],
    output:
        directory("resources/star_genome"),
    log:
        "logs/star_index_genome.log",
    params:
        extra=""
    threads: 36
    wrapper:
        "v7.2.0/bio/star/index"



def get_map_input(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return {
            "fq1": SAMPLES.at[wildcards.sample, 'fq1'], 
            "fq2": SAMPLES.at[wildcards.sample, 'fq2']
        }
    else:
        # Updated to point to fastp's compressed outputs
        return {
            "fq1": f"results/trimmed/{wildcards.sample}.1.fastq.gz", 
            "fq2": f"results/trimmed/{wildcards.sample}.2.fastq.gz"
        }


# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/star/align.html
rule star_pe_multi:
    input:
        unpack(get_map_input),
        idx=config.get("precomputed_ref_index", "resources/star_genome"), # if precomputed ref_index provided --> skip star_index
    output:
        aln="results/star/pe/{sample}/{sample}_aligned.bam",
        log="logs/pe/{sample}/Log.out",
        sj="results/star/pe/{sample}/SJ.out.tab",
        unmapped=[
            "results/star/pe/{sample}/unmapped.1.fastq.gz",
            "results/star/pe/{sample}/unmapped.2.fastq.gz"
        ],
    log:
        "logs/pe/{sample}.log",
    params:
        extra="--outSAMtype BAM Unsorted",
    threads: 12
    wrapper:
        "v9.4.2/bio/star/align"
# https://snakemake-wrappers.readthedocs.io/en/0.70.0/wrappers/subread/featurecounts.html

rule feature_counts:
    input:
        sam=expand("results/star/pe/{sample}/{sample}_aligned.bam", sample=SAMPLES.index), # list of sam or bam files
        annotation=config["annotation"],
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        multiext("results/featureCounts/all",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        2
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O -p --fracOverlap 0.2" # -p: libraries are assumed to contain paired-end reads
    log:
        "logs/featureCounts/all.log"
    wrapper:
        "0.70.0/bio/subread/featurecounts"
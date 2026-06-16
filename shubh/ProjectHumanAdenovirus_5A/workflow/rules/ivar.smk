rule ivar_consensus:
    input:
        ref_fa = config["ref"],
        ref_index = "results/index/reference.mmi",
        sorted_bam = "results/bam_sorted/{sample}_sorted.bam",
        sorted_bai = "results/bam_sorted/{sample}_sorted.bam.bai"
    output:
        "results/consensus/{sample}_consensus.fa",
        "results/consensus/{sample}_consensus.qual.txt"
    log:
        "logs/ivar/{sample}.log"
    threads: 4
    params:
        prefix = lambda wildcards, output: output[0].replace(".fa", ""),
        max_depth=config["ivar_params"].get("max_depth", 60000),
        min_freq=config["ivar_params"].get("min_freq", 0.6),
        min_qual=config["ivar_params"].get("min_qual", 20),
        min_cov=config["ivar_params"].get("min_cov", 10)

    conda:
        "../envs/ivar.yaml"
    shell:
        """
        (samtools mpileup -A -d {params.max_depth} -Q {params.min_qual} -q {params.min_qual} -f {input.ref_fa} {input.sorted_bam} | \
        ivar consensus -p {params.prefix} -q {params.min_qual} -t {params.min_freq} -m {params.min_cov} -n N) > {log} 2>&1
        """

rule concat_consensus:
    input:
        expand("results/consensus/{sample}_consensus.fa", sample=SAMPLES.index)
    output:
        "results/consensus/combined.fa"
    log:
        "logs/concat_consensus.log"
    
    run:
        from pathlib import Path

        with open(output[0], "w") as out:
            for fasta in input:
                sample = Path(fasta).name.replace("_consensus.fa", "")

                with open(fasta) as fin:
                    for line in fin:
                        if line.startswith(">"):
                            out.write(f">{sample}\n")
                        else:
                            out.write(line)
    


rule msa:
    input:
        "results/consensus/combined.fa"
    output:
        "results/phylo_tree/msa.fa"
    log:
        "logs/msa.log"
    threads: 8
    resources:
        mem_mb=8000  # Gives the job 8GB of RAM to prevent crashes
    conda:
        "../envs/phylo.yaml"
    shell:
        "mafft --thread {threads} --nomemsave --auto {input} > {output} 2> {log}"


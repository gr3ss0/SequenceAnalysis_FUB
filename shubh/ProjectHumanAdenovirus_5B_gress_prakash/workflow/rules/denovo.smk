

rule spades_assembly:
    input:
        unpack(get_minimap_input) #all the reads (potentially filtered by trimming, decontamination)
    output:
        contigs = "results/spades/{sample}/contigs.fasta" #just 1 of the full directory of outputs
    params:
        outdir = "results/spades/{sample}",
        memory = config["assembly_params"]["spades"].get("memory_gb", 16),
        extra = config["assembly_params"]["spades"].get("extra", "")
    log:
        "logs/spades/{sample}.log"
    threads: 8
    conda:
        "../envs/denovo.yaml"
    shell:
        """
        spades.py --phred-offset 33 --pe1-1 {input.r1} --pe1-2 {input.r2} \
            -o {params.outdir} -t {threads} -m {params.memory} \
            {params.extra} > {log} 2>&1 
        """

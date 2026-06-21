

rule denovo_bowtie2_index:
    input:
        contigs = "results/spades/{sample}/contigs.fasta" #this is coming from spades, index is made of the reference, in this case they are the contigs
    output:
        multiext(
            "results/denovo_index/{sample}/contigs",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
            ".rev.1.bt2", ".rev.2.bt2"
        )
    params:
        prefix = "results/denovo_index/{sample}/contigs"
    log:
        "logs/denovo_index/{sample}.log"
    threads: 4
    conda:
        "../envs/mapping.yaml" 
    shell:
        "bowtie2-build --threads {threads} {input.contigs} {params.prefix} > {log} 2>&1"


rule denovo_bowtie2_map:
    input:
        unpack(get_minimap_input), # the same reads that were assembled by spades_assembly
        index = rules.denovo_bowtie2_index.output #from previous rule
    output:
        bam = "results/denovo_mapped/{sample}_sorted.bam"
    params:
        prefix = "results/denovo_index/{sample}/contigs",
        extra = config["assembly_params"]["polishing"].get("bowtie2_extra", "")
    log:
        "logs/denovo_map/{sample}.log"
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bowtie2 -p {threads} {params.extra} -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} 2> {log} | \
        samtools sort -@ {threads} -o {output.bam} - >> {log} 2>&1
        """


rule denovo_bam_index: #needed for the pilon polishing
    input:
        "results/denovo_mapped/{sample}_sorted.bam"
    output:
        "results/denovo_mapped/{sample}_sorted.bam.bai"
    threads: 4
    log:
        "logs/denovo_bam_index/{sample}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index -@ {threads} {input} 2>{log}"


rule pilon_polish:
    input:
        assembly = "results/spades/{sample}/contigs.fasta",
        bam = "results/denovo_mapped/{sample}_sorted.bam",
        bai = "results/denovo_mapped/{sample}_sorted.bam.bai"
    output:
        "results/polishing/{sample}/{sample}_polished.fasta" #this is the one used in Quast, not the unpolished one. Also used in the blasting
    params:
        outdir = "results/polishing/{sample}",
        prefix = "{sample}_polished"
    log:
        "logs/pilon/{sample}.log"
    threads: 4
    conda:
        "../envs/pilon.yaml"
    shell:
        """
        pilon --genome {input.assembly} --frags {input.bam} \
            --output {params.prefix} --outdir {params.outdir} \
            --threads {threads} > {log} 2>&1
        """

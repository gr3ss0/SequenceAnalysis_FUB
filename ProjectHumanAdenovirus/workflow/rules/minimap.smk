configfile: "config/config.yaml"


rule minimap2_index:
    input:
        target=config["ref"]
    output:
        index="results/index/reference.mmi"
    log:
        "logs/minimap2_index/ref.log"
    threads: 4
    conda:
        "../envs/mapping.yaml"
    shell:
        "minimap2 -t {threads} -d {output.index} {input.target} > {log} 2>&1"


# # https://snakemake-wrappers.readthedocs.io/en/v7.8.0/wrappers/bio/minimap2/aligner.html
# rule minimap2_bam_sorted:
#     input:
#         target = rules.minimap2_index.output,  # REF, can be either genome index or genome fasta
#         query = ["query/reads1.fasta", "query/reads2.fasta"], 
#     output:
#         "results/mapped/{sample}_aln.sorted.bam",
#     log:
#         "logs/minimap2/{sample}.log",
#     params:
#         extra="-x map-pb",  # optional
#         sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
#         sort_extra="",  # optional: extra arguments for samtools/picard
#     threads: 3
#     wrapper:
#         "v7.6.0/bio/minimap2/aligner"




rule minimap2:
    input:
        unpack(get_map_input), 
        ref_index = rules.minimap2_index.output,  # REF, can be either genome index or genome fasta
        
    output:
        bam = "results/bam_sorted/{sample}_sorted.bam",
    log:
        "logs/minimap2/{sample}.log",
    threads: 4
    params:
        kmer_size=15,
        seed_density=5,
        mismatch_penalty=3,
        gap_open_penalty="4,16",
        gap_extension_penalty="2,1",
    conda:
        "../envs/mapping.yaml"
    # shell:
    #     """
    #     minimap2 -ax sr -k 15 -w 5 -B 3 -O 4,16 -E 2,1 -t {threads} {input.ref_index} {input.r1} {input.r2} 2>{log} | \
    #     samtools sort -@ {threads} -o {output.bam} - >>{log} 2>&1
    #     """
    
    shell:
        """
        minimap2 -t {threads} -ax sr {input.ref_index} {input.r1} {input.r2} 2>{log} | \
        samtools sort -@ {threads} -o {output.bam} - >>{log} 2>&1
        """



rule ragtag_patch: #this will give the multifasta file with lots of sequences
    input:
        target = "results/polishing/{sample}/{sample}_polished.fasta",     
        reference = "results/reference_selection/{sample}/best_reference.fasta" 
    output:
        "results/patching/{sample}/ragtag.patch.fasta" #just one of the many files produced, this is the multifasta
    params:
        outdir = "results/patching/{sample}",
        min_len = config["assembly_params"]["ragtag_patch"].get("min_len", 2000),  # -s
        max_dist = config["assembly_params"]["ragtag_patch"].get("max_dist", 0.2),       # -i
        extra = config["assembly_params"]["ragtag_patch"].get("extra", "")
    log:
        "logs/ragtag_patch/{sample}.log"
    threads: 4
    conda:
        "../envs/ragtag.yaml"
    shell:
        """
        ragtag.py patch -o {params.outdir} -t {threads} \
            -s {params.min_len} -i {params.max_dist} {params.extra} \
            {input.target} {input.reference} > {log} 2>&1
        """


rule select_patched_sequence:
    input:
        patched = "results/patching/{sample}/ragtag.patch.fasta" #multifasta
    output:
        "results/patched_reference/{sample}/patched_reference.fasta" #the longest only (single fasta)
    log:
        "logs/select_patched_sequence/{sample}.log"
    conda:
        "../envs/phylo.yaml" 
    script:
        "../scripts/select_patched_sequence.py" #script new to 5B

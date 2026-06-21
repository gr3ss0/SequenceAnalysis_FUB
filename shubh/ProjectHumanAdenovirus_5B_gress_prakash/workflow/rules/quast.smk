

rule quast_setup_gridss: #this rule takes no input, it just downloads the gridss
    output:
        touch("results/quast/.gridss_ready") #flag file which is created when rule executes
    log:
        "logs/quast/gridss_setup.log"
    conda:
        "../envs/quast.yaml"
    shell:
        "quast-download-gridss > {log} 2>&1"


rule quast:
    input:
        unpack(get_minimap_input), # same reads used for de-novo assembly/polishing
        assembly = "results/polishing/{sample}/{sample}_polished.fasta", #only doing on the polished ones
        reference = "results/reference_selection/{sample}/best_reference.fasta", #to find SV
        gridss_ready = "results/quast/.gridss_ready" #flag file needed for finding the SV
    output:
        directory("results/qc/quast/{sample}")
    params:
        extra = config["assembly_params"]["quast"].get("extra", "")
    log:
        "logs/quast/{sample}.log"
    threads: 4
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast.py {input.assembly} \
            -r {input.reference} \
            -1 {input.r1} -2 {input.r2} \
            -o {output} \
            --threads {threads} \
            {params.extra} > {log} 2>&1
        """

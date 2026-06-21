rule kraken2:
    input:
        unpack(get_map_input), #what was the input to minimap in 4B, is input to Kraken in 5A, this can be trimmed or untrimmed reads depending on the config for trimming
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
        "../envs/mapping.yaml" 
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


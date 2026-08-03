CARD_DIR= config["card_db"]
CARD_VERSION = config["card_version"]


rule set_up_card_download:
    output:
        json = CARD_DIR + "/card.json",
        pmid = CARD_DIR + "/PMID.tsv"
    log:
        "logs/card/download.log"
    params:
        version = CARD_VERSION,
        out_dir = CARD_DIR,
        tarball = CARD_DIR + "/card-data.tar.bz2"
    conda:
        "../envs/resistance_prediction.yaml"
    shell:
        """
        wget -O {params.tarball} "https://card.mcmaster.ca/download/0/broadstreet-v{params.version}.tar.bz2" > {log} 2>&1
        tar -xf {params.tarball} -C {params.out_dir} >> {log} 2>&1 
        """

rule set_up_card_after_download:
    input:
        card_json = CARD_DIR + "/card.json", 
        pmid = CARD_DIR + "/PMID.tsv"
    output:
        dummy= CARD_DIR + "/.card_loaded"
    log:
        "logs/card/setup.log"
    conda:
        "../envs/resistance_prediction.yaml"
    shell:
        """
        rgi load \
            --card_json {input.card_json} \
            --local \
            > {log} 2>&1

        touch {output.dummy}
        """

rule card_amr_detection:
    input:
        card_ready = rules.set_up_card_after_download.output.dummy,
        proteins = "results/annotation/{sample}/{sample}.faa",
    output:
        txt = "results/card_amr_report/{sample}/card_amr_report.txt",
        json = "results/card_amr_report/{sample}/card_amr_report.json"
    log:
        "logs/rgi/{sample}.log"
    threads: 4
    conda:
        "../envs/resistance_prediction.yaml"
    params:
        prefix = "results/card_amr_report/{sample}/card_amr_report",
        alignment_tool = "blast" 
    shell:
        """
        rgi main \
            --input_sequence {input.proteins} \
            --output_file {params.prefix} \
            --input_type protein \
            --alignment_tool {params.alignment_tool} \
            -n {threads} \
            --local \
            --clean \
            > {log} 2>&1
        """

rule card_amr_detection_external:
    input:
        card_ready = rules.set_up_card_after_download.output.dummy,
        proteins = "results/annotation_ext/external_genome/external_genome.faa",
    output:
        txt = "results/card_amr_report/external/card_amr_report.txt",
        json = "results/card_amr_report/external/card_amr_report.json"
    log:
        "logs/rgi/external.log"
    threads: 4
    conda:
        "../envs/resistance_prediction.yaml"
    params:
        prefix = "results/card_amr_report/external/card_amr_report",
        alignment_tool = "blast"
    shell:
        """
        rgi main \
            --input_sequence {input.proteins} \
            --output_file {params.prefix} \
            --input_type protein \
            --alignment_tool {params.alignment_tool} \
            -n {threads} \
            --local \
            --clean \
            > {log} 2>&1
        """
rule merge_amr:
    input:
        reports = get_card_amr_reports
    output:
        xlsx = "results/card_amr_report/amr_merged.xlsx"
    conda:
        "../envs/openpyxl.yaml"
    script:
        "../scripts/merge_amr_excel.py"
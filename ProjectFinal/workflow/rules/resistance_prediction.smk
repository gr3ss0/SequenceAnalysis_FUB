CARD_DIR= config["card_db"]
CARD_VERSION = config["card_version"]


CARD_AMR_REPORTS = expand(
    "results/card_amr_report/{sample}/card_amr_report.txt",
    sample=SAMPLES_LONG.index
)

if config["external_genome"]:
    CARD_AMR_REPORTS.append(
        "results/card_amr_report/external/card_amr_report.txt"
    )
rule set_up_card:
    output:
        tarball = temp(CARD_DIR+ "/card-data.tar.bz2"),
        json = CARD_DIR + "/card.json",
        #db_instance = directory(CARD_DIR),
    log:
        "logs/card/download.log"
    params:
        version = CARD_VERSION,
        out_dir = CARD_DIR
    conda:
        "../envs/resistance_prediction.yaml"
    shell:
        """
        wget -O {output.tarball} "https://card.mcmaster.ca/download/0/broadstreet-v{params.version}.tar.bz2" > {log} 2>&1
        tar -xf {output.tarball} -C {params.out_dir} >> {log} 2>&1
        rgi load --card_json {output.json} --local  >> {log} 2>&1
        """


rule card_amr_detection:
    input:
        rules.set_up_card.output.json,
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
        alignment_tool = "blast" # Or "diamond" for even faster protein alignments
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
if config["external_genome"]:
    rule card_amr_detection_external:
        input:
            rules.set_up_card.output.json,
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
            alignment_tool = "blast" # Or "diamond" for even faster protein alignments
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
        reports = CARD_AMR_REPORTS
    output:
        xlsx = "results/card_amr_report/amr_merged.xlsx"
    conda:
        "../envs/openpyxl.yaml" # Ensure pandas and openpyxl are in this env
    script:
        "../scripts/merge_amr_excel.py"
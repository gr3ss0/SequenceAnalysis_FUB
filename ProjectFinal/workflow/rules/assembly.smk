rule long_based_assembly:
    input:
        get_map_input_long
    output:
        assembly = "results/assembly/{sample}/primary_assembly.fasta",
        #out_dir = directory("results/assembly/{sample}"),
    threads: 16
    log:
        "logs/flye/{sample}.log"
    params:
        technique=config.get("flye_technique", "--pacbio-raw"),
        out_dir = directory("results/assembly/{sample}"),
    conda:
        "../envs/flye.yaml"
    shell:
        """
        flye {params.technique} {input} --out-dir {params.out_dir} --threads {threads} > {log} 2>&1
        mv "{params.out_dir}/assembly.fasta" {output.assembly}
        """

rule align_shorts:
    input:
        unpack(get_map_input_short),
        assembly = rules.long_based_assembly.output.assembly,
    output:
        index = expand("results/assembly/{{sample}}/primary_assembly.fasta.{ext}",
               ext=["amb","ann","bwt","pac","sa"]),
        aligned1 = "results/assembly/polish/{sample}_aligned_1.sam",
        aligned2 = "results/assembly/polish/{sample}_aligned_2.sam",
    threads: 16
    log:
        "logs/bwa/{sample}.log"
    conda:
        "../envs/polypolish.yaml"
    shell:
        # the reads are mapped separately according to polypolish manual
        """
        bwa index {input.assembly} > {log} 2>&1
        bwa mem -t {threads} -a {input.assembly} {input.r1} > {output.aligned1} 2>> {log}
        bwa mem -t {threads} -a {input.assembly} {input.r2} > {output.aligned2} 2>> {log}
        """
    
rule filter_by_insert_size:
    input:
        aligned1 = rules.align_shorts.output.aligned1,
        aligned2 = rules.align_shorts.output.aligned2,
    output:
        filtered1 = "results/assembly/polish/{sample}_filtered_1.sam",
        filtered2 = "results/assembly/polish/{sample}_filtered_2.sam",
    log:
        "logs/polypolish/filter/{sample}.log"
    conda:
        "../envs/polypolish.yaml"
    shell:
        """polypolish filter \
        --in1 {input.aligned1} \
        --in2 {input.aligned2} \
        --out1 {output.filtered1} \
        --out2 {output.filtered2} 2>{log}
        """

rule short_based_polish:
    input:
        r1 = rules.filter_by_insert_size.output.filtered1,
        r2 = rules.filter_by_insert_size.output.filtered2,
        primary_assembly = rules.long_based_assembly.output.assembly,

    output:
        contig = "results/assembly/{sample}/secondary_assembly.fasta",
    log:
        "logs/polypolish/polish/{sample}.log"
    threads: 1
    conda:
        "../envs/polypolish.yaml"
    shell:
        "polypolish polish {input.primary_assembly} {input.r1} {input.r2} > {output.contig} 2>{log}"

#rule get_secondary_assembly:
    #input:
        #expand("results/assembly/{sample}/secondary_assembly.fasta", sample=SAMPLES_LONG.index),

rule quast:
    input:
        contig = f"results/assembly/{{sample}}/{FINAL_ASSEMBLY}",
    output:
        html = "results/assembly/{sample}/quast/report.html",
        out_dir = directory("results/assembly/{sample}/quast"),
    log:
        "logs/quast/{sample}.log"
    threads: 16
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast.py \
            {input.contig} \
            -o {output.out_dir} \
            --threads {threads} \
            > {log} 2>&1
        """

rule busco:
    input:
        contig = f"results/assembly/{{sample}}/{FINAL_ASSEMBLY}",

    output:
        out_dir=directory("results/assembly/{sample}/busco")

    log:
        "logs/busco/{sample}.log"

    threads: 16

    params:
        lineage=config.get("busco_lineage", "enterobacterales_odb12")

    conda:
        "../envs/busco.yaml"

    shell:
        """
        busco \
            -i {input.contig} \
            -m genome \
            -l {params.lineage} \
            -o busco \
            --out_path results/assembly/{wildcards.sample} \
            --cpu {threads} \
            > {log} 2>&1
        """

rule multiqc_quast_busco:
    input:
        expand(
            "results/assembly/{sample}/quast",
            sample=SAMPLES_LONG.index
        ),
        expand(
            "results/assembly/{sample}/busco",
            sample=SAMPLES_LONG.index
        )

    output:
        report_file="results/assembly/multiqc_quast_busco.html",
        out_dir=directory("results/assembly/multiqc_quast_busco_data")

    log:
        "logs/multiqc/quast_busco.log"

    conda:
        "../envs/qc.yaml"

    shell:
        """
        multiqc {input} \
            --filename multiqc_quast_busco.html \
            --outdir results/assembly \
            > {log} 2>&1
        """
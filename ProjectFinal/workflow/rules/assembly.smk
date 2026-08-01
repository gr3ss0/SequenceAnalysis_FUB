rule flye_long_based_assembly:
    input:
        decide_contamination_long
    output:
        assembly = "results/assembly/{sample}/primary_assembly.fasta",
        #out_dir = directory("results/assembly/{sample}"),
    threads: 16
    log:
        "logs/flye/{sample}.log"
    params:
        technique=config['assembly_options']['long_seq_method'],
        out_dir = directory("results/assembly/{sample}"),
        polishing_iterations = config['assembly_options']['long_polish_iterations']
    conda:
        "../envs/flye.yaml"
    shell:
        """
        flye {params.technique} {input} --iterations {params.polishing_iterations} --out-dir {params.out_dir} --threads {threads} > {log} 2>&1
        mv "{params.out_dir}/assembly.fasta" {output.assembly}
        """

rule bwa_align_shorts:
    input:
        unpack(decide_contamination_short),
        assembly = rules.flye_long_based_assembly.output.assembly,
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
    
# TODO Polypolish filter by insert size, recommended by manual: https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish
rule polypolish_filter_by_insert_size:
    input:
        aligned1 = rules.bwa_align_shorts.output.aligned1,
        aligned2 = rules.bwa_align_shorts.output.aligned2,
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

rule polypolish_short_based_polish:
    input:
        r1 = rules.polypolish_filter_by_insert_size.output.filtered1,
        r2 = rules.polypolish_filter_by_insert_size.output.filtered2,
        primary_assembly = rules.flye_long_based_assembly.output.assembly,

    output:
        contig = "results/assembly/{sample}/secondary_assembly.fasta",
    log:
        "logs/polypolish/polish/{sample}.log"
    threads: 1
    conda:
        "../envs/polypolish.yaml"
    shell:
        "polypolish polish {input.primary_assembly} {input.r1} {input.r2} > {output.contig} 2>{log}"

rule get_secondary_assembly:
    input:
        expand("results/assembly/{sample}/secondary_assembly.fasta", sample=SAMPLES.index),

rule quast:
    input:
        contig=rules.polypolish_short_based_polish.output.contig,
    output:
        html = "results/assembly/{sample}/contiguity/report.html",
        out_dir = directory("results/assembly/{sample}/contiguity"),
    log:
        "logs/quast/{sample}.log"
    threads: 16
    conda:
        "../envs/polypolish.yaml"
    shell:
        # TODO add --sam option
        "quast.py {input.contig} -o {output.out_dir} --threads {threads} > {log} 2>&1"

rule merge_samples:
    input: 
        left = expand("results/trimmed/{sample}.1.fastq.gz", sample=SAMPLES.index),
        right = expand("results/trimmed/{sample}.2.fastq.gz", sample=SAMPLES.index),
    output:
        left= "results/merged/left.fastq.gz",
        right=  "results/merged/right.fastq.gz"
    log:
        "logs/merge.log"
    shell:
        """
        exec 2> {log}
        zcat {input.left} | pigz -p 4 -c > {output.left}
        zcat {input.right} | pigz -p 4  -c > {output.right}
        wait
        """

# https://snakemake-wrappers.readthedocs.io/en/v9.10.1/wrappers/bio/trinity.html
rule trinity:
    input:
        left=rules.merge_samples.output.left,
        right=rules.merge_samples.output.right,
    output:
        dir=temp(directory("results/trinity_out_dir/")),
        fas="results/trinity_out_dir/assembly.fasta",
        map="results/trinity_out_dir/gene_trans_map",
    log:
        'logs/trinity/trinity.log',
    params:
        extra="",
    threads: 8
    resources:
        mem_gb=10,
    wrapper:
        "v4.1.0/bio/trinity"



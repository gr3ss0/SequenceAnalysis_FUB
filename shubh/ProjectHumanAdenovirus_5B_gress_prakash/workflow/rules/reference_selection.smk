

rule build_candidate_blastdb: #building the blast database and using for all samples
    input:
        fasta = config["assembly_params"]["candidate_references"]
    output:
        multiext("results/blast_db/candidate_references", ".nhr", ".nin", ".nsq") #not including all of them
    params:
        prefix = "results/blast_db/candidate_references"
    log:
        "logs/blast_db/build.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "makeblastdb -in {input.fasta} -dbtype nucl -out {params.prefix} > {log} 2>&1"


rule blast_contigs: 
    input:
        query = "results/polishing/{sample}/{sample}_polished.fasta",
        db = rules.build_candidate_blastdb.output #making the previous rule a dependency
    output:
        "results/reference_selection/{sample}/blast.tsv" #score for each of the references
    params:
        db_prefix = "results/blast_db/candidate_references",
        evalue = config["assembly_params"]["blast"].get("evalue", "1e-10")
    log:
        "logs/blast/{sample}.log"
    threads: 4
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastn -query {input.query} -db {params.db_prefix} \
            -evalue {params.evalue} -num_threads {threads} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -out {output} > {log} 2>&1
        """


rule select_best_reference: 
    input:
        blast = "results/reference_selection/{sample}/blast.tsv",
        candidates = config["assembly_params"]["candidate_references"] #script needs candidates as input to be able to write the correct one
    output:
        fasta = "results/reference_selection/{sample}/best_reference.fasta", #the one with the best score
        info = "results/reference_selection/{sample}/best_reference.txt"
    log:
        "logs/select_best_reference/{sample}.log"
    conda:
        "../envs/phylo.yaml" 
    script:
        "../scripts/select_best_reference.py" #script new to 5B

def get_map_input_short(wildcards):

    if DECONTAMINATION_ENABLED:

        return {
            "r1": f"results/decontaminated/{wildcards.sample}.1.fastq",
            "r2": f"results/decontaminated/{wildcards.sample}.2.fastq"
        }

    elif not config["analysis_options"].get("skip_trimming", True):

        return {
            "r1": f"results/trimmed/short/{wildcards.sample}.1.fastq.gz",
            "r2": f"results/trimmed/short/{wildcards.sample}.2.fastq.gz"
        }

    else:

        return {
            "r1": SAMPLES.at[wildcards.sample,'short_fq1'],
            "r2": SAMPLES.at[wildcards.sample,'short_fq2']
        }

def get_map_input_long(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return SAMPLES.at[wildcards.sample, 'long_fq']
        
    else:
        # Updated to point to fastp's compressed outputs
        return f"results/trimmed/long/{wildcards.sample}.fastq.gz"

def get_annotation_inputs(wildcards):
    # 1. Get your local samples
    gff_inputs = expand("results/annotation/{sample}/{sample}.gff", sample=SAMPLES.index)
    
    # 2. Append the static external path if it exists in config
    if EXTERNAL_GENOME_ENABLED:
        gff_inputs.append("results/annotation_ext/external_genome/external_genome.gff")
        
    return gff_inputs

def get_decon_input_short(wildcards): #used as input for the decontamination workflow

    if not config["analysis_options"].get("skip_trimming", True):

        return {
            "r1": f"results/trimmed/short/{wildcards.sample}.1.fastq.gz",
            "r2": f"results/trimmed/short/{wildcards.sample}.2.fastq.gz"
        }

    else:

        return {
            "r1": SAMPLES.at[wildcards.sample,'short_fq1'],
            "r2": SAMPLES.at[wildcards.sample,'short_fq2']
        }

def get_final_assembly(wildcards):
    """
    Determine the final assembly file path based on whether short reads are integrated and available.
    """
    if INTEGRATE_SHORT_READS 
    and SAMPLES.at[wildcards.sample, "short_fq1"].notnull() 
    and SAMPLES.at[wildcards.sample, "short_fq2"].notnull():
        return f"results/assembly/{wildcards.sample}/secondary_assembly.fasta"
    else:
        return f"results/assembly/{wildcards.sample}/primary_assembly.fasta"
def get_map_input_short(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return {
            "r1": SAMPLES.at[wildcards.sample, 'short_fq1'], 
            "r2": SAMPLES.at[wildcards.sample, 'short_fq2']
        }
    else:
        # Updated to point to fastp's compressed outputs
        return {
            "r1": f"results/trimmed/short/{wildcards.sample}.1.fastq.gz", 
            "r2": f"results/trimmed/short/{wildcards.sample}.2.fastq.gz"
        }

def decide_contamination_short(wildcards):
    if config["analysis_options"].get("skip_kraken_screening", True):
        return decide_trimming_short
    else:
        return{
            "r1": f"results/decontaminated/short/{wildcards.sample}.1.fastq",
            "r2": f"results/decontaminated/short/{wildcards.sample}.2.fastq"
        }

def decide_trimming_short(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return {
            "r1": SAMPLES.at[wildcards.sample, 'short_fq1'], 
            "r2": SAMPLES.at[wildcards.sample, 'short_fq2']
        }
    else:
        # Updated to point to fastp's compressed outputs
        return {
            "r1": f"results/trimmed/short/{wildcards.sample}.1.fastq.gz", 
            "r2": f"results/trimmed/short/{wildcards.sample}.2.fastq.gz"
        }

def decide_contamination_long(wildcards):
    if config["analysis_options"].get("skip_kraken_screening", True):
        return decide_trimming_long
    else:
        return f"results/decontaminated/long/{wildcards.sample}.fastq"
            

def decide_trimming_long(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return SAMPLES.at[wildcards.sample, 'long_fq']
        
    else:
        # Updated to point to fastp's compressed outputs
        return f"results/trimmed/long/{wildcards.sample}.fastq.gz"

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
    if "external_genome" in config and len(config["external_genome"]) > 0:
        gff_inputs.append("results/annotation_ext/external_genome/external_genome.gff")
        
    return gff_inputs
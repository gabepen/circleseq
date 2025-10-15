# Get parameters from command line
import sys

# Parse command line arguments
if len(sys.argv) < 3:
    print("Usage: snakemake --config reference=<reference_name> target_string=<target_string>")
    sys.exit(1)

# Get config values
REFERENCE = config.get("reference", "")
TARGET_STRING = config.get("target_string", "")

if not REFERENCE or not TARGET_STRING:
    print("Error: Both 'reference' and 'target_string' must be provided in config")
    print("Usage: snakemake --config reference=<reference_name> target_string=<target_string>")
    sys.exit(1)

# Automatically discover samples from input directory
SAMPLES, = glob_wildcards("input/{sample}_R1.fq.gz")
SAMPLES = [sample for sample in SAMPLES if TARGET_STRING in sample]

if not SAMPLES:
    print(f"Warning: No samples found containing '{TARGET_STRING}' in input directory")
    print("Available samples:", [sample for sample, in glob_wildcards("input/{sample}_R1.fq.gz")])
    sys.exit(1)

print(f"Found {len(SAMPLES)} samples: {SAMPLES}")
print(f"Using reference: {REFERENCE}")

rule all:
    input:
        expand("{sample}_{reference}.png", sample=SAMPLES, reference=REFERENCE)
rule bwa_map:
    input:
        "input/{sample}_R1.fq.gz",
        "input/{sample}_R2.fq.gz",
        f"references/{REFERENCE}.fa"
    output:
        "{sample}_{reference}_initial_mapping.sam"
    benchmark:
        "{sample}_{reference}_imap_benchmark.tsv"
    threads: 20
    shell: 
        "bwa mem -t {threads} {input[2]} {input[0]} {input[1]} > {output}"
rule process_circle_alignments:
    input:
        "{sample}_{reference}_initial_mapping.sam",
        f"references/{REFERENCE}.fa.fai"
    output:
        "{sample}_{reference}_split.fa",
        "{sample}_{reference}_regions.bed"
    benchmark:
        "{sample}_{reference}_process_benchmark.tsv"
    shell:
        "python3 process_circle_alignments.py -r {input[1]} -f {output[0]} -b {output[1]} -d 10 {input[0]}"
rule makefasta:
    input:
        "{sample}_{reference}_regions.bed",
        f"references/{REFERENCE}.fa"
    output:
        "{sample}_{reference}_reference.fa"
    shell:
        "bedtools getfasta -fi {input[1]} -bed {input[0]} -fo {output} -name"
rule make_brpileup:
    input:
        "{sample}_{reference}_split.fa",
        "{sample}_{reference}_regions.bed",
        "{sample}_{reference}_reference.fa"
    output:
        "{sample}_{reference}_consensus.sam"
    benchmark:
        "{sample}_{reference}_refpileup_benchmark.tsv"
    threads: 20
    shell:
        "python3 make_brpileup_noconsensus.py -t {threads} -c {output} {input[1]} {input[2]} {input[0]}"
rule index:
    input:
        f"references/{REFERENCE}.fa"
    output:
        f"references/{REFERENCE}.fa.fai"
    shell:
        f"samtools faidx references/{REFERENCE}.fa"
rule samprocess:
    input:
        "{sample}_{reference}_consensus.sam",
        f"references/{REFERENCE}.fa.fai"
    output:
        "{sample}_{reference}.sorted.bam"
    shell:
        "samtools view -ht {input[1]} {input[0]} | samtools view -Sb | samtools sort > {output}"
rule mpileup:
    input:
        "{sample}_{reference}.sorted.bam"
    output:
        "{sample}_{reference}_variants.txt"
    shell:
        f"samtools mpileup -B -Q 17 -f references/{REFERENCE}.fa {{input}} > {{output}}"
rule count_mutations:
    input:
        "{sample}_{reference}_variants.txt"
    output:
        "{sample}_{reference}.txt"
    shell:
        "python3 count_mutations.py -t 2 -m {input[0]} > {output}"
rule graph_mutations:
    input:
        "{sample}_{reference}.txt"
    output:
        "{sample}_{reference}.png"
    shell:
        "python3 graph_mutations.py {input}"
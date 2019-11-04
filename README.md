# circleseq
 Primary analysis pipeline for ultra-accurate sequencing data

## Dependencies
### External Software
bwa v0.7.17
bedtools v2.26.0
samtools v1.7
snakemake v5.7.1
*optional:* conda 4.7.12
### Python Packages
scikit-bio v0.5.5
biopython v1.74

## Installation
Copy the git directory:

git clone github.com/jmcbroome/circleseq

Ensure external software dependences are installed and on your shell's path.

Package dependencies can be installed independently or the circleseq.yml environment may be used via conda.

## Formatting Files
The snakefile as it stands expects input files in the format of {sample}\_R1.fa and {sample}\_R2.fa under the "input" file folder. 
Reference data is expected under references/{reference_genome}.fa, replacing bracketed values with the specific values of your sample and the reference genome name.

The file structure should look like this:

    Directory with scripts
        input
             {sample}_R1.fa
             {sample}_R2.fa
        references
             {reference_genome}.fa

## Usage
First, ensure your reference of choice is indexed with bwa.

bwa index references/{reference_genome}

Then simply call:

snakemake -j {max_threads} {sample}\_{reference_genome}\_errors.txt

Again, replacing bracketed values with the name of your sample, the name of your reference genome file, and with the maximum threads value being the maximum number of threads available to the pipeline for processing. Default value for max_threads is 1.
Add the argument "--use-conda circleseq.yml" as an alternative to global installation of requisite packages, or activate the environment with conda and call snakemake from within it.

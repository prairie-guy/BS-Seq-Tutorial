from pathlib import Path

from fname import fname, expand_fnames

# Wildcard to access samples (Change to use config.yaml)
# "," indicates tuple umpacking
SAMPLES, = glob_wildcards("samples/{sample}.fastq")

# Paths and Parameters
GENOME_DIR  = Path("../Genome_reference")
GENOME_FILE = GENOME_DIR/"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOME_URL  = "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
CPUS = 32
CPUS_MIN = 4

# Rule all defines the final outputs of the workflow

rule all:
    input:
        #{GENOME_FILE},
        #{GENOME_DIR/"Bisulfite_Genome/CT_conversion/BS_CT.1.bt2"},
        expand_fnames("fastqc",     "zip",   "fastqc", SAMPLES),
        #expand_fnames("trimmed",    "fq.gz", "trimmed", SAMPLES),
        expand_fnames("fastqc_2",   "zip",   "trimmed_fastqc", SAMPLES),
        #expand_fnames("bismark",    "bam",   "trimmed_bismark_bt2", SAMPLES),
        #expand_fnames("methylation","txt",   "trimmed_bismark_bt2_splitting_report", SAMPLES)


# Download and repare the genome references
rule gemome_refs:
    output: GENOME_FILE,
            GENOME_DIR/"Bisulfite_Genome/CT_conversion/BS_CT.1.bt2"
    log: "log/genome_refs"
    shell:
        """
        wget -O {GENOME_FILE}.gz {GENOME_URL} 2> {log}
        gunzip {GENOME_FILE}.gz 2> {log}
        bismark_genome_preparation --parallel {CPUS} --verbose {GENOME_DIR} 2> {log}
        """
rule fastqc:
    input:  fname("samples", "fastq")
    output: fname("fastqc",  "zip", "fastqc")
    params:  output_dir = "fastqc"
    log:    fname("log", "log", "fastq")
    shell:  "fastqc {input} -o {params.output_dir} 2> {log}"

# Tyring to just trim ends without removing adapter
rule cutadapt:
    input:  fname("samples", "fastq")
    output: fname("trimmed", "fq.gz", "trimmed")
    params:  output_dir = "trimmed"
    log:    fname("log", "log", "trimmed")
    shell: "cutadapt -u 5 -u -24  {input} | gzip > {output} "


# rule trim_galore:
#     input:  fname("samples", "fastq")
#     output: fname("trimmed", "fq.gz", "trimmed")
#     params:  output_dir = "trimmed"
#     log:    fname("log", "log", "trimmed")
#     shell:  "trim_galore  --gzip --output_dir {params.output_dir}  {input} 2> {log}"

rule fastqc_2:
    input:  fname("trimmed",  "fq.gz", "trimmed")
    output: fname("fastqc_2", "zip",   "trimmed_fastqc")
    params: output_dir = "fastqc_2"
    log:    fname("log", "log", "trimmed_fastqc")
    shell:  "fastqc {input} -o {params.output_dir} 2> {log}"


rule bismark_alignment:
    input:  fname("trimmed", "fq.gz", "trimmed")
    output: fname("bismark", "bam",   "trimmed_bismark_bt2")
    params: output_dir = "bismark"
    log: fname("log", "log", "trimmed_bismark_bt2")
    shell: "bismark --genome {GENOME_DIR} --single_end {input} --multicore  {CPUS_MIN} --output_dir {params.output_dir} 2> {log}"


rule methylation_extraction:
    input:  fname("bismark",     "bam", "trimmed_bismark_bt2")
    output: fname("methylation", "txt", "trimmed_bismark_bt2_splitting_report")
    params: output_dir = "methylation"
    log:    fname("log","log","methylation")
    shell: "bismark_methylation_extractor --single-end --bedGraph  --multicore {CPUS} --output_dir {params.output_dir} {input}"

           # Quality control with FastQC
# rule fastqc:
#     input:  "samples/{sample}.fastq"
#     output: "fastqc/{sample}_fastqc.zip"
#     params:  output_dir = "fastqc"
#     log:    "log/fastqc.{sample}.log"
#     shell:  "fastqc {input} -o {params.output_dir} 2> {log}"

# rule trim_galore:
#     input:  "samples/{sample}.fastq"
#     output: "trimmed/{sample}_trimmed.fq.gz"
#     params:  output_dir = "trimmed"
#     log:    "log/trim_galore.{sample}.log"
#     shell:  "trim_galore  --gzip --output_dir {params.output_dir}  {input} 2> {log}"

# # Post-Trimming Quality control with FastQC
# rule post_trimmed_fastqc:
#     input:  "trimmed/{sample}_trimmed.fq.gz"
#     output: "post-trimmed_fastqc/{sample}_trimmed_fastqc.zip"
#     params:  output_dir = "post-trimmed_fastqc"
#     log:    "log/post-trimmed_fastqc.{sample}.log"
#     shell:  "fastqc {input} -o {params.output_dir} 2> {log}"

# # Aligning with Bismark
# rule bismark_alignment:
#     input: "trimmed/{sample}_trimmed.fq.gz"
#     output: "bismark/{sample}_trimmed_bismark_bt2.bam"
#     params: output_dir = "bismark"
#     log: "log/bismark_alignment.{sample}.log"
#     shell: "bismark --genome {GENOME_DIR} --single_end {input} --multicore  {CPUS_MIN} --output_dir {params.output_dir} 2> {log}"

# # Methylation Extraction
# rule methylation_extraction:
#     input: "bismark/{sample}_trimmed_bismark_bt2.bam"
#     output: "methylation/{sample}_trimmed_bismark_bt2_splitting_report.txt"
#     params: output_dir = "methylation"
#     shell: "bismark_methylation_extractor --single-end --bedGraph  --multicore {CPUS} --output_dir {params.output_dir} {input}"

#!~/miniforge3/bin/python3

import json

configfile: "snakemake_pipeline/config_file.yaml"

FILES           = json.load(open(config['SAMPLES_JSON']))
SAMPLES         = sorted(FILES.keys())
WORKING_FOLDER  = config['WORKING_FOLDER']
FASTQ_FOLDER    = config['FASTQ_FOLDER']
BWA_INDEX       = config['BWA_INDEX']
NORM_FOLDER     = config['NORM_FOLDER']
BWA_FOLDER      = config['BWA_FOLDER']
PYTHON_FOLDER   = config['PYTHON_FOLDER']

TARGETS = []
post  = expand("{WORKING_FOLDER}/{sample}_nanoseq_post.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
cov = expand("{WORKING_FOLDER}/{sample}_nanoseq_cov.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
part = expand("{WORKING_FOLDER}/{sample}_nanoseq_part.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
indels = expand("{WORKING_FOLDER}/{sample}_nanoseq_indel.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
var = expand("{WORKING_FOLDER}/{sample}_nanoseq_var.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
dsa = expand("{WORKING_FOLDER}/{sample}_nanoseq_dsa.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
TARGETS.extend(post)
TARGETS.extend(cov)
TARGETS.extend(part)
TARGETS.extend(var)
TARGETS.extend(indels)
TARGETS.extend(dsa)

localrules: all

rule all:
    input: TARGETS

rule extract_tags:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: "{WORKING_FOLDER}/data/trimmed/{sample}_1.fq.gz", "{WORKING_FOLDER}/data/trimmed/{sample}_2.fq.gz"
    threads: 1
    message: "Extract UMIs from raw nanoseq reads"
    log: "{WORKING_FOLDER}/00_log/{sample}_stdOut.extractTags","{WORKING_FOLDER}/00_log/{sample}_error.extractTags"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {PYTHON_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "python python/extract_tags.py -a {input[0]} -b {input[1]} -c {output[0]} -d {output[1]} -m 3 -s 2 -l 150 1> {log[0]} 2> {log[1]}"
"""

rule align_tumor:
    input: "{WORKING_FOLDER}/data/trimmed/{sample}_1.fq.gz", "{WORKING_FOLDER}/data/trimmed/{sample}_2.fq.gz"
    output: "{WORKING_FOLDER}/data/align/{sample}.sam"
    threads: 16
    resources:
        mem_mb=40000,
        slurm_partition="fn_long"
    message: "Align trimmed reads in {input} to the genome using {threads} threads"
    log: "{WORKING_FOLDER}/00_log/{sample}.bwa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "bwa mem \
    -t {threads} -C {BWA_INDEX} {input[0]} {input[1]} > {output}"
"""

rule align_normal:
    input:
        norm_r1 = lambda wildcards: FILES[wildcards.sample]['matched_R1'],
        norm_r2 = lambda wildcards: FILES[wildcards.sample]['matched_R2']
    output: "{WORKING_FOLDER}/data/align_norm/{sample}.bam"
    threads: 16
    resources:
        mem_mb=40000,
        slurm_partition="fn_long"
    message: "Align trimmed reads in normal samples: {input} to the genome using {threads} threads"
    log: "{WORKING_FOLDER}/00_log/{sample}_norm.bwa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "bwa mem -t {threads} {BWA_INDEX} {input[0]} {input[1]} 2> {log} | samtools sort -@ {threads} -o {output} -"
"""

rule add_rc_mc_tags:
    input: "{WORKING_FOLDER}/data/align/{sample}.sam"
    output: "{WORKING_FOLDER}/data/align/{sample}_od.bam"
    threads: 8
    message: "append rc & mc tags"
    log: "{WORKING_FOLDER}/00_log/{sample}.bamsormadup"
    shell: """
    mkdir -p /tmp/arinze
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B /tmp/arinze/ \
    --env "PATH=/home/ubuntu/miniforge3/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "bamsormadup tmpfile=/tmp/arinze/ inputformat=sam rcsupport=1 threads=8 < {input} > {output}"
    """

rule append_rb_tag_filter:
    input:  "{WORKING_FOLDER}/data/align/{sample}_od.bam"
    output: "{WORKING_FOLDER}/data/filtered/{sample}.bam"
    threads: 8
    log: "{WORKING_FOLDER}/00_log/{sample}.bamaddreadbundles"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "bamaddreadbundles  -I {input} -O {output} 2> {log}"
    """

rule index_tumor_bam:
    input:  "{WORKING_FOLDER}/data/filtered/{sample}.bam"
    output: "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    threads: 8
    log: "{WORKING_FOLDER}/00_log/{sample}.samtoolsTumIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "samtools index -@ 4 {input} -o {output} 2> {log}"
    """

rule index_normal_bam:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam"
    output: "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai"
    threads: 8
    log: "{WORKING_FOLDER}/00_log/{sample}.samtoolsNormIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "samtools index -@ 4 {input} -o {output} 2> {log}"
    """

rule nanoseq_cov:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_cov.txt")
    threads: 10
    message: "Nanoseq_coverage"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_cov"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {PYTHON_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "python {PYTHON_FOLDER}/runNanoSeq.py -t 10 -A {input[0]} -B {input[1]} -R {BWA_INDEX} cov -Q 0 --exclude "MT,GL%,NC_%,hs37d5" 2>{log}"
    """

rule nanoseq_partition:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_part.txt")
    threads: 10
    message: "Nanoseq_partition"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_partition"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {PYTHON_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "python {PYTHON_FOLDER}/runNanoSeq.py -t 1 -A {input[0]} -B {input[1]} -R {BWA_INDEX} part -n 10 2>{log}"
    """

rule nanoseq_dsa:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_dsa.txt")
    threads: 10
    message: "Nanoseq_dsa"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_dsa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {PYTHON_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "python {PYTHON_FOLDER}/runNanoSeq.py -t 10 -A {input[0]} -B {input[1]} -R {BWA_INDEX} dsa -C SNP.sorted.bed.gz -D NOISE.sorted.bed.gz -d 2 -q 30 2>{log}"
    """

rule nanoseq_var:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_var.txt")
    threads: 10
    message: "Nanoseq_var"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_var"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {PYTHON_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "python {PYTHON_FOLDER}/runNanoSeq.py -t 60 -A {input[0]} -B {input[1]} -R {BWA_INDEX} var \
    -a 50 -b 5 -c 0 -f 0.9 -i 1 -m 8 -n 3 -p 0 -q 60 -r 144 -v 0.01 -x 8 -z 12"
    """

rule nanoseq_indel:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_indel.txt")
    threads: 10
    message: "Nanoseq_indel"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_indel"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {PYTHON_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "python {PYTHON_FOLDER}/runNanoSeq.py -t 60 -A {input[0]} -B {input[1]} -R {BWA_INDEX} indel -s sample \
  --rb 2 --t3 135 --t5 10 --mc 16 2>{log}"
    """

rule nanoseq_post:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam", 
            "{WORKING_FOLDER}/{sample}_nanoseq_indel.txt", "{WORKING_FOLDER}/{sample}_nanoseq_var.txt",
            "{WORKING_FOLDER}/{sample}_nanoseq_cov.txt", "{WORKING_FOLDER}/{sample}_nanoseq_dsa.txt",
            "{WORKING_FOLDER}/{sample}_nanoseq_part.txt"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_post.txt")
    threads: 2
    message: "Nanoseq_post"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_post"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} -B {PYTHON_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.0.2.sif /bin/bash -c "python {PYTHON_FOLDER}/runNanoSeq.py -t 2 -A {input[0]} -B {input[1]} -R {BWA_INDEX} post 2>{log}"
    """
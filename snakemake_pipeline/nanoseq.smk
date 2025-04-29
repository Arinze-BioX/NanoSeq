#!~/miniforge3/bin/python3

import json

configfile: "snakemake_pipeline/config_file.yaml"

FILES           = json.load(open(config['SAMPLES_JSON']))
SAMPLES         = sorted(FILES.keys())
WORKING_FOLDER  = config['WORKING_FOLDER']
FASTQ_FOLDER    = config['FASTQ_FOLDER']
BWA_INDEX       = config['BWA_INDEX']
BWA_FOLDER      = config['BWA_FOLDER']
CODE_HOME       = config['CODE_HOME']
NOISE           = config['NOISE']
partitions      = config['partitions']
GERMLINE_TBL   = config['GERMLINE_TBL']
GERMLINE_VCF   = config['GERMLINE_VCF']
CELL_LINE_GENOTYPE_REGIONS = config['CELL_LINE_GENOTYPE_REGIONS']
VERIFYBAMID_HOME = config['VERIFYBAMID_HOME']

TARGETS = []
post  = expand("{WORKING_FOLDER}/{sample}_nanoseq_post.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
eff = expand("{WORKING_FOLDER}/{sample}_nanoseq_efficiency.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
variant_qc = expand("{WORKING_FOLDER}/{sample}_variant_qc.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
verify = expand("{WORKING_FOLDER}/{sample}_nanoseq_efficiency.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER) 
verify_cellLine = expand("{WORKING_FOLDER}/data/cell_qc/{sample}_cell_verify.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
TARGETS.extend(post)
TARGETS.extend(verify)
TARGETS.extend(variant_qc)
TARGETS.extend(eff)
TARGETS.extend(verify_cellLine)

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
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} \
    --env "PATH=/home/ubuntu/miniforge3/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/extract_tags.py -a {input[0]} -b {input[1]} -c {output[0]} -d {output[1]} -m 3 -s 2 -l 150 1> {log[0]} 2> {log[1]}"
"""

rule align_tumor:
    input: "{WORKING_FOLDER}/data/trimmed/{sample}_1.fq.gz", "{WORKING_FOLDER}/data/trimmed/{sample}_2.fq.gz"
    output: "{WORKING_FOLDER}/data/align/{sample}.bam"
    threads: 16
    resources:
        mem_mb=40000,
        slurm_partition="petljaklab"
    message: "Align trimmed reads in {input} to the genome using {threads} threads"
    log: "{WORKING_FOLDER}/00_log/{sample}.bwa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bwa mem \
    -t {threads} -C {BWA_INDEX} {input[0]} {input[1]} | \
  samtools sort -@ 4 -o {output} -"
"""

rule align_normal:
    input:
        norm_r1 = lambda wildcards: FILES[wildcards.sample]['matched_R1'],
        norm_r2 = lambda wildcards: FILES[wildcards.sample]['matched_R2']
    output: "{WORKING_FOLDER}/data/align_norm/{sample}.bam"
    threads: 16
    resources:
        mem_mb=40000,
        slurm_partition="petljaklab"
    message: "Align trimmed reads in normal samples: {input} to the genome using {threads} threads"
    log: "{WORKING_FOLDER}/00_log/{sample}_norm.bwa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bwa mem -t {threads} {BWA_INDEX} {input[0]} {input[1]} 2> {log} | samtools sort -@ {threads} -o {output} -"
"""

rule index_normal_bam:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam"
    output: "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai"
    threads: 8
    log: "{WORKING_FOLDER}/00_log/{sample}.samtoolsNormIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {NORM_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "samtools index -@ 4 {input} -o {output} 2> {log}"
    """

rule add_rc_mc_tags:
    input: "{WORKING_FOLDER}/data/align/{sample}.bam"
    output: "{WORKING_FOLDER}/data/align/{sample}_od.bam"
    threads: 16
    message: "append rc & mc tags"
    log: "{WORKING_FOLDER}/00_log/{sample}.bamsormadup"
    shell: """
    mkdir -p {WORKING_FOLDER}/tmp_{wildcards.sample}_rcmc
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bamsormadup tmpfile={WORKING_FOLDER}/tmp_{wildcards.sample}_rcmc inputformat=bam rcsupport=1 threads={threads} < {input} > {output}"
    rm -rf {WORKING_FOLDER}/tmp_{wildcards.sample}_rcmc
    """

rule mark_duplicates:
    input:  "{WORKING_FOLDER}/data/align/{sample}_od.bam"
    output: "{WORKING_FOLDER}/data/filtered/{sample}_dupMarked.bam"
    threads: 1
    log: "{WORKING_FOLDER}/00_log/{sample}.markdups"
    shell: """
    mkdir -p {WORKING_FOLDER}/tmp_{wildcards.sample}_dups
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bammarkduplicatesopt I={input} O={output} M={log} tmpfile={WORKING_FOLDER}/tmp_{wildcards.sample}_dups"
    rm -rf {WORKING_FOLDER}/tmp_{wildcards.sample}_dups
    """

rule append_rb_tag_filter:
    input:  "{WORKING_FOLDER}/data/filtered/{sample}_dupMarked.bam"
    output: "{WORKING_FOLDER}/data/filtered/{sample}.bam"
    threads: 1
    log: "{WORKING_FOLDER}/00_log/{sample}.bamaddreadbundles"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bamaddreadbundles  -I {input} -O {output} 2> {log}"
    """

rule index_tumor_bam:
    input:  "{WORKING_FOLDER}/data/filtered/{sample}.bam"
    output: "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    threads: 4
    log: "{WORKING_FOLDER}/00_log/{sample}.samtoolsTumIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "samtools index -@ 4 {input} -o {output} 2> {log}"
    """


rule deduplicate_bam:
    input:  "{WORKING_FOLDER}/data/filtered/{sample}.bam"
    output: "{WORKING_FOLDER}/data/dedup/{sample}.bam"
    threads: 1
    message: "Carrying out dedulication"
    log: "{WORKING_FOLDER}/00_log/{sample}.randomreadinbundle"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "randomreadinbundle -I {input} -O {output} 2>{log}"
    """

rule index_dedup_bam:
    input:  "{WORKING_FOLDER}/data/dedup/{sample}.bam"
    output: "{WORKING_FOLDER}/data/dedup/{sample}.bam.bai"
    threads: 4
    log: "{WORKING_FOLDER}/00_log/{sample}.samtoolsDedupIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "samtools index -@ 4 {input} -o {output} 2> {log}"
    """

rule analyse_efficiency:
    input:  "{WORKING_FOLDER}/data/dedup/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam", "{WORKING_FOLDER}/data/dedup/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_efficiency.txt")
    threads: 12
    message: "Calculating efficiency metrics"
    log: "{WORKING_FOLDER}/00_log/{sample}.efficiency_pl"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "{CODE_HOME}/perl/efficiency_nanoseq.pl -t {threads} -dedup {input[0]} -duplex {input[1]} -o {wildcards.sample} -r {BWA_INDEX} 2>{log}"
    """


rule nanoseq_cov:
    input:  "{WORKING_FOLDER}/data/dedup/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam",
            "{WORKING_FOLDER}/data/dedup/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_cov.txt")
    threads: 16
    message: "Nanoseq_coverage"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_cov"
    shell: """
    mkdir {WORKING_FOLDER}/{wildcards.sample}
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input[0]} -B {input[1]} -R {BWA_INDEX} cov -Q 0 --exclude "MT,GL%%,NC_%,hs37d5" 2>{log}"
    """

rule nanoseq_partition:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam","{WORKING_FOLDER}/{sample}_nanoseq_cov.txt",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_part.txt")
    threads: 1
    message: "Nanoseq_partition"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_partition"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t 1 -A {input[0]} -B {input[1]} -R {BWA_INDEX} part -n {partitions} 2>{log}"
    """
rule nanoseq_dsa:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam", "{WORKING_FOLDER}/{sample}_nanoseq_part.txt",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_dsa.txt")
    threads: partitions
    message: "Nanoseq_dsa"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_dsa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input[0]} -B {input[1]} -R {BWA_INDEX} dsa -C SNP.sorted.bed.gz -D {NOISE} -d 2 -q 30 2>{log}"
    """

rule nanoseq_var:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam", "{WORKING_FOLDER}/{sample}_nanoseq_dsa.txt",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_var.txt")
    threads: partitions
    message: "Nanoseq_var"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_var"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input[0]} -B {input[1]} -R {BWA_INDEX} var \
    -a 50 -b 5 -c 0 -f 0.9 -i 1 -m 8 -n 3 -p 0 -q 60 -r 144 -v 0.01 -x 8 -z 12"
    """

rule nanoseq_indel:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam", "{WORKING_FOLDER}/data/filtered/{sample}.bam", "{WORKING_FOLDER}/{sample}_nanoseq_dsa.txt",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai", "{WORKING_FOLDER}/data/filtered/{sample}.bam.bai"
    output: touch("{WORKING_FOLDER}/{sample}_nanoseq_indel.txt")
    threads: partitions
    message: "Nanoseq_indel"
    log: "{WORKING_FOLDER}/00_log/{sample}.nanoseq_indel"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input[0]} -B {input[1]} -R {BWA_INDEX} indel -s {wildcards.sample} \
  --rb 2 --t3 140 --t5 5 2>{log}"
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
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t 2 -A {input[0]} -B {input[1]} -R {BWA_INDEX} post 2>{log}"
    """

rule nanoseq_plot_variant_qc:
    input:  "{WORKING_FOLDER}/{sample}_nanoseq_post.txt"
    output: touch("{WORKING_FOLDER}/{sample}_variant_qc.txt")
    threads: 1
    message: "Plot graphs to assess quality of variant calling and filtering"
    log: "{WORKING_FOLDER}/00_log/{sample}.varaint_qc"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "Rscript {CODE_HOME}/R/variant_filtering_qc.R {WORKING_FOLDER}/{wildcards.sample}/tmpNanoSeq/post/discardedvariants.csv \
    {WORKING_FOLDER}/{wildcards.sample}/tmpNanoSeq/post/variants.csv {WORKING_FOLDER}/{wildcards.sample}/tmpNanoSeq partitions 2>{log}"
    """

rule verify_cell_line_origin:
    input:  "{WORKING_FOLDER}/data/align_norm/{sample}.bam",
            "{WORKING_FOLDER}/data/align_norm/{sample}.bam.bai"
    output: "{WORKING_FOLDER}/data/cell_qc/{sample}_cell_verify.txt"
    threads: 1
    message: "Verify cell line origin"
    log: "{WORKING_FOLDER}/00_log/{sample}.verify_ccell_line"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {CELL_LINE_GENOTYPE_REGIONS} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bcftools mpileup \
         -R {CELL_LINE_GENOTYPE_REGIONS} \
         --fasta-ref {BWA_INDEX} \
         -A {input[0]} | bcftools call -c > {WORKING_FOLDER}/{GERMLINE_VCF};
gatk VariantsToTable -V {WORKING_FOLDER}/data/cell_qc/{wildcards.sample}_{GERMLINE_VCF} \
-O {WORKING_FOLDER}/data/cell_qc/{wildcards.sample}_{GERMLINE_TBL};
cat {WORKING_FOLDER}/data/cell_qc/{wildcards.sample}_{GERMLINE_TBL} | Rscript \
./R/qc_identity.R > {ouput} 2>{log}"
    """


rule verify_cell_origin:
    input:  "{WORKING_FOLDER}/data/dedup/{sample}.bam"
    output: touch("{WORKING_FOLDER}/{sample}_verify.txt")
    threads: 1
    message: "Verify bam ID"
    log: "{WORKING_FOLDER}/00_log/{sample}.verifyBamId"
    shell: """
    singularity exec -B {BWA_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {VERIFYBAMID_HOME} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c \
    "VerifyBamID --SVDPrefix {VERIFYBAMID_HOME}/resource/1000g.100k.b38.vcf.gz.dat \
    --Reference {BWA_INDEX} --BamFile {input} 2>{log}"
    """
#!~/miniforge3/bin/python3

import json
import os # Import the os module

configfile: "snakemake_pipeline/config_file.yaml"

NORMALS         = json.load(open(config['NORMALS_JSON']))
FILES           = json.load(open(config['SAMPLES_JSON']))
SAMPLES         = sorted(FILES.keys())
WORKING_FOLDER  = config['WORKING_FOLDER']
FASTQ_FOLDER    = config['FASTQ_FOLDER']
BWA_INDEX       = config['BWA_INDEX']
BWA_FOLDER      = config['BWA_FOLDER']
CODE_HOME       = config['CODE_HOME']
NOISE           = config['NOISE']
partitions      = config['partitions']
GERMLINE_TBL    = config['GERMLINE_TBL']
GERMLINE_VCF    = config['GERMLINE_VCF']
CELL_LINE_GENOTYPE_REGIONS = config['CELL_LINE_GENOTYPE_REGIONS']
VERIFYBAMID_HOME    = config['VERIFYBAMID_HOME']

TARGETS = []
post  = expand("{WORKING_FOLDER}/{sample}_nanoseq_post.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
eff = expand("{WORKING_FOLDER}/{sample}_nanoseq_efficiency.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
variant_qc = expand("{WORKING_FOLDER}/{sample}_variant_qc.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
verify = expand("{WORKING_FOLDER}/{sample}_nanoseq_verify.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER) # This seems duplicated with eff?
# verify_cellLine = expand("{WORKING_FOLDER}/data/cell_qc/{sample}_cell_verify.txt", sample = SAMPLES, WORKING_FOLDER = WORKING_FOLDER)
TARGETS.extend(post)
TARGETS.extend(verify)
TARGETS.extend(variant_qc)
TARGETS.extend(eff)
# TARGETS.extend(verify_cellLine)

localrules: all

rule all:
    input: TARGETS

rule extract_tags:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample][wildcards.lane]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample][wildcards.lane]['R2']
    output:
        # Outputting to 'trimmed' directory
        r1_trimmed=f"{WORKING_FOLDER}/data/trimmed/{{sample}}-{{lane}}_1.fq.gz",
        r2_trimmed=f"{WORKING_FOLDER}/data/trimmed/{{sample}}-{{lane}}_2.fq.gz"
    threads: 1
    message: "Extract UMIs from raw nanoseq reads"
    log:
        stdout=f"{WORKING_FOLDER}/00_log/{{sample}}-{{lane}}_stdOut.extractTags",
        stderr=f"{WORKING_FOLDER}/00_log/{{sample}}-{{lane}}_error.extractTags"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" \
    nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/extract_tags.py -a {input.r1} \
    -b {input.r2} -c {output.r1_trimmed} -d {output.r2_trimmed} -m 3 -s 2 -l 151 1> {log.stdout} 2> {log.stderr}"
"""

rule align_tumor:
    input:
        # Correctly taking input from 'trimmed' directory
        r1=f"{WORKING_FOLDER}/data/trimmed/{{sample}}-{{lane}}_1.fq.gz",
        r2=f"{WORKING_FOLDER}/data/trimmed/{{sample}}-{{lane}}_2.fq.gz"
    output:
        # Outputting the initial alignment BAM
        bam=f"{WORKING_FOLDER}/data/align/{{sample}}-{{lane}}_aligned.bam" # Changed name slightly to avoid overwrite confusion
    threads: 16
    resources:
        mem_mb=40000,
        slurm_partition="petljaklab"
    message: "Align trimmed reads in {input} to the genome using {threads} threads"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}-{{lane}}.bwa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bwa mem \
    -R '@RG\\tID:group1\\tSM:{{sample}}\\tLB:{{lane}}\\tPL:ILLUMINA\\tPU:flowcellA.{{lane}}' \
    -t {threads} -C {BWA_INDEX} {input.r1} {input.r2} | \
    samtools view -b -o {output.bam} -"
"""
# Note: Added {{lane}} to PU field in RG string for better distinction if needed later

rule align_normal:
    input:
        norm_r1 = lambda wildcards: NORMALS[wildcards.sample]['R1'],
        norm_r2 = lambda wildcards: NORMALS[wildcards.sample]['R2']
    output:
        bam=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam"
    threads: 16
    resources:
        mem_mb=64000,
        slurm_partition="petljaklab"
    message: "Align trimmed reads in normal samples: {input} to the genome using {threads} threads"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}_norm.bwa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bwa mem \
    -R '@RG\\tID:groupN\\tSM:{{sample}}_N\\tLB:norm_lib\\tPL:ILLUMINA\\tPU:flowcellN.laneN' \
    -t {threads} {BWA_INDEX} {input.norm_r1} {input.norm_r2} 2> {log} | \
    samtools sort -@ {threads} -o {output.bam} -"
"""

rule index_normal_bam:
    input:
        bam=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam"
    output:
        bai=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam.bai"
    threads: 8
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.samtoolsNormIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "samtools index -@ {threads} {input.bam} -o {output.bai} 2> {log}"
    """

rule add_rc_mc_tags:
    input:
        bam=f"{WORKING_FOLDER}/data/align/{{sample}}-{{lane}}_aligned.bam"
    output:
        bam_od=f"{WORKING_FOLDER}/data/align/{{sample}}-{{lane}}_od.bam"
    threads: 16
    resources:
        mem_mb=50000,
        slurm_partition="petljaklab"
    message: "append rc & mc tags"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}-{{lane}}.bamsormadup"
    shell: """
    mkdir -p {WORKING_FOLDER}/tmp_{wildcards.sample}_{wildcards.lane}_rcmc
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bamsormadup tmpfile={WORKING_FOLDER}/tmp_{wildcards.sample}_{wildcards.lane}_rcmc inputformat=bam rcsupport=1 threads={threads} < {input.bam} > {output.bam_od}"
    rm -rf {WORKING_FOLDER}/tmp_{wildcards.sample}_{wildcards.lane}_rcmc
    """

rule mark_duplicates:
    input:
        bam_od=f"{WORKING_FOLDER}/data/align/{{sample}}-{{lane}}_od.bam"
    output:
        bam_dupmarked=f"{WORKING_FOLDER}/data/align/{{sample}}-{{lane}}_dupMarked.bam"
    threads: 1 # bammarkduplicatesopt is often single-threaded
    log:
        metrics=f"{WORKING_FOLDER}/00_log/{{sample}}-{{lane}}.markdups.metrics"
    shell: """
    mkdir -p {WORKING_FOLDER}/tmp_{wildcards.sample}_{wildcards.lane}_dups
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bammarkduplicatesopt I={input.bam_od} O={output.bam_dupmarked} M={log.metrics} tmpfile={WORKING_FOLDER}/tmp_{wildcards.sample}_{wildcards.lane}_dups"
    rm -rf {WORKING_FOLDER}/tmp_{wildcards.sample}_{wildcards.lane}_dups
    """

rule append_rb_tag_filter:
    input:
        bam_dupmarked=f"{WORKING_FOLDER}/data/align/{{sample}}-{{lane}}_dupMarked.bam"
    output:
        # This is the final BAM per lane before merging
        bam_final_lane=f"{WORKING_FOLDER}/data/align/{{sample}}-{{lane}}_final.bam"
    threads: 1
    log: f"{WORKING_FOLDER}/00_log/{{sample}}-{{lane}}.bamaddreadbundles"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "bamaddreadbundles -I {input.bam_dupmarked} -O {output.bam_final_lane} 2> {log}"
    """


# Function to get the list of final BAM files for all lanes of a given sample
def get_lane_bams(wildcards):
    """
    Collects all final BAM files (after add read bundles) for the lanes associated with a sample.
    """
    # It's safer to make a copy before modifying
    sample_lane_dict = FILES[wildcards.sample].copy()
    # Remove matched normal keys if they exist, handle potential KeyError
    sample_lanes = sorted(sample_lane_dict.keys()) # Get tumor lanes for this sample

    # Construct the path for each lane's final BAM file
    # Corrected path construction: no leading '/'
    return expand(
        os.path.join(WORKING_FOLDER, "data/align", "{sample}-{lane}_final.bam"),
        sample=wildcards.sample,
        lane=sample_lanes,
        WORKING_FOLDER=WORKING_FOLDER # Pass WORKING_FOLDER explicitly if needed inside expand
    )


rule merge_bams:
    input:
        bams=get_lane_bams
    output:
        merged_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam"
    threads: 8
    message: "Merging final lane BAMs for sample {wildcards.sample}"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.merge_bams"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c \
    "samtools merge -@ {threads} -o {output.merged_bam} {input.bams} 2> {log}"
"""


rule index_tumor_bam:
    input:
        merged_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam"
    output:
        bai=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam.bai"
    threads: 4 # Should be <= samtools index -@ threads
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.samtoolsTumIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "samtools index -@ 4 {input.merged_bam} -o {output.bai} 2> {log}"
    """

rule deduplicate_bam:
    input:
        merged_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam" # Takes merged BAM as input
    output:
        dedup_bam=f"{WORKING_FOLDER}/data/dedup/{{sample}}.bam"
    threads: 1
    message: "Carrying out deduplication (random read in bundle)"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.randomreadinbundle"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "randomreadinbundle -I {input.merged_bam} -O {output.dedup_bam} 2>{log}"
    """

rule index_dedup_bam:
    input:
        dedup_bam=f"{WORKING_FOLDER}/data/dedup/{{sample}}.bam"
    output:
        bai=f"{WORKING_FOLDER}/data/dedup/{{sample}}.bam.bai"
    threads: 4 
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.samtoolsDedupIndex"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "samtools index -@ {threads} {input.dedup_bam} -o {output.bai} 2> {log}"
    """

rule analyse_efficiency:
    input:
        dedup_bam=f"{WORKING_FOLDER}/data/dedup/{{sample}}.bam",
        merged_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam", 
        dedup_bai=f"{WORKING_FOLDER}/data/dedup/{{sample}}.bam.bai" 
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_efficiency.txt")
    threads: 12
    message: "Calculating efficiency metrics"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.efficiency_pl"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "{CODE_HOME}/perl/efficiency_nanoseq.pl -t {threads} -dedup {input.dedup_bam} -duplex {input.merged_bam} -o {wildcards.sample} -r {BWA_INDEX} 2>{log}"
    """


rule nanoseq_cov:
    input:
        dedup_bam=f"{WORKING_FOLDER}/data/dedup/{{sample}}.bam",
        merged_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam",
        dedup_bai=f"{WORKING_FOLDER}/data/dedup/{{sample}}.bam.bai",
        merged_bai=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam.bai"
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_cov.txt")
    threads: 16
    message: "Nanoseq_coverage"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.nanoseq_cov"
    shell: """
    mkdir -p {WORKING_FOLDER}/{wildcards.sample} # Create sample-specific output dir if needed by script
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input.dedup_bam} -B {input.merged_bam} -R {BWA_INDEX} cov -Q 0 --exclude 'MT,GL%%,NC_%,hs37d5' 2>{log}"
    """

rule nanoseq_partition:
    input:
        norm_bam=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam",
        tumor_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam", # Use the merged tumor BAM
        cov_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_cov.txt", # Trigger dependency
        norm_bai=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam.bai",
        tumor_bai=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam.bai"
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_part.txt")
    threads: 1
    message: "Nanoseq_partition"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.nanoseq_partition"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t 1 -A {input.norm_bam} -B {input.tumor_bam} -R {BWA_INDEX} part -n {partitions} 2>{log}"
    """

rule nanoseq_dsa:
    input:
        norm_bam=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam",
        tumor_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam", # Use the merged tumor BAM
        part_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_part.txt", # Trigger dependency
        norm_bai=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam.bai",
        tumor_bai=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam.bai"
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_dsa.txt")
    threads: partitions
    message: "Nanoseq_dsa"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.nanoseq_dsa"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input.norm_bam} -B {input.tumor_bam} -R {BWA_INDEX} dsa -C SNP.sorted.bed.gz -D noise_mask_gilad_nanoseq.bed.gz -d 2 -q 30 2>{log}"
    """

rule nanoseq_var:
    input:
        norm_bam=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam",
        tumor_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam", # Use the merged tumor BAM
        dsa_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_dsa.txt", # Trigger dependency
        norm_bai=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam.bai",
        tumor_bai=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam.bai"
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_var.txt")
    threads: partitions
    message: "Nanoseq_var"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.nanoseq_var"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input.norm_bam} -B {input.tumor_bam} -R {BWA_INDEX} var \
    -a 50 -b 5 -c 0 -f 0.9 -i 1 -m 8 -n 3 -p 0 -q 60 -r 144 -v 0.01 -x 8 -z 12 2>{log}"
    """ 

rule nanoseq_indel:
    input:
        norm_bam=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam",
        tumor_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam", # Use the merged tumor BAM
        dsa_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_dsa.txt", # Trigger dependency
        norm_bai=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam.bai",
        tumor_bai=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam.bai"
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_indel.txt")
    threads: partitions
    message: "Nanoseq_indel"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.nanoseq_indel"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t {threads} -A {input.norm_bam} -B {input.tumor_bam} -R {BWA_INDEX} indel -s {wildcards.sample} \
    --rb 2 --t3 140 --t5 5 2>{log}"
    """

rule nanoseq_post:
    input:
        norm_bam=f"{WORKING_FOLDER}/data/align_norm/{{sample}}.bam",
        tumor_bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam", # Use the merged tumor BAM
        # Trigger dependencies from previous steps
        indel_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_indel.txt",
        var_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_var.txt",
        cov_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_cov.txt",
        dsa_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_dsa.txt",
        part_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_part.txt"
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_post.txt")
    threads: 2
    message: "Nanoseq_post"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.nanoseq_post"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} -B {BWA_FOLDER} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "python {CODE_HOME}/python/runNanoSeq.py --out {WORKING_FOLDER}/{wildcards.sample} -t 2 -A {input.norm_bam} -B {input.tumor_bam} -R {BWA_INDEX} post 2>{log}"
    """

rule nanoseq_plot_variant_qc:
    input:
        post_trigger=f"{WORKING_FOLDER}/{{sample}}_nanoseq_post.txt" # Trigger dependency
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_variant_qc.txt")
    threads: 1
    params:
        partitions=partitions
    message: "Plot graphs to assess quality of variant calling and filtering"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.varaint_qc"
    shell: """
    singularity exec -B {FASTQ_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c "Rscript {CODE_HOME}/R/variant_filtering_qc.R {WORKING_FOLDER}/{wildcards.sample}/tmpNanoSeq/post/discardedvariants.csv \
    {WORKING_FOLDER}/{wildcards.sample}/tmpNanoSeq/post/variants.csv {WORKING_FOLDER}/{wildcards.sample}/tmpNanoSeq {params.partitions} 2>{log}"
    """


rule verify_cell_origin:
    input:
        bam=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam",
        bai=f"{WORKING_FOLDER}/data/filtered/{{sample}}_merged.bam.bai"
    output:
        touch(f"{WORKING_FOLDER}/{{sample}}_nanoseq_verify.txt")
    threads: 1
    message: "Verify bam ID"
    log: f"{WORKING_FOLDER}/00_log/{{sample}}.verifyBamId"
    shell: """
    singularity exec -B {BWA_FOLDER} -B {WORKING_FOLDER} -B {CODE_HOME} \
    --env "PATH=/home/ubuntu/miniforge3/bin:/opt/wtsi-cgp/bin:${{PATH}}" nanoseq_1.1.0.sif /bin/bash -c \
    "VerifyBamID --SVDPrefix {VERIFYBAMID_HOME}/resource/1000g.phase3.100k.b38.vcf.gz.dat \
    --Reference {BWA_INDEX} --BamFile {input.bam} --Output {WORKING_FOLDER}/{wildcards.sample}.verifyBamID 2>{log}"
    # Note: VerifyBamID usually outputs to --Output prefix. The touch() might be just a trigger.
    # The actual verify output file is {WORKING_FOLDER}/{wildcards.sample}.verifyBamID.selfSM
    # Adjusted the output in TARGETS/rule all if needed, or adjust this rule's output.
    # For now, keeping touch() as it was in TARGETS.
    """
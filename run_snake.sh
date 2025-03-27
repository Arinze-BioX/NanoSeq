#! /bin/bash

timestamp=$(date +%Y-%m-%d_%H-%M-%S)
filename1="nanoseq_snakemake_run_out_${timestamp}.log"
filename2="nanoseq_snakemake_run_error${timestamp}.log"

snakemake --snakefile snakemake_pipeline/nanoseq.smk --profile snakemake_pipeline/slurm_singularity --rerun-incomplete --keep-going 1> "${filename1}" 2> "${filename2}"
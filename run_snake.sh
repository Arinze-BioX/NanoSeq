#! /bin/bash

timestamp=$(date +%Y-%m-%d_%H-%M-%S)
filename2="nanoseq_snakemake_run_error${timestamp}.log"

snakemake --snakefile snakemake_pipeline/nanoseq.smk --profile snakemake_pipeline/slurm_singularity --rerun-incomplete --keep-going 2> "${filename2}"
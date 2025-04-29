#!/bin/bash


snakemake -p -s snakemake_pipeline/nanoseq.smk --forceall --rulegraph | dot -Tpdf > snakemake_pipeline/dag.pdf
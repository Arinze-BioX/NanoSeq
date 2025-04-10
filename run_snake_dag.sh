#!/bin/bash


snakemake -p -s snakemake_pipeline/nanoseq.smk --forceall --rulegraph | dot -Tpng > snakemake_pipeline/dag.png
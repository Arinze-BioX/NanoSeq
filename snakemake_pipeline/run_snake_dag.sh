#!/bin/bash


snakemake -p -s nanoseq.smk --forceall --rulegraph | dot -Tpng > dag.png
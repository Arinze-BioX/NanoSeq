# NanoSeq
[![Nanoseq Docker](https://img.shields.io/docker/pulls/mtlynch/logpaste.svg?maxAge=604800)](https://hub.docker.com/repository/docker/arinzeokafor/nanoseq/general)

Nanorate sequencing (NanoSeq) is a DNA library preparation and sequencing protocol based on Duplex Sequencing ([Schmitt et al, 2012](https://doi.org/10.1073/pnas.1208715109)) and BotSeqS ([Hoang et al, 2016](https://doi.org/10.1073/pnas.1607794113)). NanoSeq allows calling mutations with single molecule resolution and extremely low error rates ([Abascal et al, 2021](https://doi.org/10.1038/s41586-021-03477-4)). The pipeline and code in this repository cover the preprocessing of NanoSeq sequencing data, the assessment of data quality and efficiency, and the calling of mutations (substitutions and indels) and the estimation of mutation burdens and substitution profiles.

### Upstream/Acknowledgement
This repository was forked from https://github.com/cancerit/NanoSeq and modified (see committed changes).

### Dependencies
- singularity=3.9.8
- snakemake=9.1.1
- snakemake-executor-plugin-slurm=1.1.0

### Usage
This pipeline can be used from the already created container image, by:
1. Cloning or downloading this repository
2. Pulling [the image](https://hub.docker.com/repository/docker/arinzeokafor/nanoseq) from dockerhub. If using singularity, let the .sif file be in the parent NanoSeq directory after cloning, by running the below:
```
git clone git@github.com:Arinze-BioX/NanoSeq.git #clones the directory
cd NanoSeq
singularity pull docker://arinzeokafor/nanoseq:1.1.0
```
Of course, this would require locally installing git as well as whatever software is being used to pull the image and run the container. I used singularity=3.9.8.

3. Modify the pipeline as follows:
    1. modify the snakemake_pipeline/config_file.yaml file to ensure the variables point to the correct directories in your local system
    2. modify the snakemake_pipeline/slurm_singularity/config.yaml file according to the job scheduler you intend to use for parallel processing. Also install the necessary snakemake plug-ins to parallelize the pipeline or run on hpc. I used slurm (https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
    3. modify the snakemake_pipeline/samples.json file to include your desired input file paths. Note that the first-level keys in the json file are for identifying biological samples and the second-level keys are lanes or sequencing runs per biological sample. The alignment files for different lanes per sample are eventually merged, after optical duplicate marking. This is crucial for the correct estimation of PCR duplicates. The final outputs will be per biological sample. The pipeline will run without errors even if each sample has only one sequencing run.
4. Run the run_snake.sh (after making it executable) file from within the cloned parent NanoSeq directory.

#### Note: 
The run_snake.sh file would need to be run in a slurm context if running on a hpc, due to local memory requirement. The terminal where run_snake.sh will also need to be protected from time-out disconnection either running it through sbatch or in a tmux session.

Please find more details about this pipeline in the README of the original pipeline and the associated papers.
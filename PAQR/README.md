# Run the PAQR pipeline

PAQR is a tool that allows the quantification of transcript 3' ends (or poly(A) sites) based on standard RNA-seq data.

The method is actually a combination of scripts that depend on each other. To ease the execution of the scripts they are assembled into a pipeline. The execution of the pipeline is very simple. 

## Requirements
Depending on how deep your original material was sequenced (i.e. how big your input files are in terms of the memory usage), your computer might require **50 GB of available memory (RAM)** to process these files. [Snakemake]{https://snakemake.readthedocs.io/en/stable/} is used to facilitate the automated execution of all analysis steps. The easiest way to make use of the pipeline is to set up a python3 and a python2 virtual environment and run the pipeline in the python3 environment. In the following section, instructions on how to install the virtual environment for the analysis are given.

### snakemake
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires python3 and can be most easily installed via the bioconda package of the python anaconda distribution. You're setup right after the following steps:

### activate the python 3 environment
Please see the README file one level up. There you can find information how to install miniconda and how to create a python3 environment that is used now to run the pipeline.

Activate the environment:
  ```bash
  source activate paqr_kapac
  ```


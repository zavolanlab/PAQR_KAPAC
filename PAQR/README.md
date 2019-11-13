# Run the PAQR pipeline

PAQR is a tool that allows the quantification of transcript 3' ends (or poly(A) sites) based on standard RNA-seq data. As input it requires alignment files in BAM format and it returns a table of quantified tandem poly(A) sites (i.e. poly(A) sites that belong to the same gene and change the length of the 3' UTR).

The method is actually a combination of scripts that depend on each other. To ease the execution of the scripts they have been assembled into a pipeline. The use PAQR, you simply have to run the pipeline. 

## Requirements
For normal samples (e.g. the samples that are used as test cases here) the pipeline will run on a laptop or stationary working machine with at least 10 GB of memory (RAM). Be aware though that the memory consumptions depends on the sequencing depth of your samples (i.e. how large your input files are). In rare cases, we found a peak memory usage of around 50 GB for very deeply sequenced samples (more than 400 million reads). [Snakemake](https://snakemake.readthedocs.io/en/stable/) is used to facilitate the automated execution of all analysis steps. The easiest way to make use of the pipeline is to set up a python3 and a python2 virtual environment and run the pipeline in the python3 environment. Instructions on how to install the virtual environments for the analysis are given in the README of the main directory of this repository.

### snakemake
[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a workflow management system that helps to create and execute data processing pipelines. It is already part of the virtual environment which was introduced in the introduction README.

### activate the python 3 environment
Please see the README file one level up. There you can find information how to install miniconda and how to create a python3 environment that is used now to run the pipeline.

Activate the environment:
  ```bash
  source activate paqr_kapac
  ```

### Configure the input parameters
The PAQR subdirectory (you should be in this directory while you go through this README) contains a file called "config.yaml". This files contains all information about used parameter values, data locations, file names and so on. During a run, all steps of the pipeline will retrieve their paramter values from this file. It follows the yaml syntax (find more information about yaml and it's syntax [here](http://www.yaml.org/)) what makes it easy to read and edit. The main principles are:
  - everything that comes after a `#` symbol is considered as comment and will not be interpreted
  - paramters are given as key-value pair, with `key` being the name and `value` the value of any paramter

Before starting a run of the pipeline, open "config.yaml" in a text editor of your choice, e.g.
on Linux:
  ```bash
  gedit config.yaml
  ```

on MaxOS X:

  ```bash
  open -e config.yaml
  ```

Some entries require your editing while most of them you can leave unchanged. However, this config file contains all parameters used in the pipeline and the comments should give you the information about their meaning. If you need to change path names please ensure to **use relative instead of absolute path names**.

The first value you have to change is "py2_env_path" and the value should be set to the path where you installed your python2 environment before. Since you're in the activated python3 environment at this moment, it is easy to figure out the path. Just type
  ```bash
  which python
  ```

The output of this command should look similar to `~/miniconda3/envs/paqr_kapac/bin/python`. From this output, copy the first part up to "envs" and replace it with the value in the config file given for "py2_env_path".

Next, define a name for the study your samples belong to under "studies". This name is used for the output directory of the run. It is possible to process the samples of multiple studies at the same time (even if they share input files, e.g. a bam file from a control sample). In this case, provide a comma-separated list of study names.

Connect the name of your samples with the study they belong to: The config file should have one entry per study with the study name as key. Follow the example in the config file and adopt the list for "samples". This list should contain the names you use for your input samples (make sure the sample names are unique, though). Again, each sample name should have its own entry in the config file with the name of the BAM and the condition/type of the sample given as a dictionary. The value for "bam" must be the name of the BAM file (without ".bam" extension and without the pathname). All BAM files are required to be stored in the same directory. Provide the pathname to this directory as value of "dir.input".


## Start the first part of the pipeline
Before you run the pipeline, ensure that your current working directory is PAQR.
Once you prepared your config file, you can start the pipeline by simply executing:
```bash
max_cores=4 # maximum number of threads that will run in parallel
snakemake -s part_one.Snakefile -p --cores ${max_cores}
```
or, if you like to save the snakemake comments in a log-file (recommended):
```bash
max_cores=4 # maximum number of threads that will run in parallel
snakemake -s part_one.Snakefile -p --cores ${max_cores} &> log_output.log
```

It is recommend to set the `max_cores` parameter so that mutliple steps of the pipeline can run in parallel. Moreover, some scripts of the pipeline make use of a parallel environment and defining a maximum number of available cores prevents overloading the machine.

After the first part is finished, the transcript integrity was evaluated and only samples above the cutoff for the median TIN per sample (mTIN) will be processed in the second part.

## Second part of the pipeline
The second part can be started immediately after successful finishing the first part. However, please consider the following caveat: It might occur that your study (or a single study of multiple ones that you run in parallel) did not yield any valid sample. In this case, manually check results directory of all studies for a file called `DUMMY_USED.out`; delete all study names from the "studies" entry in the config for which such a dummy file is available; afterwards, proceed with part two and it will only run on those studies that have valid files.

```bash
max_cores=8 # maximum number of threads that will run in parallel
snakemake -s part_two.Snakefile -p --cores ${max_cores} &>> log_output.log
```

## Use case tutorial
Let's run PAQR on an RNA-seq data set from a study of HNRNPC in HEK cells (you can find the publication [here](https://www.ncbi.nlm.nih.gov/pubmed/25719671) and the data is deposited [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56010)). 

The data was downloaded and mapped to the human genome with [STAR v2.5.2a](https://github.com/alexdobin/STAR). First, download the mapping files (in BAM format) and an index for each mapping file from <http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/>, for example with the following commands:
```bash
mkdir data/bam_files
# control replicate 1
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/CTL_rep1.bam
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/CTL_rep1.bam.bai
# control replicate 2
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/CTL_rep2.bam
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/CTL_rep2.bam.bai
# HNRNPC knock-down replicate 1
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/KD_rep1.bam
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/KD_rep1.bam.bai
# HNRNPC knock-down replicate 2
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/KD_rep2.bam
wget -P data/bam_files/ http://www.clipz.unibas.ch/RNAseq_HNRNPC_KD_study/data/KD_rep2.bam.bai
```
If not yet done, open the configuration file `config.yaml` and adjust the value for "py2_env_path" (as described above). All other values are already set up for this test case. Start the first part of the analysis with:
```bash
max_cores=8 # maximum number of threads that will run in parallel
snakemake -s part_one.Snakefile -p --cores ${max_cores} &> log_output.log
```
When the first part is finished, start the second part with:
```bash
max_cores=8 # maximum number of threads that will run in parallel
snakemake -s part_two.Snakefile -p --cores ${max_cores} &>> log_output.log
```
Adjust the value max_cores according to your architecture. The number of cores of your machine can be obtained with `sysctl -n hw.ncpu` for Max OS X and with `nproc` for Linux systems. 
After successfully finishing the pipeline, the results directory ("HNRNPC_KD") contains a file called "tandem_pas_expressions.rpm.filtered.tsv", which is the starting file for a KAPAC analysis.


## Detailed description of the single steps
The following notes should provide more detailed information about the single steps/scripts of the pipeline.

### create_log_dir:
**input:** -

**output:** *temp, HNRNPC_KD/created_log.tmp

**details:** Create a directory for log-entries written by HPC jobs (only used when the pipeline is run on a HPC cluster). *local

**shell command:** pipeline inherent python code

### bam_index:
**input:** data/bam_files/KD_rep1.bam

**output:** data/bam_files/KD_rep1.bam.bai

**details:** If not available, create an index for each input BAM file.

**shell command:** samtools index {input.bam} {output.bam_index}

### assess_TIN_bias_transcript_wide:
**input:** data/bam_files/KD_rep1.bam, data/bam_files/KD_rep1.bam.bai

**parameters:**
  - transcripts: set of transcripts in BED format for which the TIN score is calculated
  - min_raw_reads: minimum number of raw reads to overlap with a transcripts; otherwise the transcript is not evaluated
  - sample_size: number of randomly selected positions within the transcript body used together with all exon start and stop positions for the TIN calculation

**output:** data/samplewise_TIN/KD_rep1.tsv

**shell command:**
```bash
transcripts="data/annotation/full_transcripts.hg38.canonical_chr.tandem.noOverlap_strand_specific.bed"
min_raw_reads=10
sample_size=100

python calculate-TIN-values.py \
  -i data/bam_files/KD_rep1.bam \
  -r ${transcripts} \
  -c ${min_raw_reads} \
  --names KD_rep1 \
  -n ${sample_size} \
  > data/samplewise_TIN/KD_rep1.tsv
```

## merge_TIN_vals:
**input:** data/samplewise_TIN/KD_rep1.tsv, data/samplewise_TIN/CTL_rep1.tsv

**output:** HNRNPC_KD/bias.transcript_wide.TIN.tsv

**details:** Combine the TIN values from all samples into one table; *local

**shell command:**
```bash
python merge-TIN-tables.py \
  --verbose \
  --input data/samplewise_TIN/KD_rep1.tsv data/samplewise_TIN/CTL_rep1.tsv \
  > HNRNPC_KD/bias.transcript_wide.TIN.tsv
```

## TIN_assessment:
**input:** HNRNPC_KD/bias.transcript_wide.TIN.tsv

**parameters:**
  - bias_cutoff: median TIN (mTIN) score cutoff; samples below this cutoff are excluded from further analysis

**output:** HNRNPC_KD/bias.TIN.median_per_sample.tsv, HNRNPC_KD/no_bias_samples.out

**details:** The mTIN value for each sample is obtained and written to a list. The names of all samples above the mTIN cutoff are stored in HNRNPC_KD/no_bias_samples.out. When all samples are below the cutoff, another file, "HNRNPC_KD/DUMMY_USED.out", is written to indicate that the samples in HNRNPC_KD/no_bias_samples.out are only dummy placeholders. *local

**shell command:** pipeline inherent python code

## plot_full_trans_TIN_distributions:
**input:** HNRNPC_KD/bias.TIN.median_per_sample.tsv

**output:** HNRNPC_KD/bias.transcript_wide.TIN.boxplots.pdf

**details:** Creates a boxplot for the TIN value distribution of each sample; *local

**shell command:**
```bash
Rscript boxplots-TIN-distributions.R \
  --file HNRNPC_KD/bias.TIN.median_per_sample.tsv\
  --pdf HNRNPC_KD/bias.transcript_wide.TIN.boxplots.pdf
```

## infer_valid_samples:
**input:** HNRNPC_KD/no_bias_samples.out

**output:** HNRNPC/valid_samples/KD_rep1.tsv, HNRNPC/valid_samples/CTL_rep1.tsv

**details:** Creates softlinks for those *bam and *bai files that belong to valid samples. Otherwise, dummy bam files are linked.

**shell command:** pipeline inherent python code

---

*temp: indicates temporary files that are deleted as soon as the pipeline was run successfully. So, under normal conditions you will never find these files in your results directory.

*local: these steps are executed on the local machine, even when the pipeline is run in cluster mode (i.e. snakemake would submit the rule execution as job to a cluster computer

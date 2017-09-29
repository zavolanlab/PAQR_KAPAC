# KAPAC (K-mer Activity on Poly(A) site Choice)

KAPAC, standing for K-mer Activity on Polyadenylation Site Choice, aims to identify sequence motifs (k-mers) that can explain changes in poly(A) site usage across conditions (e.g. in samples in which the expression of a potential regulator has been perturbed). It models changes in poly(A) site usage as a linear function of the occurrence of specific k-mers and the unknown regulatory impact (also called "activity") of k-mers.

## Installation requirements and dependencies
KAPAC was developed and tested in R version 3.3.1 and needs an installation of the R "optparse" library. However, other R versions might work as well. Especially when you are interested to run KAPAC in combination with PAQR we recommend to run KAPAC within the virtual environment that needs to be set up for PAQR. Please see the README file one level up. There you can find information how to install miniconda and how to create a python3 environment that is used now to run the pipeline.

Activate the environment:
  ```bash
  source activate paqr_kapac
  ```

## Memory and runtime requirements
The memory and runtime scales with the input data. A detault run for all k-mers of length 3 to 6 on typical A-seq experiments (example data provided) can be performed using less than 8 GB RAM. The chosen number of runs on randomized expression data (number_of_randomized_runs, see below) will have a big influence on the runtime of KAPAC. Using default settings (number_of_randomized_runs=1000) a typical KAPAC run will take several hours. In order to get a first glimpse on data within less than an hour of runtime one can run KAPAC with a smaller number of randomized runs (e.g. number_of_randomized_runs=100). However, please note that the results might be less reproducible because of the small number of randomizations that are used for calculating the adjusted p-value for mean difference z-scores ("mean_diff_zscores_PADJ" column in the output file). Finally, even thouth it is recommended to consider adjusted p-values for mean difference z-scores, it is possible to run KAPAC without calculating them by setting the number of randomized runs below the minimum number required (=30). Running KAPAC without randomizing the expression data (number_of_randomized_runs=0) takes only a few minutes to finish. Please note that KAPAC will drop a warning message when running it without randomized runs.

## How to run KAPAC
Simple example for running KAPAC on the example data provided in the "DATA" directory:

```bash
Rscript --vanilla KAPAC.R \
  --sample_design DATA/kapac_design.tsv \
  --rowname_to_group_mapping_file DATA/pas2exon.tsv \
  --expression_matrix DATA/pas_expression.tsv \
  --pas_overlap PAS_overlap \
  --sitecount_matrix DATA/kmer_counts_all.tsv \
  --considered_region_length 50 \
  --results_dir RESULTS
```

Detailed example on how to run KAPAC on a list of selected k-mers (using the example data in the "DATA" directory):

```bash
Rscript --vanilla KAPAC.R \
  --sample_design DATA/kapac_design.tsv \
  --expression_matrix DATA/pas_expression.tsv \
  --sitecount_matrix DATA/kmer_counts.tsv \
  --rowname_to_group_mapping_file DATA/pas2exon.tsv \
  --selected_motifs DATA/selected_kmers.tsv \
  --results_dir RESULTS \
  --create_plots_for_each_motif TRUE \
  --create_files_for_each_motif TRUE \
  --row_center_expression TRUE \
  --treat_minus_ctrl FALSE \
  --expression_pseudocount 1 \
  --consider_excess_counts_only TRUE \
  --considered_region_length 50 \
  --min_kmer_abundance_fraction 0.01 \
  --number_of_randomized_runs 1000 \
  --pas_overlap PAS_overlap \
  --verbose TRUE \
&> RESULTS.log
```


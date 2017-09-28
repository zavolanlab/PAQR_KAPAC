# K-mer Activity on Polyadenylation Site Choice (KAPAC)

KAPAC, standing for K-mer Activity on Polyadenylation Site Choice, aims to identify sequence motifs (k-mers) that can explain changes in poly(A) site usage across conditions (e.g. in samples in which the expression of a potential regulator has been perturbed). It models changes in poly(A) site usage as a linear function of the occurrence of specific k-mers and the unknown regulatory impact (also called "activity") of k-mers.

## Installation requirements and dependencies
KAPAC was developed and tested in R version 3.3.1 and needs an installation of the R "optparse" library. However, other R versions might work as well. Especially when you are interested to run KAPAC in combination with PAQR we recommend to run KAPAC within the virtual environment that needs to be set up for PAQR. Please see the README file one level up. There you can find information how to install miniconda and how to create a python3 environment that is used now to run the pipeline.

Activate the environment:
  ```bash
  source activate paqr_kapac
  ```

## Run KAPAC
Once all input files have been prepared KAPAC can be ran as follows:

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
  --number_of_randomized_runs 100 \
  --pas_overlap PAS_overlap \
  --verbose TRUE \
&> RESULTS.log
```


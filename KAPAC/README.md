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
Please get help on KAPAC options by running:
```bash
Rscript --vanilla KAPAC.R --help
```

### Required input files
KAPAC requires the following four input files:

1. __Sample design file (option: --sample_design)__: A tabulator separated values (TSV) file that contains information about the samples and contrasts that will be used in the analysis. See also the example sample design file provided on GIT (DATA/kapac_design.tsv). The file must contain one line per sample providing (at least) the following columns:
  * __column: "sample_id"__: Contains a unique id for each sample.
  * __column: "group"__: The name of the group to which a sample belongs (e.g. CNTRL for the control group and KD for the knock-down group).
  * __column: "contrast"__: Must contain the character string "CNTRL" for samples that will be used as contrast controls. For contrast treatment samples this column must contain the unique sample id (see column "sample_id") of the sample that should be used as control for the treatment sample.

2. __Poly(A) site / exon associations file (option: --pas_exon_associations)__: A tabulator separated values (TSV) file that contains unique poly(A) site identifiers in a column named "pas" and unique exon identifiers in a column named "exon". See also the example poly(A) site / exon associations file provided on GIT (DATA/pas2exon.tsv)

3. __Poly(A) site expression file (option: --expression_matrix)__: A tabulator separated values (TSV) file that contains unique poly(A) site identifiers in a column named "pas". Another column named "PAS_overlap" must contain 'OK' for non overlapping poly(A) sites and 'OVERLAP' otherwise. Overlapping means, that they are not distant from each other. It is recommended to use at least twice the length of the region that is investigated with KAPAC in order to prevent confounding effects from k-mers that belong to other poly(A) sites located on the same exon. In case all poly(A) sites should be considered this columns should contain for each row 'OK' (not recommended). Further columns must contain the expression of the poly(A) sites within a specific sample, whereat the column name must be the unique sample identifier used in the sample design file (option: --sample_design). See also the example expression file provided on GIT (DATA/pas_expression.tsv)

4. __Sitecounts file (option: --sitecount_matrix)__: A tabulator separated values (TSV) file that contains unique poly(A) site identifiers in a column named "pas". All other columns contain a unique motif/k-mer name and the corresponding number of sites obsvered within a defined region relative to the poly(A) site in the row. See also the example sitecounts file provided on GIT (DATA/kmer_counts_all.tsv).

### KAPAC options
KAPAC takes the following options:

* __--sample_design [REQUIRED]__: The sample design file (see above).

* __--treat_minus_ctrl [DEFAULT=FALSE]__: If true, KAPAC will consider treatment minus control for the fit and control minus treatment otherwise.

* __--pas_exon_associations [REQUIRED]__: The file that contains the associations of poly(A) sites to the exon on which they are located (see above).

* __--expression_matrix [REQUIRED]__: The poly(A) site expression matrix (see above).

* __--expression_pseudocount [DEFAULT=1.0]__: The pseudocount that should be add to the expression values prior going to log-space.

* __--row_center_expression [DEFAULT=TRUE]__: If set, the expression matrix (see option --expression_matrix) will be row-centered. That is, changes in relative usage across samples will be explained.

* __--sitecount_matrix [REQUIRED]__: The site count matrix (see above).

* __--min_kmer_abundance_fraction [DEFAULT=0.01]__: The fraction of poly(A) sites that needs to contain counts for a specific k-mer in order to be considered. That is, k-mers that have counts for a smaller fraction of poly(A) sites in the sitecount matrix (see option --sitecount_matrix) will not be considered in the KAPAC run. 

* __--consider_excess_counts_only [DEFAULT=TRUE]__: If set, background correction will be performed. That is, site counts will be corrected by the number of counts that are expected to be found per chance given the considered region length (see option --considered_region_length).

* __--considered_region_length [REQUIRED]__: The length of the region in which the k-mers have been counted (e.g. needed for background correction, see option --consider_excess_counts_only).

* __--selected_motifs [DEFAULT=all]__: Per default, all k-mers present in the sitecount matrix (see option --sitecount_matrix) found in a high enough fraction of poly(A) sites (see option --min_kmer_abundance_fraction) will be considered in the KAPAC run. Alternatively, a file can be provided which contains the k-mers that should be considered in a column named "motif".

* __--create_plots_for_each_motif [DEFAULT=FALSE]__: If set, for each kmer plots will be created (NOTE: dependent on the sitecount matrix (see option --sitecount_matrix) and the selected k-mers (see option --selected_motifs), thousands of directories and files might be created).

* __--create_files_for_each_motif [DEFAULT=FALSE]__: If set, for each kmer detailed files (activities, errors, z-scores)  will be created (NOTE: dependent on the sitecount matrix (see option --sitecount_matrix) and the selected k-mers (see option --selected_motifs), thousands of directories and files might be created).

* __--number_of_randomized_runs [DEFAULT=1000]__: The number of runs done with randomized expression to k-mer count associations in order to calculate adjusted p-values for the reported activity difference z-scores.

* __--verbose [DEFAULT=TRUE]__: If set, the script will be verbose (reporting detailed infos on what is done).

* __--results_dir [REQUIRED]__: The directory to which the results will be written.

* __--help__: Show the KAPAC help message and exit

### KAPAC example runs

Example on how to run KAPAC on the data provided in the "DATA" directory using default options:

```bash
Rscript --vanilla KAPAC.R \
  --sample_design DATA/kapac_design.tsv \
  --pas_exon_associations DATA/pas2exon.tsv \
  --expression_matrix DATA/pas_expression.tsv \
  --sitecount_matrix DATA/kmer_counts.tsv \
  --considered_region_length 50 \
  --results_dir RESULTS \
&> RESULTS.log
```

Example on how to run KAPAC on the data provided in the "DATA" directory using a list of selected k-mers and customized options:

```bash
Rscript --vanilla KAPAC.R \
  --sample_design DATA/kapac_design.tsv \
  --expression_matrix DATA/pas_expression.tsv \
  --sitecount_matrix DATA/kmer_counts.tsv \
  --pas_exon_associations DATA/pas2exon.tsv \
  --selected_motifs DATA/selected_kmers.tsv \
  --results_dir RESULTS \
  --create_plots_for_each_motif TRUE \
  --create_files_for_each_motif TRUE \
  --row_center_expression TRUE \
  --treat_minus_ctrl FALSE \
  --expression_pseudocount 1.0 \
  --consider_excess_counts_only TRUE \
  --considered_region_length 50 \
  --min_kmer_abundance_fraction 0.01 \
  --number_of_randomized_runs 100 \
  --verbose TRUE \
&> RESULTS.log
```


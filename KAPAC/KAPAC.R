# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# KAPAC (K-mer activity on poly(A) site choice)
# -----------------------------------------------------------------------------
# Developed in R version: R version 3.3.1
# -----------------------------------------------------------------------------
rm(list=ls())

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read options
# -----------------------------------------------------------------------------
# load libraries
library(optparse)

option_list <- list(
  make_option(c("--sample_design"), action="store", type="character", help="The sample design file."),
  make_option(c("--expression_matrix"), action="store", type="character", help="The poly(A) site expression matrix."),
  make_option(c("--sitecount_matrix"), action="store", type="character", help="The site count matrix."),
  make_option(c("--pas_exon_associations"), action="store", type="character", help="The file that contains the associations of poly(A) sites to the exon on which they are located."),
  make_option(c("--selected_motifs"), default="all", action="store", type="character", help="Per default, all k-mers present in the sitecount matrix (see option --sitecount_matrix) found in a high enough fraction of poly(A) sites (see option --min_kmer_abundance_fraction) will be considered in the KAPAC run. Alternatively, a file can be provided which contains the k-mers that should be considered in a column named 'motif'."),
  make_option(c("--results_dir"), action="store", type="character", help="The directory to which the results will be written."),
  make_option(c("--create_plots_for_each_motif"), default=FALSE, action="store", type="logical", help="If set TRUE, for each kmer plots will be created (NOTE: dependent on the sitecount matrix (see option --sitecount_matrix) and the selected k-mers (see option --selected_motifs), thousands of directories and files might be created)."),
  make_option(c("--create_files_for_each_motif"), default=FALSE, action="store", type="logical", help="If set TRUE, for each kmer detailed files (activities, errors, z-scores)  will be created (NOTE: dependent on the sitecount matrix (see option --sitecount_matrix) and the selected k-mers (see option --selected_motifs), thousands of directories and files might be created)."),
  make_option(c("--row_center_expression"), default=TRUE, action="store", type="logical", help="If set TRUE, the expression matrix (see option --expression_matrix) will be row-centered. That is, changes in relative usage across samples will be explained."),
  make_option(c("--treat_minus_ctrl"), default=FALSE, action="store", type="logical", help="If true, KAPAC will consider treatment minus control for the fit and control minus treatment otherwise."),
  make_option(c("--expression_pseudocount"), default=1.0, action="store", type="double", help="The pseudocount that should be add to the expression values prior going to log-space."),
  make_option(c("--consider_excess_counts_only"), default=TRUE, action="store", type="logical", help="If set TRUE, background correction will be performed. That is, site counts will be corrected by the number of counts that are expected to be found per chance given the considered region length (see option --considered_region_length)."),
  make_option(c("--considered_region_length"), action="store", type="double", help="The length of the region in which the k-mers have been counted (e.g. needed for background correction, see option --consider_excess_counts_only)."),
  make_option(c("--min_kmer_abundance_fraction"), default=0.01, action="store", type="double", help="The fraction of poly(A) sites that needs to contain counts for a specific k-mer in order to be considered. That is, k-mers that have counts for a smaller fraction of poly(A) sites in the sitecount matrix (see option --sitecount_matrix) will not be considered in the KAPAC run."),
  make_option(c("--number_of_randomized_runs"), default=1000, action="store", type="double", help="The number of runs done with randomized expression to k-mer count associations in order to calculate adjusted p-values for the reported activity difference z-scores."),
  make_option(c("--verbose"), default=TRUE, action="store", type="logical", help="If set TRUE, the script will be verbose (reporting detailed infos on what is done)."))
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS]",
                           option_list = option_list, add_help_option=TRUE)
opt <- parse_args(opt_parser)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# for debugging
# -----------------------------------------------------------------------------
# do you want to run in debug mode?
debugging_mode = FALSE

# in case we want to run in debug mode, we set all variables
if (debugging_mode == TRUE)
{
  opt = list()
  opt$sample_design = "DATA/kapac_design.tsv"
  opt$expression_matrix = "DATA/pas_expression.tsv"
  opt$sitecount_matrix = "DATA/kmer_counts.tsv"
  opt$pas_exon_associations = "DATA/pas2exon.tsv"
  opt$selected_motifs = "DATA/selected_kmers.tsv"
  #opt$selected_motifs = "all"
  opt$results_dir = "RESULTS"
  opt$create_files_for_each_motif = TRUE
  opt$create_plots_for_each_motif = TRUE
  opt$row_center_expression = TRUE
  opt$treat_minus_ctrl = FALSE
  opt$expression_pseudocount = 1
  opt$consider_excess_counts_only = TRUE
  opt$considered_region_length = 50
  opt$min_kmer_abundance_fraction = 0.01
  opt$number_of_randomized_runs = 30
  opt$verbose = TRUE
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# variables for writing messages, warnings and errors
# -----------------------------------------------------------------------------
docs_80_unerlines = "_______________________________________________________________________________"
docs_80_dashes = "-------------------------------------------------------------------------------"

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Create the nucleotides frequency vector
# -----------------------------------------------------------------------------
# In case we want to background correct, we will use the following
# frequencies and report it to the user
nt_frequency_vector = c(0.2973, 0.1935, 0.2007, 0.3084)
names(nt_frequency_vector) = c("A","C","G","T")

# -----------------------------------------------------------------------------
# we require to perform at least 30 random runs in order to get p-values 
# calculated based on a z-statistic (see: https://en.wikipedia.org/wiki/Z-test)
min_zstatistic_sample_size = 30
recommended_zstatistic_sample_size = 1000

# -----------------------------------------------------------------------------
# what should be used to rank final results
top_kmer_candidates_selection_criterion = "mean_diff_zscores"

# -----------------------------------------------------------------------------
# The threshold that will be used in order to consider a distribution to be
# normal.
not_norm_pval_threshold = 0.05

# -----------------------------------------------------------------------------
# The name of the pas overlapping column in the expression file.
pas_overlap_col = "PAS_overlap"

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# create the results dir
# -----------------------------------------------------------------------------
dir.create(opt$results_dir, showWarnings = FALSE, recursive = TRUE)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# determine if we should plot treatment minus control or
# control minus treatment
# -----------------------------------------------------------------------------
if (opt$treat_minus_ctrl == TRUE) {
  treat_idx = 1
  ctrl_idx = 2
  ylabel = 'activity difference (treat - ctrl)'
} else {
  treat_idx = 2
  ctrl_idx = 1
  ylabel = 'activity difference (ctrl - treat)'
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: create_design_table
# ARGUMENTS: sample_design_file_path
# DESCRIPTION: Reads in the design table from a given path          
# -----------------------------------------------------------------------------
create_design_table <- function(sample_design_file_path)
{
  # read the table and create row names
  design_table = 
    read.table(sample_design_file_path, 
               h=TRUE, 
               as.is=T, 
               sep="\t",
               comment.char="",
               check.names=FALSE)
  design_table$name = apply(design_table, 1, function(x) {paste(x['group'], x['sample_id'], sep='_')} )
  rownames(design_table) = design_table$name
  
  # return the design table
  return(design_table)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: create_contrast_pairs
# ARGUMENTS: design_table
# DESCRIPTION: From a given design table, the function creates a list of
#              vectors, each of which contains a control and a treatment 
#              sample name. The list names are created from the two sample
#              names.
# -----------------------------------------------------------------------------
create_contrast_pairs <- function(design_table)
{
  # create the contrast pairs
  contrast_pairs = list()
  for (TREAT_name in design_table[(design_table$contrast != 'CNTRL'),'name']) 
  {
    if (design_table[TREAT_name,'contrast'] %in% design_table$sample_id) 
    {
      CNTRL_name = design_table[(design_table$sample_id == design_table[TREAT_name,'contrast']),'name']
      contrast_pairs[[paste(TREAT_name, 'vs', CNTRL_name, sep="_")]] = c(TREAT_name, CNTRL_name)
    } else {
      stop(paste("[ERROR] Control sample", 
                 design_table[TREAT_name,'contrast'], "specified as contrast for sample", 
                 design_table[TREAT_name,'sample_id'], "could not be found in the design table. Sample skipped from analysis!", sep=" "))
    }
  }
  return(contrast_pairs)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: write_incl_rownames
# ARGUMENTS: data
#            col_name
#            filename
# DESCRIPTION: Write a dataframe or matrix including its rownames.
# -----------------------------------------------------------------------------
write_incl_rownames <- function(data, col_name, filename, verbose=FALSE)
{
  # merge the rownames to the matrix
  data.incl_rownames =
    cbind(rownames(data),
          data)
  
  # create a column name
  colnames(data.incl_rownames)[1] = col_name
  
  # create the full filename
  full_filename = paste(filename, '.tsv', sep='')
  
  # write the data to the file
  write.table(data.incl_rownames, file=full_filename,
              quote=F, row.names=F, col.names=T, sep="\t")
  
  # give some feedback to the user
  if (verbose) {
    message(paste('[INFO] Wrote successfully: ', full_filename, '\n', sep=''))
  }
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: correct_by_sitecounts_expected_per_chance
# ARGUMENTS: sitecount_matrix.raw_counts
#            considered_region_length
#            nt_frequency_vector
# DESCRIPTION: Removes for each motiv the number of counts that are expected
#              to be found within a region of length 
#              'considered_region_length' given the nucleotide frequencies 
#              'nt_frequency_vector'.
#              Remove motif counts expected per chance and ensure that depleted 
#              k-mers are set to 0 counts. (Hint: it would not make sense to set them to
#              a minus value, since after centering the sitecounts per exon, they would
#              again have more or less the same distribution. That is, the k-mer counts
#              would still count, even though each poly(A) site on an exon was depleted in 
#              the motif. For instance, if we find TGT at one poly(A) site but in no other
#              one, then TGT is depleted for ALL the poly(A) sites, since we expect to 
#              find a TGT all 4^3=64 nucleotides and we expect to find 1-2 TGTs per chance
#              for each poly(A) site).
# -----------------------------------------------------------------------------
correct_by_sitecounts_expected_per_chance <- function(sitecount_matrix.raw_counts, 
                                                      considered_region_length,
                                                      nt_frequency_vector)
{
  # create a new matrix that we can correct
  sitecount_matrix.background_corrected = 
    as.data.frame(sitecount_matrix.raw_counts,
                  stringsAsFactors=FALSE)
  ##kmer = "TCTC"
  for (kmer in colnames(sitecount_matrix.raw_counts))
  {
    # calculate the probability to see the current k-mer
    kmer_nt_probs = numeric()
    for (nt in strsplit(kmer, "")[[1]]) {
      kmer_nt_probs = c(kmer_nt_probs,nt_frequency_vector[nt])
    }
    prob_to_see_kmer = prod(kmer_nt_probs)
    
    # TGTA: (0.3006*0.1996*0.3006*0.3004)*97 = 0.5255453
    nr_kmers_expected_in_region = prob_to_see_kmer*(considered_region_length-nchar(kmer)+1)
    
    # remove the expected counts from each region in the sitecount matrix
    sitecount_matrix.background_corrected[,kmer] = 
      (sitecount_matrix.background_corrected[,kmer] - nr_kmers_expected_in_region)
  }
  
  # set negative counts to zero and round the sitecounts
  sitecount_matrix.background_corrected[(sitecount_matrix.background_corrected < 0)] = 0
  
  # return the background corrected counts
  return(sitecount_matrix.background_corrected)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: center_matrix_rows_by_groups
# ARGUMENTS: matrix_to_center
#            rowname_to_group_mappings
#            group_col
#            scale_by_sd
# DESCRIPTION: centers the rows of a matrix ('matrix_to_center') according
#              to the group to which they belong to as specified by 
#              another matrix ('rowname_to_group_mappings') that has the 
#              same rownames as 'matrix_to_center' and the group to which 
#              the rownames belong to specified in a column with the 
#              colum nname specified by 'group_col'. The results can be
#              scaled by the standard deviation if wanted ('scale_by_sd').
# -----------------------------------------------------------------------------
center_matrix_rows_by_groups <- function(matrix_to_center, 
                                         rowname_to_group_mappings,
                                         group_col,
                                         scale_by_sd=FALSE)
{  
  # map the rownames to the groups
  matrix_to_center.incl_groups = 
    merge(x=matrix_to_center, 
          y=rowname_to_group_mappings[,group_col,drop=F],
          by.x="row.names",
          by.y="row.names",
          all=FALSE)
  rownames(matrix_to_center.incl_groups) = matrix_to_center.incl_groups$Row.names
  matrix_to_center.incl_groups$Row.names <- NULL
  
  # convert the factor back to character
  factor_cols = sapply(matrix_to_center.incl_groups, is.factor)
  matrix_to_center.incl_groups[factor_cols] = 
    lapply(matrix_to_center.incl_groups[factor_cols], as.character)
  
  # center all read data by group
  kmers.aggregated = 
    do.call(rbind,
            lapply(
              split(matrix_to_center.incl_groups, matrix_to_center.incl_groups[,group_col], drop=F),
              function(group_df) { return(scale(group_df[,!(colnames(group_df) %in% group_col), drop=F], center = TRUE, scale=scale_by_sd)) } ))
  
  # return the group centered matrix
  return(kmers.aggregated)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: fit_simple_linear_model_for_each_motif
# ARGUMENTS: contrast_pairs
#            treat_idx
#            ctrl_idx
#            sitecount_matrix.group_centered
#            expr_table
#            PAS_overlap_col
#            polyA2exon_mapping
#            exon_col
#            center_data = TRUE
#            compare_kmer_results_to_full_model = FALSE
#            randomize_data = FALSE
#            verbose = FALSE
# DESCRIPTION: fits a linear model to the given data
# -----------------------------------------------------------------------------
fit_simple_linear_model_for_each_motif <- function(contrast_pairs,
                                                   treat_idx,
                                                   ctrl_idx,
                                                   sitecount_matrix.group_centered,
                                                   expr_table,
                                                   PAS_overlap_col,
                                                   polyA2exon_mapping,
                                                   exon_col,
                                                   results_dir,
                                                   create_files_for_each_motif=FALSE,
                                                   create_plots_for_each_motif=FALSE,
                                                   center_data = TRUE,
                                                   randomize_data = FALSE,
                                                   verbose = FALSE)
{
  # ---------------------------------------------------------------------------
  # randomize the expression table in case wanted
  if (randomize_data == TRUE)
  {
    # randomize the rownames
    rows.original = rownames(expr_table)
    rows.randomized = sample(rows.original)
    
    # get the rows in a very different order
    expr_table.for_model_fitting = expr_table[rows.randomized,]
    
    # however, give back the original names
    rownames(expr_table.for_model_fitting) = rows.original
    
    # give also back the original overlapping info
    # (so that we conserve how many exons will be filterd out
    #  and the only thing that changed are the expression values)
    expr_table.for_model_fitting[rows.original,PAS_overlap_col] = 
      expr_table[rows.original,PAS_overlap_col]
    
  } else {
    # use the expression table as it is (no randomization)
    expr_table.for_model_fitting = expr_table
  }
  
  # ---------------------------------------------------------------------------
  # get the expression table
  expr_table.centered.filtered = 
    create_relative_usage_table(expr_table = expr_table.for_model_fitting,
                                PAS_overlap_col = PAS_overlap_col,
                                polyA2exon_mapping = polyA2exon_mapping,
                                exon_col = exon_col)
  
  # ---------------------------------------------------------------------------
  # get only those entries for which we also have expression values 
  # and ensure that the matrices match
  expr_table.contrast_pairs = 
    expr_table.centered.filtered[,(colnames(expr_table.centered.filtered) %in% unlist(contrast_pairs))]
  
  # ---------------------------------------------------------------------------
  # get the sitecounts for the poly(A) sites we consider
  sitecount_matrix.incl_zero_cols = 
    sitecount_matrix.group_centered[rownames(expr_table.contrast_pairs),]
  
  # ---------------------------------------------------------------------------
  # filter out k-mers having 0 counts (since we cannot fit a model for such a
  # k-mer)
  sitecount_matrix = 
    sitecount_matrix.incl_zero_cols[,apply(sitecount_matrix.incl_zero_cols, 2, 
                                           function(x) {
                                             !(all(x == 0.0))
                                           })]
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Prepare the matrixes
  # ---------------------------------------------------------------------------
  # center the columns of the sitecount matrix
  Ns = sitecount_matrix # should already be centered
  
  # column and row center the expression matrix
  Es = expr_table.contrast_pairs # should already be centered
  
  # center the rows, if wanted
  if (center_data) {
    Es = center.rows(Es)
    center = 'centered'
  } else {
    center = 'NOT_centered'
  }
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Prepare the matrix for the end results
  # ---------------------------------------------------------------------------
  # get the k-mers we are interested in
  kmers = colnames(Ns)
  
  # create a matrix for the final results
  result_cols = c('z_score', 
                  'mean_diff_zscores', 
                  'mean_diff_activity', 
                  'mean_diff_delta')
  results.single_kmer_per_run =
    matrix(ncol=length(result_cols),
           nrow=length(kmers))
  rownames(results.single_kmer_per_run) = kmers
  colnames(results.single_kmer_per_run) = result_cols
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Get the fit for each k-mer
  # ---------------------------------------------------------------------------
  # kmer = colnames(Ns)[1]
  for (kmer in colnames(Ns))
  {
    if (verbose) {
      message(paste("[INFO] FITTING LINEAR MODEL FOR MOTIF: ", kmer, sep=""))
    }
    
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # fit the linear model to the current k-mer
    results = fit_linear_model(N=Ns[,kmer,drop=F],E=Es,lambda=0.0)
    
    # calculate activity differences as well as the corresponding deltas and 
    # z-scores
    diff_results = calculate_activity_differences(activities=t(results$Ahat),
                                                  deltas=t(results$AhatSE),
                                                  contrast_pairs=contrast_pairs,
                                                  treat_idx=treat_idx,
                                                  ctrl_idx=ctrl_idx)
    
    # _________________________________________________________________________
    # -------------------------------------------------------------------------   
    # Store the results
    # -------------------------------------------------------------------------   
    # combined z-score (over all samples, contrast independent)
    results.single_kmer_per_run[kmer,'z_score'] = results$combined.Zscore
    
    # contrast dependent results
    results.single_kmer_per_run[kmer,'mean_diff_zscores'] = 
      diff_results$activity_difference_zscores[kmer,'mean_diff_zscore']
    
    results.single_kmer_per_run[kmer,'mean_diff_activity'] = 
      mean(diff_results$activity_differences[kmer,])
    
    results.single_kmer_per_run[kmer,'mean_diff_delta'] = 
      mean(diff_results$activity_difference_deltas[kmer,])
    
    # ___________________________________________________________________________
    # ---------------------------------------------------------------------------  
    # Create files and/or plots if wanted
    # ---------------------------------------------------------------------------  
    if (create_files_for_each_motif | create_plots_for_each_motif)
    {
      # create the results directory for the current k-mer
      path_kmer_results_dir = paste(results_dir, kmer, sep="/")
      dir.create(path_kmer_results_dir, showWarnings=FALSE, recursive=TRUE)
      
      # Create FILES, if wanted
      if (create_files_for_each_motif)
      {
        # write out the zscores
        write_incl_rownames(data = results[["combined.Zscore"]],
                            col_name = 'motif', 
                            filename = paste(path_kmer_results_dir, 
                                             'combined_zscore', 
                                             sep='/'))
        
        # write out the activities
        write_incl_rownames(data = t(results[["Ahat"]]),
                            col_name = 'sample_id', 
                            filename = paste(path_kmer_results_dir, 
                                             'activities', 
                                             sep='/'))
        
        # write out the activities
        write_incl_rownames(data = t(results[["AhatSE"]]),
                            col_name = 'sample_id', 
                            filename = paste(path_kmer_results_dir, 
                                             'deltas', 
                                             sep='/'))
      }
      
      # Create PLOTS, if wanted
      if (create_plots_for_each_motif)
      {
        # create activity difference plots
        activity_profile_plots(motifs=kmer,
                               activities.mat=t(diff_results$activity_differences),
                               deltas.mat=t(diff_results$activity_difference_deltas),
                               results_dir=path_kmer_results_dir,
                               mar=c(15,4,4,2),
                               mgp=c(3,1,0),
                               colors=NULL,
                               x_labels.vec=NULL,
                               y_label=paste("activity difference"),
                               cex_value=1.0,
                               cex_axis=1.0,
                               cex_lab=1.0,
                               cex_main=1.0,
                               lwd=1.0,
                               pdf_height=6,
                               pdf_width=1+length(contrast_pairs),
                               main.text=paste("activity difference profile\nof motif ", kmer, sep=""),
                               dot_type=1,
                               additional_code_to_be_executed=NULL)
      }
    }
  }
  
  # return the result
  results_list = list()
  results_list[["kapac_results_matrix"]] = results.single_kmer_per_run
  results_list[["expression_table_not_row_centered"]] = expr_table.contrast_pairs
  results_list[["Es"]] = Es
  results_list[["Ns"]] = Ns
  
  # return the results
  return(results_list)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: create_realative_usage_table
# ARGUMENTS: expr_table
#            PAS_overlap_col
#            polyA2exon_mapping
#            exon_col
# DESCRIPTION: creates the relative usage table
# -----------------------------------------------------------------------------
create_relative_usage_table <- function(expr_table,
                                        PAS_overlap_col,
                                        polyA2exon_mapping,
                                        exon_col)
{
  # remove all exons that are overlapping with other PAS
  expr_table.non_overlapping = 
    keep_only_non_overlapping_multipas_exons(expr_table=expr_table,
                                             polyA2exon_mapping=polyA2exon_mapping,
                                             exon_col=exon_col,
                                             PAS_overlap_col=PAS_overlap_col)
  
  # add a pseudocounts
  expr_table.incl_pseudocount = add_pseudocount(expr_table=expr_table.non_overlapping,
                                                cols_to_remove=PAS_overlap_col,
                                                pseudocount=1)
  # calculate the relative usage
  rel_usage_table = calculate_relative_usage(expr_table=expr_table.incl_pseudocount,
                                             cols_to_remove=NULL,
                                             polyA2exon_mapping=polyA2exon_mapping,
                                             exon_col=exon_col,
                                             in_percent=TRUE)
  
  # add a pseudocount and then go into log space
  expr_table.log2 = get_log_expression(expr_table=rel_usage_table,
                                       cols_to_remove=NULL,
                                       pseudocount=0.0)
  
  # center the matrix columns by groups (=exons in our case)
  expr_table.centered = 
    center_matrix_rows_by_groups(matrix_to_center=expr_table.log2, 
                                 rowname_to_group_mappings=polyA2exon_mapping,
                                 group_col=exon_col)
  
  # return the expression table
  return(expr_table.centered)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: add_pseudocount
# ARGUMENTS: expr_table
#            cols_to_remove=NULL
#            pseudocount=1
# DESCRIPTION: adds a pseudocount to a dataframe or matrix within the 
#              columns that are not in 'cols_to_remove'.
# -----------------------------------------------------------------------------
add_pseudocount <- function(expr_table,
                            cols_to_remove=NULL,
                            pseudocount=1)
{
  # add the pseudocount to the columns we are interested in
  expr_table.incl_pseudo_count = 
    expr_table[,!(colnames(expr_table) %in% cols_to_remove)] + pseudocount
  
  # return the expression table
  return(expr_table.incl_pseudo_count)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_relative_usage
# ARGUMENTS: expr_table
#            cols_to_remove
#            polyA2exon_mapping
#            exon_col
#            in_percent=TRUE
# DESCRIPTION: Given a poly(A) site expression table, the function calculates
#              the relative usage for all columns except 'cols_to_remove' 
#              in percent (if in_percent=TRUE) or as fractions otherwise.
#              The usage is relative to other poly(A) sites located on the 
#              same exon as specified by the matrix 'polyA2exon_mapping', 
#              which contains poly(A) site ids as rownames and associated
#              exon ids in the 'exon_col' column.
# -----------------------------------------------------------------------------
calculate_relative_usage  <- function(expr_table,
                                      cols_to_remove,
                                      polyA2exon_mapping,
                                      exon_col,
                                      in_percent=TRUE)
{
  # map the rownames to the groups
  expr_table.incl_groups = 
    merge(x=expr_table[,!(colnames(expr_table) %in% cols_to_remove)], 
          y=polyA2exon_mapping[,exon_col,drop=F],
          by.x="row.names",
          by.y="row.names",
          all=FALSE)
  rownames(expr_table.incl_groups) = expr_table.incl_groups$Row.names
  expr_table.incl_groups$Row.names <- NULL
  
  # convert the factor back to character
  factor_cols = sapply(expr_table.incl_groups, is.factor)
  expr_table.incl_groups[factor_cols] = 
    lapply(expr_table.incl_groups[factor_cols], as.character)
  
  # center all read data by group
  relative_usage = 
    do.call(rbind,
            lapply(
              split(expr_table.incl_groups, expr_table.incl_groups[,exon_col], drop=F),
              function(group_df) { 
                return(as.matrix(sweep(group_df[,!(colnames(group_df) %in% exon_col)], 
                                       2, 
                                       colSums(group_df[,!(colnames(group_df) %in% exon_col)]), FUN='/'))) } ))
  
  # calculate in percent if wanted
  if (in_percent)
  {
    relative_usage = relative_usage * 100
  }
  
  # return the expression table
  return(relative_usage)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: get_log_expression
# ARGUMENTS: expr_table
#            cols_to_remove=NULL
#            pseudocount=1
# DESCRIPTION: Get the log2 expression of theexpression table except the
#              'cols_to_remove'. Add a pseudocount before going to log-space.
# -----------------------------------------------------------------------------
get_log_expression <- function(expr_table,
                               cols_to_remove=NULL,
                               pseudocount=1)
{
  # add the pseudocount to the columns we are interested in
  expr_table.incl_pseudo_count = 
    expr_table[,!(colnames(expr_table) %in% cols_to_remove)] + pseudocount
  
  # get the log2 values
  expr_table.log2 = log2(expr_table.incl_pseudo_count)
  
  # return the expression table
  return(expr_table.log2)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: get_log_expression
# ARGUMENTS: expr_table
#            cols_to_remove=NULL
#            pseudocount=1
# DESCRIPTION: Get the log2 expression of theexpression table except the
#              'cols_to_remove'. Add a pseudocount before going to log-space.
# -----------------------------------------------------------------------------
keep_only_non_overlapping_multipas_exons  <- function(expr_table,
                                                      polyA2exon_mapping,
                                                      exon_col,
                                                      PAS_overlap_col)
{
  expr_table.incl_groups = 
    merge(x=expr_table, 
          y=polyA2exon_mapping[,exon_col,drop=F],
          by.x="row.names",
          by.y="row.names",
          all=FALSE)
  rownames(expr_table.incl_groups) = expr_table.incl_groups$Row.names
  factor_cols = sapply(expr_table.incl_groups, is.factor)
  expr_table.incl_groups[factor_cols] = 
    lapply(expr_table.incl_groups[factor_cols], as.character)
  expr_table.incl_groups.filtered = 
    do.call(rbind,
            lapply(
              split(expr_table.incl_groups, expr_table.incl_groups[,exon_col], drop=F),
              function(group_df) { 
                if (dim(group_df)[1] > 1) {
                  if (length(unique(group_df[,PAS_overlap_col])) == 1) {
                    if (unique(group_df[,PAS_overlap_col]) == "OK") {
                      return(group_df)
                    }
                  }
                }
              }))
  rownames(expr_table.incl_groups.filtered) = expr_table.incl_groups.filtered$Row.names
  expr_table.incl_groups.filtered$Row.names <- NULL
  return(expr_table.incl_groups.filtered[,!(colnames(expr_table.incl_groups.filtered) %in% exon_col)])
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: center.rows
# ARGUMENTS: X
# DESCRIPTION: Center the rows of a matrix X.
# -----------------------------------------------------------------------------
center.rows <- function(X) 
{
  return( t( scale( t(X), scale=FALSE) ) )
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: fast.svd
# ARGUMENTS: M
# DESCRIPTION: Perform an SVD (fast) on matrix M and return the result.
# -----------------------------------------------------------------------------
perform_svd_fast <- function(M) 
{
  if (nrow(M) > 2*ncol(M)) {
    s = svd(crossprod(M)) 
    s$d = sqrt(s$d)
    s$u = M %*% sweep(s$v,2,s$d,FUN='/') 
  } else if (ncol(M) > 2*nrow(M)) {
    s = svd(tcrossprod(M)) 
    s$d = sqrt(s$d)
    s$v = sweep(crossprod(M,s$u),2,s$d,FUN='/') 
  } else {
    s = svd(M)
  }
  
  return(s)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: fit_linear_model
# ARGUMENTS: N
#            E
#            lambda
# DESCRIPTION: Fit a linear model and return a list
#              containing a matrix with coefficients ('Ahat') 
#              and errors ('AhatSE').
# -----------------------------------------------------------------------------
fit_linear_model <- function(N,E,lambda)
{
  
  # check whether the matrixes fit
  if ( (!identical(rownames(N),rownames(E))) 
       | (is.null(rownames(N))) 
       | (is.null(rownames(E)))) {
    stop('Rownames of N and E matrixes are not equal. Please ensure that the matrices match!')
  }
  
  # do a SVD on matrix N
  Ns = perform_svd_fast(N)
  
  # calculate the right hand side
  rhs = crossprod(Ns$u,E)
  
  # calculate diagonal matrix
  dia = Ns$d/(Ns$d^2 + nrow(N)*lambda)
  
  # calculate the activities
  Ahat = sweep(Ns$v,2,dia,FUN='*') %*% rhs
  dimnames(Ahat) = list(colnames(N),colnames(E)) 
  
  # calculate Chi2
  Chi2 = colSums((E - N %*% Ahat)^2)
  
  # calculate matrix C
  C = tcrossprod(sweep(Ns$v,2,1/(Ns$d^2 + nrow(N)*lambda),FUN='*'),Ns$v)
  
  # calculate the errors
  AhatSE = sqrt(diag(C) %x% t(Chi2/nrow(E)))
  colnames(AhatSE) = colnames(Ahat)
  rownames(AhatSE) = rownames(Ahat)
  
  # calculate the z-scores
  Zscore = Ahat/AhatSE
  
  # calculate the combined z-score
  combined.Zscore = sqrt(rowMeans(Zscore^2))
  
  # return the results in form of a list
  fit = list(Ahat=Ahat,
             Zscore=Zscore,
             combined.Zscore=combined.Zscore,
             AhatSE=AhatSE,
             Chi2=Chi2)
  return(fit)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_activity_differences
# ARGUMENTS: activities
#            deltas
#            contrast_pairs
#            treat_idx
#            ctrl_idx
# DESCRIPTION: Calculate activity differences, errors and z-scores.
# -----------------------------------------------------------------------------
calculate_activity_differences <- function(activities,
                                           deltas,
                                           contrast_pairs,
                                           treat_idx,
                                           ctrl_idx)
{  
  # specify the order of kmers we want to work with
  kmers = colnames(activities)
  
  # create matrices for the activities and the deltas
  activities.mat = matrix(nrow=length(names(contrast_pairs)),ncol=length(kmers),
                          dimnames=list(names(contrast_pairs), kmers))
  deltas.mat = matrix(nrow=length(names(contrast_pairs)),ncol=length(kmers),
                      dimnames=list(names(contrast_pairs), kmers))
  zscores.mat = matrix(nrow=length(names(contrast_pairs)),ncol=length(kmers),
                       dimnames=list(names(contrast_pairs), kmers))
  
  # calculate activity differences for the current contrast
  for (contrast in names(contrast_pairs)) 
  {
    activities.mat[contrast, kmers] = 
      calculate_activity_difference(activity_ctrl=activities[contrast_pairs[[contrast]][ctrl_idx],kmers,drop=F],
                                    activity_treat=activities[contrast_pairs[[contrast]][treat_idx],kmers,drop=F])
    deltas.mat[contrast, kmers] = 
      calculate_delta_difference(delta_ctrl=deltas[contrast_pairs[[contrast]][ctrl_idx],kmers,drop=F],
                                 delta_treat=deltas[contrast_pairs[[contrast]][treat_idx],kmers,drop=F])
    zscores.mat[contrast, kmers] = 
      activities.mat[contrast, kmers] / deltas.mat[contrast, kmers]
  }
  
  # calculate the mean z-score
  zscores.mat = 
    rbind(apply(zscores.mat, 2, mean), zscores.mat)
  rownames(zscores.mat)[1] = c("mean_diff_zscore")
  
  # finally create a zorted data frame from it
  mean_diff_zscore_df = as.data.frame(t(zscores.mat))
  mean_diff_zscore_df.sorted = mean_diff_zscore_df[order(abs(mean_diff_zscore_df$mean_diff_zscore), decreasing=T), ]
  
  # write all the results into one list  
  results_list = list()
  results_list[["activity_differences"]] = t(activities.mat)
  results_list[["activity_difference_deltas"]] = t(deltas.mat)
  results_list[["activity_difference_zscores"]] = t(zscores.mat)
  results_list[["activity_difference_zscores.sorted"]] = mean_diff_zscore_df.sorted
  
  # return the results
  return(results=results_list)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_activity_difference
# ARGUMENTS: activity_ctrl
#            activity_treat
# DESCRIPTION: Calculate the activity difference of two given activities.
# -----------------------------------------------------------------------------
calculate_activity_difference <- function(activity_ctrl,
                                          activity_treat)
{
  activity_difference = activity_treat - activity_ctrl
  return(activity_difference)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_delta_difference
# ARGUMENTS: delta_ctrl
#            delta_treat
# DESCRIPTION: Calculate the delta difference of two given deltas.
# -----------------------------------------------------------------------------
calculate_delta_difference <- function(delta_ctrl,
                                       delta_treat)
{ 
  delta_difference = sqrt(delta_treat^2.0 + delta_ctrl^2.0)
  return(delta_difference)
} 

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calcultate_zvalue
# ARGUMENTS: coefficient
#            error
# DESCRIPTION: Calculate a z-value for a given coefficient and error.
# -----------------------------------------------------------------------------
calcultate_zvalue <- function(coefficient,
                              error)
{ 
  zvalue = coefficient / error
  return(zvalue)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: activity_profile_plots
# ARGUMENTS: motifs
#            activities.mat
#            deltas.mat
#            results_dir=NULL
#            mar=c(10,10,10,10)
#            mgp=c(6,2,0)
#            colors=NULL
#            x_labels.vec=NULL
#            y_label=NULL
#            cex_value=1.0
#            cex_axis=1.0
#            cex_lab=1.0
#            cex_main=1.0
#            lwd=1.0
#            pdf_height=10
#            pdf_width=10
#            main.text=''
#            dot_type=19
#            additional_code_to_be_executed=NULL
# DESCRIPTION: Create a plot for multiple motifs.
# -----------------------------------------------------------------------------
activity_profile_plots <- function(motifs,
                                   activities.mat,
                                   deltas.mat,
                                   results_dir=NULL,
                                   mar=c(10,10,10,10),
                                   mgp=c(6,2,0),
                                   colors=NULL,
                                   x_labels.vec=NULL,
                                   y_label=NULL,
                                   cex_value=1.0,
                                   cex_axis=1.0,
                                   cex_lab=1.0,
                                   cex_main=1.0,
                                   lwd=1.0,
                                   pdf_height=10,
                                   pdf_width=10,
                                   main.text='',
                                   dot_type=19,
                                   additional_code_to_be_executed=NULL)
{
  # ---------------------------------------------------------------------------
  # create one plot per motif
  for (motif in motifs)
  {
    # create the plots
    activity_profile_plot(motifs_to_plot.vec=motif,
                          activities.mat=activities.mat,
                          deltas.mat=deltas.mat,
                          results_dir=results_dir,
                          file_name=paste("activity_difference_profile", sep="_"),
                          mar=mar,
                          mgp=mgp,
                          colors=c("blue"),
                          x_labels.vec=NULL,
                          y_label=y_label,
                          cex_value=cex_value,
                          cex_axis=cex_axis,
                          cex_lab=cex_lab,
                          cex_main=cex_main,
                          lwd=lwd,
                          pdf_height=pdf_height,
                          pdf_width=pdf_width,
                          main.text=main.text,
                          dot_type=dot_type,
                          additional_code_to_be_executed=NULL)        
  }
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: activity_profile_plot
# ARGUMENTS: motifs_to_plot.vec
#            activities.mat
#            deltas.mat
#            results_dir=NULL
#            file_name=NULL
#            mar=c(10,10,10,10)
#            mgp=c(6,2,0)
#            colors=NULL
#            x_labels.vec=NULL
#            y_label=NULL
#            cex_value=1.0
#            cex_axis=1.0
#            cex_lab=1.0
#            cex_main=1.0
#            lwd=1.0
#            pdf_height=10
#            pdf_width=10
#            main.text=''
#            dot_type=1
#            additional_code_to_be_executed=NULL
# DESCRIPTION: Create a motif activity profile plot.
# -----------------------------------------------------------------------------
activity_profile_plot <- function(motifs_to_plot.vec,
                                  activities.mat,
                                  deltas.mat,
                                  results_dir=NULL,
                                  file_name=NULL,
                                  mar=c(10,10,10,10),
                                  mgp=c(6,2,0),
                                  colors=NULL,
                                  x_labels.vec=NULL,
                                  y_label=NULL,
                                  cex_value=1.0,
                                  cex_axis=1.0,
                                  cex_lab=1.0,
                                  cex_main=1.0,
                                  lwd=1.0,
                                  pdf_height=10,
                                  pdf_width=10,
                                  main.text='',
                                  dot_type=1,
                                  additional_code_to_be_executed=NULL)
{
  # ---------------------------------------------------------------------------
  # check if the requirements that we need to do the plots are fulfilled
  # (we cannot put more than two activity profiles into one plot)
  if (length(motifs_to_plot.vec) > 2) {
    stop(paste("This function cannot put more than two motifs into one plot! ",
               "Please choose < 3 motifs!", sep=''))
  }
  
  # ---------------------------------------------------------------------------
  # create some colors, if we did not specify them
  if (is.null(colors))
    colors = rainbow(length(motifs_to_plot.vec))
  
  # ---------------------------------------------------------------------------
  # create the plotting dirctory if it does not exist already
  # if there is given a results dir, we write the plot into
  # a file
  if (!is.null(file_name))
  {
    if (!is.null(results_dir))
    {
      # create the directory if it does not exist
      dir.create(results_dir, 
                 showWarnings = FALSE)
    }
    # split the filename
    file_name.split = unlist(strsplit(x=file_name, split='.', fixed=TRUE))
    # create the file
    if ((length(file_name.split) > 1) & (tail(file_name.split, n=1) == "png")) {
      png(paste(results_dir, '/', file_name, sep=''), 
          height=pdf_height, width=pdf_width, res=300, units="in") # we use inches, because pdf does this also per default
    } else if ((length(file_name.split) > 1) & (tail(file_name.split, n=1) == "pdf")) {
      pdf(paste(results_dir, '/', file_name, sep=''), 
          height=pdf_height, width=pdf_width, useDingbats=F)
    } else {
      pdf(paste(results_dir, '/', file_name, '.pdf', sep=''), 
          height=pdf_height, width=pdf_width, useDingbats=F)
    }
  }
  
  # ---------------------------------------------------------------------------
  # check if there are lables given
  if (is.null(x_labels.vec)) {
    labels = rownames(activities.mat)
  } else {
    if (length(x_labels.vec) == nrow(activities.mat))
    {
      labels = x_labels.vec
    } else {
      stop('[ERROR] given lables vector length differs from activities matrix length!\n')
    }
  }
  
  # create the x_axis datapoints
  x_axis_datapoints = seq(from=1,to=length(labels),by=1)
  
  # ---------------------------------------------------------------------------
  # I think this was done by Piotr
  profileHeight = 960;
  if( max(x_axis_datapoints) > 50 )
  {
    profileHeight = 120 + max(x_axis_datapoints) * 15;
  }
  width = 1000;
  width = width/100;
  profileHeight = profileHeight/100;
  
  # ---------------------------------------------------------------------------
  # define the x-axis boarders
  xmin <- min(x_axis_datapoints)
  xmax <- max(x_axis_datapoints) #*1.05
  
  # ---------------------------------------------------------------------------
  # create the plot for the first motif
  for (motif_idx in 1:min(c(length(motifs_to_plot.vec),2)))
  {
    motif = motifs_to_plot.vec[motif_idx]
    y_axis_activities = activities.mat[,motif]
    y_axis_deltas = deltas.mat[,motif]
    
    # -------------------------------------------------------------------------
    # define the y-axis boarders
    scaling_factor = 1000.0
    ymax = ceiling(max(activities.mat[,motif] + deltas.mat[,motif])*scaling_factor) / scaling_factor
    ymin = 0.0 - ymax
    if (ymin > min(activities.mat[,motif] - deltas.mat[,motif]))
    {
      ymin = floor(min(activities.mat[,motif] - deltas.mat[,motif])*scaling_factor) / scaling_factor
      ymax = 0.0 - ymin
    }  
    
    # -------------------------------------------------------------------------
    # create the plot
    if (motif_idx == 1) {
      par(mar=mar, mgp=mgp)
      side=2
      padj=1
    } else {
      par(new=T)
      side=4
      padj=-1
    }
    plot(x_axis_datapoints, 
         y_axis_activities,
         xlim = c(xmin, xmax),
         ylim = c(ymin, ymax),
         type = "o",
         lwd=lwd,
         #cex.axis=cex_axis,
         cex.axis=cex_axis,
         cex.lab=cex_lab,
         cex.main=cex_main,
         main = main.text, 
         col = colors[motif_idx],
         pch = dot_type,
         yaxt="n", ylab="", xaxt="n", xlab="", bty="n",
         axes = F);
    
    # create the error bars
    arrows(x_axis_datapoints, y_axis_activities - y_axis_deltas,
           x_axis_datapoints, y_axis_activities + y_axis_deltas,
           code = 3, angle = 90, length= 0.05,
           col = colors[motif_idx],
           lty = 1, lwd=lwd)
    
    # add the y axis (try also the adj=0.1 param)
    axis(side=side, line = 0, padj=padj, cex.axis=cex_axis)
    
    # add the y-axis text
    if (is.null(y_label)) {
      ylab_text = paste(unlist(strsplit(x=motif, split='.', fixed=TRUE))[1], "activities", sep=' ')
    } else {
      ylab_text = y_label
    }
    mtext(side=side, line = 3, col=colors[motif_idx], cex=cex_lab,
          text = ylab_text)
  }
  
  # ---------------------------------------------------------------------------
  # add the rest to the plot
  
  # add a line at 0
  abline(0, 0, lty=2, col = "black")
  
  # add the x-axis
  axis(1, 
       at=c(xmax, xmin),
       labels=FALSE,
       lwd=1,
       lwd.ticks=0)
  
  # add the tick marks for the x-axis
  axis(1, 
       at=unique(x_axis_datapoints), 
       labels=FALSE,
       tcl=-0.5,        # how long are the ticks lines and in which direction do they go
       lwd.ticks=1)
  
  # add the labels for the x-axis
  axis(1, 
       at=unique(x_axis_datapoints), #c(1.0,0.5,0.0,-0.5,-1.0,-1.5,-2.0),
       lwd=0,
       cex.axis=cex_axis,
       lwd.ticks=0,       
       labels=labels,
       las=2)
  
  # add a box to the plot
  box()
  
  # if we should add some additional staff to the plot, we do that now
  if (!is.null(additional_code_to_be_executed))
  {
    eval(parse(text=additional_code_to_be_executed))
  }
  
  # if we have created a new file, we have also to close it now
  if (!is.null(file_name))
  {
    dev.off();
  }  
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# read in the design table and create contrasts
# -----------------------------------------------------------------------------
# create the design table from the design file
design_table = create_design_table(sample_design_file_path = opt$sample_design)

# create contrast pairs
contrast_pairs = create_contrast_pairs(design_table = design_table)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read in the PAS to exon mapping so that we are able to center by exons
# later on
# -----------------------------------------------------------------------------
polyA2exon_mapping = 
  as.matrix(read.table(opt$pas_exon_associations, 
                       h=TRUE, 
                       as.is=T, 
                       sep="\t",
                       row.names="pas",
                       comment.char="",
                       stringsAsFactors=FALSE,
                       check.names=FALSE))

# get the name of the column with the exons
exon_col="exon"

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read in the expression table
# -----------------------------------------------------------------------------
expr_table = 
  read.table(opt$expression_matrix, 
             h=TRUE, 
             as.is=T, 
             sep="\t",
             row.names="pas",
             comment.char="",
             check.names=FALSE)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read in the full sitecount matrix
# -----------------------------------------------------------------------------
sitecount_matrix.all = 
  as.matrix(read.table(opt$sitecount_matrix, 
                       h=TRUE, 
                       as.is=T, 
                       sep="\t",
                       row.names="pas",
                       comment.char="",
                       check.names=FALSE))

# -----------------------------------------------------------------------------
# Determine how many motifs are considered in total (can be used for 
# Bonferroni correction).
total_nr_considered_kmers = ncol(sitecount_matrix.all)

# -----------------------------------------------------------------------------
# Select specific k-mers, if wanted
if (opt$selected_motifs != "all")
{
  # read in the selected kmers
  selected_kmers = 
    read.table(opt$selected_motifs, 
               h=TRUE, 
               as.is=T, 
               sep="\t",
               comment.char="",
               check.names=FALSE)[,1]
  
  # drop all other k-mers
  sitecount_matrix.all = sitecount_matrix.all[,selected_kmers,drop=F]
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Filter out k-mers that have counts in less than x [%] of all poly(A) sites
# in the genome.
# In case we want to run on selected k-mers, we want to get estimated
# -----------------------------------------------------------------------------
cols_to_keep = 
  apply(sitecount_matrix.all, 2, 
        function(x) { (sum(x > 0) / length(x)) >= opt$min_kmer_abundance_fraction } )

# -----------------------------------------------------------------------------
# get the columns we are interested in
sitecount_matrix.raw_counts = 
  sitecount_matrix.all[, cols_to_keep, drop = FALSE]

# -----------------------------------------------------------------------------
# get the k-mers that are filtered out
filtered_out_kmers = 
  colnames(sitecount_matrix.all[ ,!(colnames(sitecount_matrix.all) %in% colnames(sitecount_matrix.raw_counts))])

# -----------------------------------------------------------------------------
# calculate the fraction of filtered out k-mers
fraction_of_filtered_out_kmers = 
  length(filtered_out_kmers) / ncol(sitecount_matrix.all)

if (opt$verbose) {
  message(paste(docs_80_unerlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] FILTERING OUT K-MERS HAVING NO/LOW ABUNDANCE.", sep=""))
  
  # report how many k-mers will be dropped
  message(paste("[INFO] Minimum k-mer abundance fraction to be considered (--min_kmer_abundance_fraction): ", 
                opt$min_kmer_abundance_fraction, "\n", sep=""))
  
  # report how many k-mers will be dropped
  message(paste("[INFO] Percentage of filtered out k-mers: ", 
                (fraction_of_filtered_out_kmers*100), "[%] (=", length(filtered_out_kmers), " k-mers). ", sep=""))
}

# background correct, if wanted
if (opt$consider_excess_counts_only == TRUE) {
  
  if (opt$verbose) {
    message(paste(docs_80_unerlines, sep=""))
    message(paste(docs_80_dashes, sep=""))
    message(paste("[INFO] BACKGROUND CORRECTION: 'Excess' motif counts will be removed.", sep=""))
    message(paste("[INFO] The following nucleotide frequencies will be used:", sep=""))
    message(paste("[INFO] ", names(nt_frequency_vector), " ", nt_frequency_vector, "\n", sep=" "))
  }
  
  # perform the background correction
  sitecounts_fit = 
    correct_by_sitecounts_expected_per_chance(sitecount_matrix.raw_counts=sitecount_matrix.raw_counts, 
                                              considered_region_length=opt$considered_region_length,
                                              nt_frequency_vector=nt_frequency_vector)
} else {
  
  # give some feedback to the user
  if (opt$verbose) {
    message(paste(docs_80_unerlines, sep=""))
    message(paste(docs_80_dashes, sep=""))
    message(paste("[INFO] NO BACKGROUND CORRECTION will be applied.", sep=""))
    message(paste("[INFO] Raw motif counts will be used.", sep=""))
  }
  
  # in case we do not want to perform background correction, we use the raw counts
  sitecounts_fit = sitecount_matrix.raw_counts
}

# -----------------------------------------------------------------------------
# center the motif counts per exon
# -----------------------------------------------------------------------------
if (opt$verbose) {
  message(paste(docs_80_unerlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] CENTERING MOTIF COUNTS PER EXON.", sep=""))
}

# center the matrix
sitecount_matrix.group_centered = 
  center_matrix_rows_by_groups(matrix_to_center=sitecounts_fit, 
                               rowname_to_group_mappings=polyA2exon_mapping,
                               group_col=exon_col)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Now fit the model to the actual data (not randomized)
# -----------------------------------------------------------------------------

if (opt$verbose) {
  message(paste(docs_80_unerlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] FITTING MOTIF ACTIVITIES.", sep=""))
}

# run kapac on the data
kapac.result = 
  fit_simple_linear_model_for_each_motif(contrast_pairs = contrast_pairs,
                                         treat_idx = treat_idx,
                                         ctrl_idx = ctrl_idx,
                                         sitecount_matrix.group_centered = sitecount_matrix.group_centered, 
                                         expr_table = expr_table,
                                         PAS_overlap_col = pas_overlap_col,
                                         polyA2exon_mapping = polyA2exon_mapping,
                                         exon_col = exon_col,
                                         results_dir = opt$results_dir,
                                         create_files_for_each_motif = opt$create_files_for_each_motif,
                                         create_plots_for_each_motif = opt$create_plots_for_each_motif,
                                         center_data = opt$row_center_expression,
                                         randomize_data = FALSE,
                                         verbose = opt$verbose)

# get the result matrix
kapac.result_matrix = 
  kapac.result[["kapac_results_matrix"]]

# get the rownames, so that we can bind afterwards the random result to it
kapac.result_matrix = as.data.frame(kapac.result_matrix)
kapac.result_matrix.rownames = rownames(kapac.result_matrix)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Fit the model to randomized data (if wanted) 
# -----------------------------------------------------------------------------
if (opt$verbose) {
  message(paste(docs_80_unerlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] FITTING MOTIF ACTIVITIES FOR RANDOMIZED EXPRESSION DATA.", sep=""))
}

# create a vector for the column names of the random runs
kapac.result_matrix.rand_run_cols = character(0)

# if required, we create random runs
if (opt$number_of_randomized_runs < min_zstatistic_sample_size) {

  # in case the number of randomized runs is smaller than the recommended
  # number, drop a warning.
  if (opt$number_of_randomized_runs < recommended_zstatistic_sample_size) {
    warning(paste("[WARNING] The chosen number of randomized runs ",
                  "'--number_of_randomized_runs' (=", 
                  opt$number_of_randomized_runs, ") ",
                  "is smaller than the minimum number required (=",
                  min_zstatistic_sample_size, ") for calculating ",
                  "z-scores on activity differences. ",
                  "Therefore z-scores and p-values can not be provided ",
                  "in the final results. ",
                  sep=""))
  }
  
  # specify the columns that should be written to the results file
  cols_to_write = 
    c("mean_diff_zscores",
      "mean_diff_activity",
      "mean_diff_delta")
  
} else {

  # in case the number of randomized runs is smaller than the recommended
  # number, drop a warning.
  if (opt$number_of_randomized_runs < recommended_zstatistic_sample_size) {
    warning(paste("[WARNING] The chosen number of randomized runs ",
                  "'--number_of_randomized_runs' (=", 
                  opt$number_of_randomized_runs, ") ",
                  "is relatively small and might lead to variable results. ",
                  "The recommended minimum number of randomized runs is ",
                  recommended_zstatistic_sample_size, ".",
                  sep=""))
  }
  
  rand_runs_prefix = "Random_run"
  #rand_run_nr = 1
  for (rand_run_nr in (1:opt$number_of_randomized_runs))
  {
    message(paste("[INFO] Fitting activities for randomized expression data: ", rand_run_nr, sep=""))
    
    # run kapac on randomized data
    kapac.rand_result = 
      fit_simple_linear_model_for_each_motif(contrast_pairs = contrast_pairs,
                                             treat_idx = treat_idx,
                                             ctrl_idx = ctrl_idx,
                                             sitecount_matrix.group_centered = sitecount_matrix.group_centered, 
                                             expr_table = expr_table,
                                             PAS_overlap_col = pas_overlap_col,
                                             polyA2exon_mapping = polyA2exon_mapping,
                                             exon_col = exon_col,
                                             results_dir = NULL,
                                             create_files_for_each_motif = FALSE,
                                             create_plots_for_each_motif = FALSE,
                                             center_data = opt$row_center_expression,
                                             randomize_data = TRUE,
                                             verbose = FALSE)
    
    # get the results matrix  
    kapac.rand_result_matrix = 
      kapac.rand_result[["kapac_results_matrix"]]
    
    # bind the result to the original result
    kapac.result_matrix[kapac.result_matrix.rownames,paste(rand_runs_prefix, rand_run_nr, top_kmer_candidates_selection_criterion, sep="_")] = 
      kapac.rand_result_matrix[kapac.result_matrix.rownames,top_kmer_candidates_selection_criterion, drop=FALSE]
  }
  
  # _____________________________________________________________________________
  # -----------------------------------------------------------------------------
  # Calculate mean, stdev, z-scores for the runs on random data
  # -----------------------------------------------------------------------------
  # get the colnames of the random run results
  kapac.result_matrix.rand_run_cols = 
    grep(x=colnames(kapac.result_matrix), pattern = rand_runs_prefix, value=TRUE)
  
  # calculate the mean of the random runs
  kapac.result_matrix[,paste(rand_runs_prefix, "mean", sep="_")] = 
    apply(kapac.result_matrix[,kapac.result_matrix.rand_run_cols], 1, mean)
  
  # calculate the standard deviation of the random runs
  kapac.result_matrix[,paste(rand_runs_prefix, "stdev", sep="_")] = 
    apply(kapac.result_matrix[,kapac.result_matrix.rand_run_cols], 1, sd)
  
  # calculate the z-score for the top_kmer_candidates_selection_criterion
  kapac.result_matrix[,paste(top_kmer_candidates_selection_criterion, "ZSCORE", sep="_")] = 
    apply(kapac.result_matrix, 1, 
          function(x) {
            diff_from_mean = (x[top_kmer_candidates_selection_criterion] - x[paste(rand_runs_prefix, "mean", sep="_")])
            zscore = diff_from_mean / x[paste(rand_runs_prefix, "stdev", sep="_")]
            return( zscore )
          })
  
  # calculate the pval for the top_kmer_candidates_selection_criterion
  kapac.result_matrix[,paste(top_kmer_candidates_selection_criterion, "PVAL", sep="_")] = 
    apply(kapac.result_matrix, 1, 
          function(x) {
            pval.two_sided = 2 * pnorm(-abs( x[paste(top_kmer_candidates_selection_criterion, "ZSCORE", sep="_")] ))
            return( pval.two_sided )
          })
  
  # calculate the adjusted pval for the "mean_diff_zscores"
  kapac.result_matrix[,paste(top_kmer_candidates_selection_criterion, "PADJ", sep="_")] = 
    p.adjust(p = kapac.result_matrix[,paste(top_kmer_candidates_selection_criterion, "PVAL", sep="_")], 
             method = "bonferroni", 
             n = total_nr_considered_kmers)
  
  # test whether the values are normally distributed
  not_norm_pval_colname = paste("not_normal_PVAL", sep="_")
  kapac.result_matrix[,not_norm_pval_colname] = 
    apply(kapac.result_matrix[,kapac.result_matrix.rand_run_cols], 1, 
          function(x) {
            return(shapiro.test(x)$p.value)
          })
  
  # calculate the adjusted pval for not normal distribution.
  not_norm_padj_colname = paste("not_normal_PADJ", sep="_")
  kapac.result_matrix[,not_norm_padj_colname] = 
    p.adjust(p = kapac.result_matrix[,not_norm_pval_colname], 
             method = "bonferroni", 
             n = total_nr_considered_kmers)
  
  # _____________________________________________________________________________
  # -----------------------------------------------------------------------------
  # Check how many of the z-score distributions are not normally distributed
  # (in this case a z-statistics might not be good to use...)
  # -----------------------------------------------------------------------------
  kapac.result_matrix.not_norm_distributed = 
    kapac.result_matrix[(kapac.result_matrix[,not_norm_padj_colname] < not_norm_pval_threshold),]
  
  # calculate how much percent of the k-mer z-scores are not normal distributed
  percent_of_not_norm_distributed_kmers = 
    (nrow(kapac.result_matrix.not_norm_distributed) / nrow(kapac.result_matrix))*100
  
  # give some user feedback
  if (percent_of_not_norm_distributed_kmers > 0.0)
  {
    warning(paste("[WARNING] ", percent_of_not_norm_distributed_kmers, 
                  " percent of the z-scores are not normal distributed ",
                  " (using a p-value threshold of ", not_norm_pval_threshold, ").", 
                  " Reported z-scores and p-values for activity differences ",
                  " are only valid for motifs (k-mers) that are ",
                  " normally distributed (as indicated by the '",
                  not_norm_pval_colname, "' column in the results).",
                  sep=""))
  }

  # specify the columns that should be written to the results file
  cols_to_write = 
    c("mean_diff_zscores",
      "mean_diff_activity",
      "mean_diff_delta",
      "mean_diff_zscores_PADJ",
      not_norm_padj_colname)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Write out the k-mers sorted by the top_kmer_candidates_selection_criterion
# -----------------------------------------------------------------------------
if (opt$verbose) {
  message(paste(docs_80_unerlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] SORTING RESULTS AND WRITING THEM TO THE RESULTS FILE.", sep=""))
}

# sort it also by zscore
results.single_kmer_per_run.sorted = 
  kapac.result_matrix[(order(abs(kapac.result_matrix[,top_kmer_candidates_selection_criterion]), decreasing=T)),
                              !(colnames(kapac.result_matrix) %in% kapac.result_matrix.rand_run_cols)]

# write out the zscores
write_incl_rownames(data=results.single_kmer_per_run.sorted[, cols_to_write, drop=F],
                    col_name='kmer', 
                    filename=paste(opt$results_dir, '/', top_kmer_candidates_selection_criterion, sep=''))

if (opt$verbose) {
  message(paste(docs_80_unerlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] KAPAC run finished.", sep=""))
}

# ___________________________________________________________________________
# ---------------------------------------------------------------------------
# Since in the sitecount matrices we only consider k-mers that have been found
# in the defined region, it can be that we do not have every k-mer represented
# in the sitecount matrix. This can also happen due to filtering out k-mers 
# based on their abundance (--min_kmer_abundance_fraction). For such k-mers 
# we still might want to have output files. One can create them as shown 
# below.
# ---------------------------------------------------------------------------
# if (opt$create_files_for_each_motif & (opt$selected_motifs != "all"))
# {
#   # get k-mers for which a result is requested, but which was not considered
#   # in the fit (maybe also because there do not exist counts in the currently
#   # investigated region).
#   missing_kmers = setdiff(selected_kmers, colnames(kapac.result[["Ns"]]))
#   
#   if (length(missing_kmers) > 0)
#   {
#     if (opt$verbose) {
#       message(paste(docs_80_unerlines, sep=""))
#       message(paste(docs_80_dashes, sep=""))
#       message(paste("[INFO] WRITING FILES FOR K-MERS THAT ARE REQUESTED (--selected_kmers) BUT DO NOT HAVE COUNTS.", sep=""))
#     }
#     
#     # -----------------------------------------------------------------------------
#     # create a matrix that contains only 0.0 as activities and also as deltas
#     # HINT: if a k-mer is not found at all within a region, it should have an 
#     #       activity of 0.0 and we are also 100% about this (delta=0.0).
#     template_matrix = t(kapac.result$Ahat)
#     template_matrix[,1] = 0.0
#     
#     ##kmer = missing_kmers[1]
#     for (kmer in missing_kmers)
#     {
#       # create the k-mer matrix
#       kmer_matrix = template_matrix
#       colnames(kmer_matrix) = kmer
#       
#       # -----------------------------------------------------------------------------
#       # create the results directory
#       path_kmer_results_dir = paste(opt$results_dir, kmer, sep="/")
#       dir.create(path_kmer_results_dir, showWarnings=FALSE, recursive=TRUE)
#       
#       # prepare a matrix containing the z-scores (sorted)
#       overall_z_scores = as.matrix(c(0.0),nrow=1, ncol=1)
#       rownames(overall_z_scores) = kmer
#       colnames(overall_z_scores) = 'z_score'
#       
#       # -----------------------------------------------------------------------------
#       # write out the zscores
#       write_incl_rownames(data=overall_z_scores,
#                           col_name='pwm', 
#                           filename=paste(path_kmer_results_dir, '/zscore', sep=''))
#       
#       # write out the activities
#       write_incl_rownames(data=kmer_matrix,
#                           col_name='sample_id', 
#                           filename=paste(path_kmer_results_dir, '/activities', sep=''))
#       
#       # write out the activities
#       write_incl_rownames(data=kmer_matrix,
#                           col_name='sample_id', 
#                           filename=paste(path_kmer_results_dir, '/deltas', sep=''))
#     }
#   }
# }

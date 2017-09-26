configfile: "config.yaml"

from snakemake.utils import makedirs
from snakemake.utils import listfiles

import numpy as np
import os
import yaml

localrules: design_file, tpm_normalize, get_filtered_rel_usages, rel_pos_of_pAs, weighted_average_exon_length, plot_average_exon_length, all

################################################################################
# define the target rule
################################################################################

rule all:
    ## LOCAL ##
    '''
    Defines the target rule by specifying all final output files.
    Additionally, the cluster_logs dir is deleted if
    snakemake was run locally.
    '''
    input:
        cdf_plots = expand( "{study}/weighted_avg_exonLength.CDFs.pdf", study = config['studies']),
        order_of_samples = expand( "{study}/INFO.tsv", study = config['studies']),
        expressions = expand( "{study}/tandem_pas_expressions.rpm.filtered.tsv", study = config['studies'])
    run:
        # check if cluster_logs dir for first sample is empty
        # if yes, delete the cluster_logs dir
        for s in config['studies']:
            file_list = os.listdir("{log}/cluster_logs/".format(log = s + "/logs") )
            if len(file_list) == 0:
                shell("rm -r {s}/logs/cluster_logs")
            break

################################################################################
# chronological processing of the reads
################################################################################

#-------------------------------------------------------------------------------
# get the terminal exon read coverage  profiles
#-------------------------------------------------------------------------------
rule create_coverages:
    input:
        bam = "{study}/valid_samples/{sample_na}.bam"
    output:
        "{study}/coverages/{sample_na}.pkl"
    params:
        script_dir = config['dir.scripts'],
        py2_env_path = config['py2_env_path'],
        clusters= config['clusters'],
        ds_extend=config['cvg.ds_extend'],
        start2prox_minDist=config['cvg.start2prox_minDist'],
        unstranded = lambda wildcards: "--unstranded" if config['cvg.unstranded'] == "yes" else "",
        cluster_log = "{study}/logs/cluster_logs/{sample_na}.log"
    resources:
        mem = lambda wildcards: int( np.ceil(50 / config['cvg.threads'] ) ),
        time = 6
    threads:
        config['cvg.threads']
    log:
        "{study}/logs/{sample_na}/create_coverages.log"
    shell:
        '''
        {params.py2_env_path}/python {params.script_dir}/mRNAseq-coverage-objects-and-wiggleFiles.py \
        --verbose \
        --bam {input.bam} \
        --cluster {params.clusters} \
        --ds_extension {params.ds_extend} \
        --min_dist_exStart2prox {params.start2prox_minDist} \
        --processors {threads} \
        --pickle_out {output} {params.unstranded} \
        2> {log}
        '''

#-------------------------------------------------------------------------------
# infer used poly(A) sites based on the coverage profiles
# and determine for those the relative usage
#-------------------------------------------------------------------------------

rule infer_relative_usage:
    input:
        coverages = lambda wildcards: expand(wildcards.study + "/coverages/{sample_na}.pkl", sample_na = [line.rstrip().split("\t")[0] for line in open(wildcards.study + "/no_bias_samples.out", "r").readlines()])
    output:
        relUse = "{study}/relative_usages.tsv",
	expressions = "{study}/tandem_pas_expressions.tsv",
        header = "{study}/relative_usages.header.out"
    params:
        script_dir = config['dir.scripts'],
        py2_env_path = config['py2_env_path'],
        clusters = config['clusters'],
        read_len = config['relUse.read_length'],
        min_region_length_to_infer_meanCov = config['relUse.minLength.meanCvg'],
        rawCounts_minMeanCoverage_per_sample = config['relUse.minMeanCvg.perSample'],
        downstream2search4noCov = config['relUse.distal_ds'],
        max_downstream_cov_relative_to_start = config['relUse.distal_ds.maxCvg'],
        min_cluster_distance = config['relUse.min_cluster_distance'],
        valid_region_upstream_of_cluster_to_locate_globalBestBreakPoint_in = config['relUse.us_reg_for_best_breakPoint'],
        mse_ratio_threshold = config['relUse.mse_ratio_threshold'],
        cluster_log = "{study}/logs/cluster_logs/infer_relative_usage.log"
    threads:
        config['relUse.threads']
    resources:
        mem = lambda wildcards: int( np.ceil(32 / config['relUse.threads'] ) ),
        time = 6
    log:
        "{study}/logs/infer_relative_usage.log"
    run:
        conds = []
        extens = []
        with open( output.header, "w" ) as header_out:
            for cvg in input.coverages:
                conds.append( config[ cvg.split("/")[-1].split(".")[0] ]['type'] )
                extens.append( cvg.replace("pkl", "extensions.tsv") )
                header_out.write("%s\n" % cvg.split("/")[-1].split(".")[0])
        conds_string = " ".join(conds)
        extens_string = " ".join(extens)
        shell('''
        {params.py2_env_path}/python {params.script_dir}/deNovo-used-sites-and-usage-inference.py \
        --verbose \
	--expressions_out {output.expressions} \
        --clusters {params.clusters} \
        --coverages {input.coverages} \
        --conditions {conds_string} \
	--ex_extension_files {extens_string} \
        --read_length {params.read_len} \
        --min_coverage_region {params.min_region_length_to_infer_meanCov} \
        --min_mean_coverage {params.rawCounts_minMeanCoverage_per_sample} \
        --ds_reg_for_no_coverage {params.downstream2search4noCov} \
        --min_cluster_distance {params.min_cluster_distance} \
        --mse_ratio_threshold {params.mse_ratio_threshold} \
        --best_break_point_upstream_extension {params.valid_region_upstream_of_cluster_to_locate_globalBestBreakPoint_in} \
        --processors {threads} \
        --max_downstream_coverage {params.max_downstream_cov_relative_to_start} \
        > {output.relUse} \
        2> {log}
        ''')


#-------------------------------------------------------------------------------
# create design file for polyA-MARA
#-------------------------------------------------------------------------------
rule design_file:
    ##LOCAL##
    input:
        header = "{study}/relative_usages.header.out"
    output:
        "{study}/INFO.tsv"
    run:
        with open (output[0], "w") as outfile:
            outfile.write("sample_id\tcell_type\ttreatment\treplicate\tgroup\tcontrast\n")
            study_abbr = "{study}"
            for i in  [line.rstrip() for line in open( input.header, "r").readlines() ] :
                if 'control' not in config[i]:
                    outfile.write("%s\t%s\t%s\tmRNAseq\tCNTRL\tCNTRL\n"
                          % (i,
                             study_abbr,
                             "control"
                          )
                    )
                else:
                    outfile.write("%s\t%s\t%s\tmRNAseq\tKD\t%s\n"
                          % (i,
                             study_abbr,
                             "kd",
                             config[i]['control']
                             )
                    )

#-------------------------------------------------------------------------------
# tpm normalize the expression values by the number of mapped reads
#-------------------------------------------------------------------------------
rule tpm_normalize:
    ##LOCAL##
    input:
        expression = "{study}/tandem_pas_expressions.tsv"
    output:
        tpm_expr = "{study}/tandem_pas_expressions.rpm.tsv"
    run:
        libsizes = []
        with open(input.expression, "r") as ifile:
            for line in ifile:
                F = line.rstrip().split("\t")
                if len(libsizes) == 0:
                    for i in range(10, len(F)):
                        libsizes.append(0)
                for i in range(10, len(F)):
                    if float(F[i]) == -1:
                        continue
                    libsizes[ i - 10 ] += float(F[i])

        with open(output.tpm_expr, "w") as ofile:
            with open(input.expression, "r") as ifile:
                for line in ifile:
                    F = line.rstrip().split("\t")
                    for i in range(10, len(F)):
                        if float(F[i]) == -1:
                            continue
                        F[i] = "{0:.6f}".format( float(F[i]) / libsizes[ i - 10 ] * 1000000 )
                    ofile.write("%s\n" % "\t".join(F))

#-------------------------------------------------------------------------------
# create a version of the tables with tandem PAS expressions that does only
# contain exons for which all poly(A) sites were quantified in all samples (per study)
#-------------------------------------------------------------------------------
rule get_filtered_rel_usages:
    ##LOCAL##
    input:
	rel_use = "{study}/relative_usages.tsv",
        expressions = "{study}/tandem_pas_expressions.rpm.tsv"
    output:
        rel_use_filtered = "{study}/relative_usages.filtered.tsv",
	expressions_filtered = "{study}/tandem_pas_expressions.rpm.filtered.tsv"
    run:
        exons = {}
        with open(input.rel_use, "r") as rel_use_in:
            for line in rel_use_in:
                selected_entries = []
                F = line.rstrip().split("\t")
                if F[8] not in exons:
                    exons[ F[8] ] = []
                selected_entries += F[0:10]
                consider = True
                for idx in range(10, len(F)):
                    curr_val = float(F[idx])
                    selected_entries.append( curr_val )
                    if curr_val == -1.0:
                        consider = False
                        break
                if consider:
                    exons[ F[8] ].append(selected_entries)
                else:
                    exons[ F[8] ].append( None )

        with open(output.rel_use_filtered, "w") as out_table:
            for ex in exons:
                if None in exons[ex]:
                    continue
                for site_list in exons[ex]:
                    out_table.write("%s\n" % "\t".join([str(i) for i in site_list]))

        # same for expression table
        exons = {}
        with open(input.expressions, "r") as expr_in:
            for line in expr_in:
                selected_entries = []
                F = line.rstrip().split("\t")
                if F[8] not in exons:
                    exons[ F[8] ] = []
                selected_entries += F[0:10]
                consider = True
                for idx in range(10, len(F)):
                    curr_val = float(F[idx])
                    selected_entries.append( curr_val )
                    if curr_val == -1.0:
                        consider = False
                        break
                if consider:
                    exons[ F[8] ].append(selected_entries)
                else:
                    exons[ F[8] ].append( None )

        with open(output.expressions_filtered, "w") as out_table_expr:
            for ex in exons:
                if None in exons[ex]:
                    continue
                for site_list in exons[ex]:
                    out_table_expr.write("%s\n" % "\t".join([str(i) for i in site_list]))

#-------------------------------------------------------------------------------
# get the relative position of poly(A) sites within the exon
#-------------------------------------------------------------------------------
rule rel_pos_of_pAs:
    ##LOCAL##
    input:
        "{study}/relative_usages.filtered.tsv"
    output:
        "{study}/relative_usages.relPos_per_pA.out"
    params:
        script_dir = config['dir.scripts'],
        py2_env_path = config['py2_env_path']
    shell:
        '''
        {params.py2_env_path}/python {params.script_dir}/relative-pas-position-within-exon.py \
        --relUsage {input} \
        > {output}
        '''

#-------------------------------------------------------------------------------
# get weighted average exon length
#-------------------------------------------------------------------------------
rule weighted_average_exon_length:
    ##LOCAL""
    input:
        relUse = "{study}/relative_usages.filtered.tsv",
        relPos = "{study}/relative_usages.relPos_per_pA.out",
        header = "{study}/relative_usages.header.out"
    output:
        "{study}/weighted_avg_exonLength.tsv"
    params:
        script_dir = config['dir.scripts'],
        py2_env_path = config['py2_env_path']
    run:
        sample_names = " ".join([curr_id.rstrip() for curr_id in open(input.header, "r").readlines()])
        shell('''
        {params.py2_env_path}/python {params.script_dir}/calculate-average-3pUTR-length.py \
        --samples {sample_names} \
        --relativePos={input.relPos} \
        --relUsage {input.relUse} \
        > {output}
        ''')

#-------------------------------------------------------------------------------
# plot the cumulative distribution for the weighted average exon lengths
# per sample
#-------------------------------------------------------------------------------
rule plot_average_exon_length:
    ##LOCAL##
    input:
        "{study}/weighted_avg_exonLength.filtered.tsv"
    output:
        "{study}/weighted_avg_exonLength.CDFs.pdf"
    params:
        script_dir = config['dir.scripts'],
        py2_env_path = config['py2_env_path'],
        study = "{study}"
    shell:
        '''
        Rscript {params.script_dir}/rs-plot-ecdfs.R \
        --main={params.study}_samples \
        --pdf={output} \
        --file={input}
        '''

configfile: "config.yaml"

from snakemake.utils import makedirs
from snakemake.utils import listfiles

import os
import numpy as np
import yaml

localrules: create_log_dir, create_log_dir_TIN_calc, bam_index, merge_TIN_vals, TIN_assessment, plot_full_trans_TIN_distributions, infer_valid_samples, all

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
        expand("{study}/valid_samples/inferred_valid_samples.out", study = config['studies'])



################################################################################
# chronological order of the processing rules
################################################################################

rule create_log_dir:
    ## LOCAL ##
    ''' This step creates the log directory, if necessary.
    This is required when jobs are submitted and the
    job output should be written to these files.
    '''
    output:
        "{study}/logs/created_log.tmp"
    run:
        makedirs( wildcards.study + "/logs/cluster_logs" )
        shell("touch {output}")

rule create_log_dir_TIN_calc:
    ## LOCAL ##
    ''' This step creates the log directory, if necessary.
    This is required when jobs are submitted and the
    job output should be written to these files.
    '''
    output:
        config["dir.samplewise_TIN"] + "/logs/created_log.tmp"
    run:
        makedirs( config["dir.samplewise_TIN"] + "/logs/cluster_logs" )
        shell("touch {output}")


#-------------------------------------------------------------------------------
# create a bam index
#-------------------------------------------------------------------------------

rule bam_index:
    ##LOCAL##
    input:
        bam = config["dir.input"] + "/{sample}.bam"
    output:
        bam_index = config["dir.input"] + "/{sample}.bam.bai"
    shell:
        '''
        samtools index {input.bam} {output.bam_index}
        '''

#-------------------------------------------------------------------------------
# assess coverage bias transcript wide (via TIN measure)
#-------------------------------------------------------------------------------
rule assess_TIN_bias_transcript_wide:
    input:
        created_log = config["dir.samplewise_TIN"] + "/logs/created_log.tmp",
        bam = lambda wildcards: config["dir.input"] + "/" + config[wildcards.sample]["bam"] + ".bam",
        bam_index = lambda wildcards: config["dir.input"] + "/" + config[wildcards.sample]["bam"] + ".bam.bai"
    output:
        config["dir.samplewise_TIN"] + "/{sample}.tsv"
    params:
        script_dir = config['dir.scripts'],
        py2_env_path = config['py2_env_path'],
        transcripts = config['transcripts'],
        min_raw_reads = config['bias.rawReads.min'],
        sample_size = config['bias.transcripts.sampleSize'],
        cluster_log = config["dir.samplewise_TIN"] + "/logs/cluster_logs/{sample}.log"
    log:
        config["dir.samplewise_TIN"] + "/logs/{sample}.log"
    shell:
        '''
        {params.py2_env_path}/py2_paqr/bin/python {params.script_dir}/calculate-TIN-values.py \
        -i {input.bam} \
        -r {params.transcripts} \
        -c {params.min_raw_reads} \
        --names {wildcards.sample} \
        -n {params.sample_size} \
        > {output} \
        2> {log}
        '''

#-------------------------------------------------------------------------------
# merge single sample TIN-vals tables
#-------------------------------------------------------------------------------
rule merge_TIN_vals:
    ##LOCAL##
    input:
        create_log_dir = "{study}/logs/created_log.tmp",
        tables = lambda wildcards: expand( config["dir.samplewise_TIN"] + "/{na}.tsv", na = config[wildcards.study]['samples'])
    output:
        "{study}/bias.transcript_wide.TIN.tsv"
    params:
        script_dir = config['dir.scripts'],
        py2_env_path = config['py2_env_path']
    log:
        "{study}/logs/merge_TIN_vals.log"
    shell:
        '''
        {params.py2_env_path}/py2_paqr/bin/python {params.script_dir}/merge-TIN-tables.py \
        --verbose \
        --input {input.tables} \
        > {output} \
        2> {log}
        '''

#-------------------------------------------------------------------------------
# create a table with the medianTIN values per sample
# and
# create a table with samples that are above the medTIN cutoff
#-------------------------------------------------------------------------------
rule TIN_assessment:
    ## LOCAL ##
    input:
        table = "{study}/bias.transcript_wide.TIN.tsv"
    output:
        overview = "{study}/bias.TIN.median_per_sample.tsv",
        no_bias = "{study}/no_bias_samples.out"
    params:
        script_dir = config['dir.scripts'],
        bias_cutoff = config['bias.median.cutoff'],
        dummy_file = "{study}/DUMMY_USED.out"
    run:
        transcript_vals = {}
        names = []
        medians = []
        with open(input.table, "r") as f:
            for line in f:
                F = line.rstrip().split("\t")
                if line.startswith("transcript"):
                    # header line; store names
                    for idx in range(1, len(F)):
                        transcript_vals[ idx - 1 ] = []
                        names.append( F[idx] )
                else:
                    for idx in range(1, len(F)):
                        val = float( F[idx])
                        if val != 0.0:
                            transcript_vals[ idx -1 ].append( val )

        for idx in sorted(transcript_vals):
            medians.append( np.median( transcript_vals[idx]) )

        no_bias_samples = 0
        with open(output.no_bias, "w") as noBias:
            with open(output.overview, "w") as out:
                out.write("sample\tmedian_TIN\n")
                for idx in range(len(names)):
                    curr_name = names[idx]
                    curr_median = medians[idx]
                    out.write("%s\t%f\n" % (curr_name, curr_median))

                    if "control" in config[curr_name]:
                        ctl_idx = names.index( config[curr_name]['control'] )
                        ctl_median = medians[ctl_idx]
                        if curr_median > params.bias_cutoff and ctl_median > params.bias_cutoff:
                            no_bias_samples += 2
                            noBias.write("%s\n" % config[curr_name]['control'])
                            noBias.write("%s\n" % curr_name)

            if no_bias_samples == 0:
            # create an artificial pair of samples
                with open(params.dummy_file, "w") as dummy:
                    dummy.write("No valid samples\nOnly dummy samples are used\n")
                for idx in range(len(names)):
                    curr_name = names[idx]
                    if "control" in config[curr_name]:
                        noBias.write("%s\n" % config[curr_name]['control'])
                        noBias.write("%s\n" % curr_name)
                        break

#-------------------------------------------------------------------------------
# plot the boxplots of full-transcript TIN vals per sample
#-------------------------------------------------------------------------------
rule plot_full_trans_TIN_distributions:
    ##LOCAL##
    input:
        table = "{study}/bias.transcript_wide.TIN.tsv"
    params:
        script_dir = config['dir.scripts']
    output:
        "{study}/bias.transcript_wide.TIN.boxplots.pdf"
    shell:
        '''
        Rscript {params.script_dir}/boxplots-TIN-distributions.R \
	--file {input.table} \
	--pdf {output}
        '''

#-------------------------------------------------------------------------------
# create symbolic links for all valid samples
#-------------------------------------------------------------------------------
rule infer_valid_samples:
    ##LOCAL##
    input:
        no_bias = "{study}/no_bias_samples.out",
        prev_plot = "{study}/bias.transcript_wide.TIN.boxplots.pdf"
    output:
        touch("{study}/valid_samples/inferred_valid_samples.out")
    params:
        dummy_file = "{study}/DUMMY_USED.out"
    run:
        import os.path
        # iterate over samples in input and create symbolic links
        with open(input.no_bias, 'r') as ifile:
            line_cnt = 0
            for line in ifile:
                line_cnt += 1
                F = line.rstrip().split("\t")
                curr_uuid = F[0]
                if os.path.isfile(params.dummy_file):
                    # bam_file = "../../" + config["dir.input"] + "/dummy_" + str(line_cnt) + ".bam"
                    pass
                else:
                    bam_file = "../../" + config["dir.input"] + "/" + config[curr_uuid]['bam'] + ".bam"
                bam_index_file = bam_file + ".bai"
                res_dir = wildcards.study
                # create a soft link for the bam and for the bai
                # shell("cd {res_dir}/valid_samples/;ln -fs {bam_file} {curr_uuid}.bam")
                # shell("cd {res_dir}/valid_samples/;ln -fs {bam_index_file} {curr_uuid}.bam.bai")
                if not os.path.isfile(params.dummy_file):
                    shell("cd {res_dir}/valid_samples/;ln -fs {bam_file} {curr_uuid}.bam")
                    shell("cd {res_dir}/valid_samples/;ln -fs {bam_index_file} {curr_uuid}.bam.bai")

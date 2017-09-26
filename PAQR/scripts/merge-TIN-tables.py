# -*- coding: utf-8 -*-
'''
Created on Thu Mar 09 09:30:36 2017

Use tables with TIN values per transcript for single samples and merge them to one table.

@author: schmiral
'''
__date__ = "Thu Mar 09 09:30:36 2017"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# import
import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter
import pymongo
import datetime
import time


# parse input arguments

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
                    
parser.add_argument("--input",
                    dest="input",
                    nargs="*",
                    help="file names of sample-specific tables")

syserr = sys.stderr.write
sysout = sys.stdout.write

def main(options):

    nr_input = len(options.input)
    header = "transcript"
    vals_per_transcript = {}
    for f in options.input:
        with open(f, "r") as ifile:
            for line in ifile:
                F = line.rstrip().split("\t")
                if F[0] == "transcript":
                    # append header
                    header += "\t%s" % F[1]
                else:
                    if F[0] not in vals_per_transcript:
                        # create entry for current transcript
                        vals_per_transcript[F[0]] = []
                    # append current value to list
                    vals_per_transcript[F[0]].append( F[1] )

    # output merged table
    sysout("%s\n" % header)
    for trans in vals_per_transcript:
        if len( vals_per_transcript[trans] ) != nr_input:
            syserr("[INFO] Could not find the expected number of entries (%i) for transcript %s; did not include this transript into final table\n" % (nr_input, trans))
            continue
        sysout("%s\t%s\n" % (trans, "\t".join(vals_per_transcript[trans])))


if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception, e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" %
                   start_date)

        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)

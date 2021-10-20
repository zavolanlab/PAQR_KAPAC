#!/usr/bin/env python
"""
My template:
"""

__date__ = "2016-07-07"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("--relativePos",
                    dest="relativePos",
                    help="path and name of the file with the positions of the poly(A) sites relative to the exon")

parser.add_argument("--relUsage",
                    dest="rel_usages",
                    help="path and name of the file with the relative usages")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""

    posPerPas = {}
    terminalExons = {}

    # parse the relative poly(A) site positions
    with open( options.relativePos, "r") as relPosIn:
        for line in relPosIn:
            F = line.rstrip().split("\t")
            posPerPas[ F[0] ] = F[1]

    n = 0
    with open( options.rel_usages, "r") as relUseIn:
        for line in relUseIn:
            n += 1
            F = line.rstrip().split("\t")
            if n == 1:
                # for first row, only get samples from header
                samples = F[10:]
                continue

            pas = F[3]
            relPos = float(posPerPas[ pas ])
            nr = 0
            for sample_usage in F[10:len(F)]:
                sample_usage = float(sample_usage)
                if F[8] not in terminalExons:
                    terminalExons[ F[8] ] = []
                    for k in range(10,len(F)):
                        terminalExons[ F[8] ].append(0)
                if terminalExons[ F[8] ][nr] != -1:
                    if sample_usage == -1:
                        terminalExons[ F[8] ][nr] = -1
                    else:
                        terminalExons[ F[8] ][nr] += sample_usage * relPos
                nr += 1

    # output
    sysout("exon\t%s\n" % "\t".join(samples))
    for termExon in terminalExons:
        sysout("%s" % termExon)
        for relative_length in terminalExons[ termExon ]:
            if relative_length == -1:
                sysout("\t-1")
            else:
                sysout("\t%.2f" % (relative_length / 100.0) )
        sysout("\n")

if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception:
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

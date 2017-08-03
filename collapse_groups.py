#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import gzip

__author__ = "Alexandra Scott (ajscott@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2017-07-19 11:41 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
    	gs_cluster_linkage.py\n\
    	author: " + __author__ + "\n\
    	version: " + __version__ + "\n\
    	description: collapse groups so that each variant only present in one group")
    parser.add_argument('-i', '--input',
                        metavar='FILE', dest='input_path',
                        type=str, default=None,
                        help='Input file containing comma separated cluster variant IDs')

    # parse the arguments
    args = parser.parse_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)

    # send back the user input
    return args

# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data

def make_chains(group, chains):
    to_merge=group
    merged = False
    remove = set()
    if len(chains) != 0:
        for x in range(0,len(chains)):
            for index in group:
                if index in chains[x]:
                    to_merge = to_merge | chains[x]
                    merged = True
                    remove.add(x)
    chains.append(to_merge)
    if merged==True:
        for y in remove:
            chains.pop(y)
    return chains    



# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        input_file = sys.stdin
    else:
        input_file = get_file(args.input_path)

    # run functions
    chains = []
    for line in input_file:
        chains = make_chains(set(line.rstrip('\n').split(',')), chains)

    # print results
    for x in chains:
        print(",".join(x))

    # close files
    input_file.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:
        # ignore SIGPIPE
            raise
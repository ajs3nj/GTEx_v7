#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import gzip

__author__ = "Alexandra Scott (ajscott@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-12-06 12:55 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
    	gs_cluster_linkage.py\n\
    	author: " + __author__ + "\n\
    	version: " + __version__ + "\n\
    	description: check that overlapping gs varaints are strongly linked and should all be included in cluster (only gives one cluster per set of variants)")
    parser.add_argument('-i', '--input',
                        metavar='FILE', dest='input_path',
                        type=str, default=None,
                        help='Input file containing comma separated cluster variant IDs')
    parser.add_argument('-l', '--linkage',
                        metavar='FILE', dest='linkage_path',
                        required=True,
                        type=str, default=None,
                        help='File containing variant linkage values')
    parser.add_argument('-r2', '--r2-threshold',
                        dest='r2_threshold',
                        required=False,
                        type=float,default=0.97)

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

# get linkage information
def get_linkage(linkage_file):
    linkage = {}
    for l in linkage_file:
        l_split = l.rstrip('\n').split('\t')
        linkage[l_split[0]] = {}
        linkage[l_split[0]][l_split[1]]=l_split[2]
    return linkage

# check if all variants in cluster are linked, if not choose one to remove
def check_cluster(cluster,linkage,r2):
    unlinked_counts = {}
    for test_variant in cluster:
        unlinked_counts[test_variant] = 0
        for compare_variant in linkage[test_variant]:
            if float(linkage[test_variant][compare_variant]) < r2:
                unlinked_counts[test_variant] += 1
    if all(x==0 for x in unlinked_counts.values()):
        return "NA"
    else:
        return max(unlinked_counts, key=unlinked_counts.get)


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

    # get input data
    linkage_file = get_file(args.linkage_path)
    r2 = args.r2_threshold

    # run functions
    linkage = get_linkage(linkage_file)

    for line in input_file:
        cluster = line.rstrip('\n').split(',')
        delete = check_cluster(cluster,linkage,r2)
        while delete != "NA":
            cluster.remove(delete)
            # check if we need to find additional clusters for deleted variant
            if any(x>r2 for x in linkage[delete]):
                clust2 = set()
                clust2.add(delete)
                for var in linkage[delete]:
                    if linkage[delete][var] < r2:
                        clust2.add(var)
                del2=check_cluster(clust2,linkage,r2)
                while del2 != 'NA':
                    clust2.remove(del2)
                    del2=check_cluster(clust2,linkage,r2)
                print(",".join(clust2))
            delete = check_cluster(cluster,linkage,r2)
        print(",".join(cluster))

    # close files
    input_file.close()
    linkage_file.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:
        # ignore SIGPIPE
            raise



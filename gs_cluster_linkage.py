#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

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
                        type=float,default=0.9)

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
        if l_split[0] not in linkage:
            linkage[l_split[0]] = {}
        linkage[l_split[0]][l_split[1]]=l_split[2]
    return linkage

def subset_linked(cluster,linkage,r2):
    linked = set()
    for var1 in cluster:
        for var2 in cluster:
            if var1 != var2:
                if linkage[var1][var2] >= r2:
                    linked.add(var1)
                    linked.add(var2)                
    return linked


# check if all variants in cluster are linked, if not choose one to remove
def check_cluster(cluster,linkage,r2):
    unlinked_counts = {}
    for test_variant in cluster:
        unlinked_counts[test_variant] = 0
        for compare_variant in cluster:
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
        # get input group
        group = line.rstrip('\n').split(',')

        # only create clusters with variants that are tightly linked to another variant
        linked_cluster = subset_linked(group,linkage,r2)

        # make clusters
        while len(linked_cluster) > 1:
            cluster = linked_cluster.copy()
            delete = check_cluster(cluster,linkage,r2)
            while delete != "NA":
                cluster.remove(delete)
                delete = check_cluster(cluster,linkage,r2)
            if len(cluster) > 1:
                print(",".join(cluster))
                linked_cluster = subset_linked(linked_cluster-cluster,linkage,r2)
            else:
                break

   

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



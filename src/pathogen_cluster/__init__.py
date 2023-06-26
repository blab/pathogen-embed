"""
pathogen-cluster.

"""

__version__ = "1.0.0"
__author__ = 'Sravani Nanduri, John Huddleston'
__credits__ = 'Bedford Lab, Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA'

import argparse
import sys
from .pathogen_cluster import pathogen_cluster

def make_parser():
    parser = argparse.ArgumentParser(description = "Reduced dimension embeddings for pathogen sequences", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # correct
    parser.add_argument("--embedding", help="The embedding to assign clustering labels to via HDBSCAN (https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html)")
    parser.add_argument("--random-seed", default = 314159, type=int, help="an integer used for reproducible results.")
    parser.add_argument("--min-size", type=int, default=5, help="minimum cluster size for HDBSCAN")
    parser.add_argument("--min-samples", type=int, default=5, help="minimum number of sample to seed a cluster for HDBSCAN. Lowering this value reduces number of samples that do not get clustered.")
    parser.add_argument("--distance-threshold", type=float, help="The float value for the distance threshold by which to cluster data in the embedding and assign labels via HDBSCAN. If no value is given in cluster-data or cluster-threshold, the default distance threshold of 0.0 will be used.")
    parser.add_argument("--output-dataframe", help="a csv file outputting the embedding with the strain name and its components.")
    parser.add_argument("--output-figure", help="outputs a PDF with a plot of the embedding colored by cluster")

    return parser


def run(argv):
    args = make_parser().parse_args(argv)
    try:
        return pathogen_cluster(args)
    except Exception as error:
        print(error, file=sys.stderr)
        sys.exit(1)

"""
pathogen-distance.

"""

__version__ = "1.0.0"
__author__ = 'Sravani Nanduri, John Huddleston'
__credits__ = 'Bedford Lab, Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA'

import argparse
import sys
from .pathogen_distance import pathogen_distance

def make_parser():
    parser = argparse.ArgumentParser(description = "Reduced dimension embeddings for pathogen sequences", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--alignment", help="an aligned FASTA file to create a distance matrix with. Make sure the strain order in this file matches the order in the distance matrix.")
    parser.add_argument("--indel-distance", action="store_true", help="include insertions/deletions in genetic distance calculations")
    parser.add_argument("--output", help="a csv file outputting the distance matrix annotated with strain names as the columns")
   
    return parser


def run(argv):
    args = make_parser().parse_args(argv)
    try:
        return pathogen_distance(args)
    except Exception as error:
        print(error, file=sys.stderr)
        sys.exit(1)

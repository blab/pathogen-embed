
"""
Stub function and module used as a setuptools entry point.
"""

import pathogen_cluster
from sys import argv, exit

# Entry point for setuptools-installed script and bin/embed dev wrapper.
def main():
    return pathogen_cluster.run( argv[1:] )

# Run when called as `python -m embed`, here for good measure.
if __name__ == "__main__":
    exit( main() )
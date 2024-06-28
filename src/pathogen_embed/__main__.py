import argparse
import sys
from sys import argv
from .pathogen_embed import embed, distance, cluster

def autoOrFloat(values):
    if values == "auto":
        return "auto"
    else:
        try:
            val = float(values)
            return val
        except ValueError:
            raise argparse.ArgumentTypeError(f"Invalid value: {values}. Must be a float or 'auto'.")


def make_parser_embed():
    parser = argparse.ArgumentParser(description = "Reduced dimension embeddings for pathogen sequences", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--alignment", nargs='+', help="one or more aligned FASTA files used for PCA, for PCA initialization of t-SNE, or to create distance matrices on the fly when distances are needed and none are provided.")
    parser.add_argument("--distance-matrix", nargs='+', help="one or more distance matrix CSVs with sequence names in a header row and the first column.")
    parser.add_argument("--separator", default=",", help="separator between columns in the given distance matrix")
    parser.add_argument("--indel-distance", action="store_true", help="include insertions/deletions in genetic distance calculations")
    parser.add_argument("--random-seed", default = 314159, type=int, help="an integer used for reproducible results.")
    parser.add_argument("--output-dataframe", help="a csv file outputting the embedding with the strain name and its components.")
    parser.add_argument("--output-figure", help="outputs a plot of the embedding")
    parser.add_argument("--embedding-parameters", help="The file containing the parameters by which to tune the embedding. The values from the first record of this file will override default values or values provided by the command line arguments.")
    parser.add_argument("--output-pairwise-distance-figure", help="a scatterplot correlating the genetic vs Euclidean distances")

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    pca = subparsers.add_parser("pca", description="Principal Component Analysis", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pca.add_argument(
        "--encoding",
        default="integer",
        choices=["integer", "genotype", "simplex", "biallelic"],
        help="""method to use to encode the given sequence alignment as a matrix for input to PCA.
        The "integer" encoding maps each ACGT nucleotide character to an integer (A to 1, G to 2, C to 3, T to 4) and all other characters to 5.
        The "genotype" encoding maps each ACGT nucleotide to a one-hot binary encoding (A to [1, 0, 0, 0], C to [0, 1, 0, 0], etc.) and all other characters to [0, 0, 0, 0].
        The "simplex" encoding maps each ACGT nucleotide to a simplex space as described in Stormo 2011 (A to [1, -1, -1], C to [-1, 1, -1], etc.) and all other characters to [0, 0, 0].
        The "biallelic" encoding identifies all biallelic positions in the input alignment, encodes each allele that matches the first ("reference") sequence as a 0 and each alternate allele as a 1.
        """
    )
    pca.add_argument("--components", default=10, type=int, help="the number of components for PCA")
    pca.add_argument("--explained-variance", help="the path for the CSV explained variance for each component")

    tsne = subparsers.add_parser("t-sne", description="t-distributed Stochastic Neighborhood Embedding", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    tsne.add_argument("--pca-encoding", default="integer", choices=["integer", "genotype", "simplex", "biallelic"], help="method to use to encode the given sequence alignment as a matrix for input the PCA embedding that initializes the t-SNE embedding. See help for the PCA embedding subcommand for more details.")
    tsne.add_argument("--components", default=2, type=int, help="the number of components for t-SNE")
    tsne.add_argument("--perplexity", default=30.0, type=float, help="The perplexity is related to the number of nearest neighbors. Because of this, the size of the dataset is proportional to the best perplexity value (large dataset -> large perplexity). Values between 5 and 50 work best. The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")
    tsne.add_argument("--learning-rate", default="auto", type=autoOrFloat, help="The learning rate for t-SNE is usually between 10.0 and 1000.0. Values out of these bounds may create innacurate results. The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")

    umap = subparsers.add_parser("umap", description="Uniform Manifold Approximation and Projection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    umap.add_argument("--components", default=2, type=int, help="the number of components for UMAP")
    umap.add_argument("--nearest-neighbors", default=200, type=int, help="Nearest neighbors controls how UMAP balances local versus global structure in the data (finer detail patterns versus global structure). This value is proportional to the size of the data (large dataset -> large nearest neighbors). The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")
    umap.add_argument("--min-dist", default=.5, type=float, help="Minimum Distance controls how tightly packed the UMAP embedding is. While it does not change the structure of the data, it does change the embedding's shape. The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")

    mds = subparsers.add_parser("mds", description="Multidimensional Scaling", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    mds.add_argument("--components", default=10, type=int, help="the number of components for MDS")
    mds.add_argument("--stress", help="the path for the CSV stress for the embedding")

    return parser

def make_parser_distance():
    parser = argparse.ArgumentParser(description = "Hamming distance (optionally indel sensitive) similarity matrix for pathogen sequences", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--alignment", required = True,  help="an aligned FASTA file to create a distance matrix with. Make sure the strain order in this file matches the order in the distance matrix.")
    parser.add_argument("--indel-distance", action="store_true", help="include insertions/deletions in genetic distance calculations")
    parser.add_argument("--output", required = True, help="a csv file outputting the distance matrix annotated with strain names as the columns")

    return parser

def make_parser_cluster():
    parser = argparse.ArgumentParser(description = "HDBSCAN clustering for reduced dimension embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--embedding", required = True, help="The embedding to assign clustering labels to via HDBSCAN (https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html)")
    parser.add_argument("--label-attribute", help="the name of the cluster used to label the column in the resulting dataframe")
    parser.add_argument("--random-seed", default = 314159, type=int, help="an integer used for reproducible results.")
    parser.add_argument("--min-size", type=int, default=5, help="minimum cluster size for HDBSCAN")
    parser.add_argument("--min-samples", type=int, default=5, help="minimum number of sample to seed a cluster for HDBSCAN. Lowering this value reduces number of samples that do not get clustered.")
    parser.add_argument("--distance-threshold", type=float, help="The float value for the distance threshold by which to cluster data in the embedding and assign labels via HDBSCAN. If no value is given in distance-threshold, the default distance threshold of 0.0 will be used.")
    parser.add_argument("--output-dataframe", required = True, help="a csv file outputting the embedding with the strain name and its components.")
    parser.add_argument("--output-figure", help="outputs a PDF with a plot of the embedding colored by cluster")

    return parser

def run_embed():
    args = make_parser_embed().parse_args(argv[1:])
    return embed(args)

def run_distance():
    args = make_parser_distance().parse_args(argv[1:])
    return distance(args)

def run_cluster():
    args = make_parser_cluster().parse_args(argv[1:])
    return cluster(args)

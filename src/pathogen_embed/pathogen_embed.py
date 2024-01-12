# Ignore warnings from Numba deprecation:
# https://numba.readthedocs.io/en/stable/reference/deprecation.html#deprecation-of-object-mode-fall-back-behaviour-when-using-jit
# Numba is required by UMAP.
from numba.core.errors import NumbaDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)

import argparse
import Bio.SeqIO
from collections import OrderedDict
import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from scipy.spatial.distance import squareform, pdist
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS
import sys
from umap import UMAP


def find_ranges(positions):
    """
    Find ranges of adjacent integers in the given list of integers and
    return a dictionary of gap lengths indexed by start position.

    >>> find_ranges([])
    {}
    >>> find_ranges([0])
    {0: 1}
    >>> find_ranges([0, 2, 3, 4])
    {0: 1, 2: 3}
    >>> find_ranges([2, 3, 4])
    {2: 3}
    >>> find_ranges([2, 3, 4, 6, 7, 9])
    {2: 3, 6: 2, 9: 1}

    """
    ranges = {}
    start = 0
    end = 0
    for i in range(len(positions)):
        # If the next position is one greater than the current position, update
        # the end point to the next position.
        if i < len(positions) - 1 and positions[i] + 1 == positions[i + 1]:
            end = i + 1
        # Otherwise, if the next position is more than one away or we're at the
        # end of the list, save the current range and set the next range to
        # start and end at the next position.
        else:
            # If the range is a singleton, output only that value. Otherwise,
            # output the range from start to end.
            ranges[positions[start]] = end - start + 1

            start = i + 1
            end = i + 1

    return ranges


def get_hamming_distances(genomes, count_indels=False):
    """Calculate pairwise Hamming distances between the given list of genomes
    and return the nonredundant array of values for use with scipy's squareform function.
    Bases other than standard nucleotides (A, T, C, G) are ignored. Treat indels as a single event.

    Parameters
    ----------
    genomes : list
        a list of strings corresponding to genomes that should be compared
    count_indels : boolean
        true means indels are counted in the distance calculation, false if not.
        the default value is false.

    Returns
    -------
    list
        a list of distinct Hamming distances as a vector-form distance vector

    >>> genomes = ["ATGCT", "ATGCT", "ACGCT"]
    >>> get_hamming_distances(genomes, True)
    [0, 1, 1]
    >>> get_hamming_distances(["AT--CT", "AC--CT"], True)
    [1]
    >>> genomes = ["AT-GCT", "AT--CT", "AC--CT"]
    >>> get_hamming_distances(genomes, True)
    [1, 2, 1]
    >>> genomes = ["ACTGG", "A--GN", "A-NGG"]
    >>> get_hamming_distances(genomes, True)
    [1, 1, 1]

    When counting indels, we ignore leading and trailing gaps that indicate
    different sequence lengths and not biological events.

    >>> genomes = ["ACTGTA", "A--CCA", "A--GT-"]
    >>> get_hamming_distances(genomes, True)
    [3, 1, 2]
    >>> genomes = ["ACTGTA", "A--CCA", "---GT-"]
    >>> get_hamming_distances(genomes, True)
    [3, 0, 2]

    When not counting indels, we ignore gaps altogether.

    >>> genomes = ["ATGCT", "ATGCT", "ACGCT"]
    >>> get_hamming_distances(genomes)
    [0, 1, 1]
    >>> genomes = ["AT-GCT", "AT--CT", "AC--CT"]
    >>> get_hamming_distances(genomes)
    [0, 1, 1]
    >>> genomes = ["ACTGG", "A--GN", "A-NGG"]
    >>> get_hamming_distances(genomes)
    [0, 0, 0]
    >>> genomes = ["ACTGTA", "A--CCA", "A--GT-"]
    >>> get_hamming_distances(genomes)
    [2, 0, 2]

    """

    # Define an array of valid nucleotides to use in pairwise distance calculations.
    # Using a numpy array of byte strings allows us to apply numpy.isin later.
    nucleotides = np.array([b'A', b'T', b'C', b'G'])

    # Convert genome strings into numpy arrays to enable vectorized comparisons.
    genome_arrays = [
        np.frombuffer(genome.encode(), dtype="S1")
        for genome in genomes
    ]

    # Precalculate positions of valid bases (A, T, C, and G) in each genome to speed up later comparisons.
    valid_bases = [
        np.isin(genome_array, nucleotides)
        for genome_array in genome_arrays
    ]

    # Calculate Hamming distance between all distinct pairs of genomes at valid bases.
    # The resulting list is a reduced representation of a symmetric matrix that can be
    # converted to a square matrix with scipy's squareform function:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
    total_genomes = len(genomes)
    alignment_length = len(genomes[0])
    hamming_distances = []
    for i in range(total_genomes):
        # Only compare the current genome, i, with all later genomes.
        # This avoids repeating comparisons or comparing each genome to itself.
        if count_indels:
            i_gaps = np.where(genome_arrays[i] == b"-")
            i_gap_ranges = find_ranges(i_gaps[0])

        for j in range(i + 1, total_genomes):
            # Find all mismatches at valid nucleotide bases.
            distance = ((genome_arrays[i] != genome_arrays[j]) & valid_bases[i] & valid_bases[j]).sum()

            if count_indels:
                j_gaps = np.where(genome_arrays[j] == b"-")
                j_gap_ranges = find_ranges(j_gaps[0])

                # Mismatched gaps include total number of different start
                # positions plus the number of the same start positions with
                # different lengths. Note that we ignore gaps that start at the
                # beginning of the alignment which occur because of differing
                # sequence lengths and not necessarily a biological event.
                num_indel = 0
                for gap_start, gap_length  in j_gap_ranges.items():
                    # Skip leading gaps.
                    if gap_start == 0:
                        continue

                    # Skip trailing gaps.
                    if gap_start + gap_length == alignment_length:
                        continue

                    if gap_start not in i_gap_ranges or i_gap_ranges[gap_start] != gap_length:
                        num_indel += 1

                distance += num_indel

            hamming_distances.append(distance)

    return hamming_distances


def embed(args):
    # Setting Random seed for numpy
    np.random.seed(seed=args.random_seed)

    if args.output_dataframe is None and args.output_figure is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)

    if args.alignment is None and args.command == "pca":
        print("You must specify an alignment for pca, not a distance matrix", file=sys.stderr)
        sys.exit(1)

    # getting or creating the distance matrix
    distance_matrix = None
    if args.distance_matrix is not None:
        distance_matrix  = pd.read_csv(args.distance_matrix, index_col=0)

    if args.alignment is not None:
        sequences_by_name = OrderedDict()

        for sequence in Bio.SeqIO.parse(args.alignment, "fasta"):
            sequences_by_name[sequence.id] = str(sequence.seq)

        sequence_names = list(sequences_by_name.keys())
        if args.command != "pca" and distance_matrix is None:
            # Calculate Distance Matrix
            hamming_distances = get_hamming_distances(
                list(sequences_by_name.values()),
                args.indel_distance,
            )
            distance_matrix = pd.DataFrame(squareform(hamming_distances))
            distance_matrix.index = sequence_names

    # Load embedding parameters from an external CSV file, if possible.
    external_embedding_parameters = None
    if args.embedding_parameters is not None:
        external_embedding_parameters_df = pd.read_csv(args.embedding_parameters)

        # Get a dictionary of additional parameters provided by the external
        # file to override defaults for the current method.
        external_embedding_parameters = external_embedding_parameters_df.to_dict("records")[0]

    if external_embedding_parameters is not None and "components" in external_embedding_parameters:
        n_components = int(external_embedding_parameters["components"])
        external_embedding_parameters["n_components"] = n_components
    else:
        n_components = args.components

    # Use PCA as its own embedding or as an initialization for t-SNE.
    if args.command == "pca" or args.command == "t-sne":
        sequence_names = list(sequences_by_name.keys())

        numbers = list(sequences_by_name.values())[:]
        for i in range(0,len(list(sequences_by_name.values()))):
            numbers[i] = re.sub(r'[^AGCT]', '5', numbers[i])
            numbers[i] = list(numbers[i].replace('A','1').replace('G','2').replace('C', '3').replace('T','4'))
            numbers[i] = [int(j) for j in numbers[i]]

        genomes_df = pd.DataFrame(numbers)
        genomes_df.columns = ["Site " + str(k) for k in range(0,len(numbers[i]))]


        #performing PCA on my pandas dataframe
        pca = PCA(n_components=n_components, svd_solver='full')
        principalComponents = pca.fit_transform(genomes_df)

        # Create a data frame from the PCA embedding.
        embedding = principalComponents
        embedding_df = pd.DataFrame(principalComponents)
        embedding_df.index = sequence_names

    if args.command == "t-sne":
        embedding_class = TSNE
        embedding_parameters = {
            "n_components": n_components,
            "metric": "precomputed",
            "init": principalComponents[:, :n_components],
            "perplexity": args.perplexity,
            "learning_rate": args.learning_rate,
            "random_state" : args.random_seed,
            "square_distances": True,
        }
    elif args.command == "umap":
        embedding_class = UMAP
        embedding_parameters = {
            "n_neighbors": args.nearest_neighbors,
            "min_dist": args.min_dist,
            "n_components": n_components,
            "init": "spectral",
            "random_state" : args.random_seed
        }
    elif args.command == "mds":
        embedding_class = MDS
        embedding_parameters = {
            "dissimilarity": "precomputed",
            "n_components": n_components,
            "n_jobs": 1,
            "n_init": 2
        }

    # Override defaults with parameter values passed through embedding parameters, if
    # possible.
    if external_embedding_parameters is not None and args.command != "pca":
        for key, value in external_embedding_parameters.items():
            if key in embedding_parameters:
                value_type = type(embedding_parameters[key])
                print(
                    f"INFO: Replacing embedding parameter {key} value of '{embedding_parameters[key]}' with '{value_type(value)}' provided by '{args.embedding_parameters}'.",
                    file=sys.stderr
                )
                embedding_parameters[key] = value_type(value)

    if args.command != "pca":
        embedder = embedding_class(**embedding_parameters)
        embedding = embedder.fit_transform(distance_matrix)

        # Output Embedding
        # create dictionary to be "wrapped" by write_json

        embedding_df = pd.DataFrame(embedding)
        embedding_df.index = list(distance_matrix.index)

    if args.command == "mds" or args.command == "pca":
        embedding_df.columns=[args.command + str(i) for i in range(1, n_components + 1)]
    else:
        embedding_df.columns = [args.command.replace('-', '') + "_x" , args.command.replace('-', '') + "_y"]

    if args.command == "pca":

        #add explained variance as the first row of the dataframe
        explained_variance = pd.DataFrame([round(pca.explained_variance_ratio_[i],4) for i in range(0,len(pca.explained_variance_ratio_))], columns=["explained variance"])
        explained_variance["principal components"] = [i for i in range(1, n_components + 1)]
        explained_variance.to_csv(args.explained_variance, index=False)

    if args.command == "mds":
        if args.stress:
            with open(args.stress, "w", encoding="utf-8") as oh:
                print(embedder.stress_, file=oh)

    if args.output_dataframe is not None:
        embedding_df.to_csv(args.output_dataframe, index_label="strain")

    if args.output_figure:
        plot_data = {
            "x": embedding[:, 0],
            "y": embedding[:, 1],
        }

        plot_df = pd.DataFrame(plot_data)
        ax = sns.scatterplot(
            data=plot_df,
            x="x",
            y="y",
            alpha=0.5,
        )
        plt.savefig(args.output_figure)
        plt.close()

def cluster(args):

    embedding_df = pd.read_csv(args.embedding, index_col=0)

    clustering_parameters = {
        **({"min_cluster_size": args.min_size} if args.min_size is not None else {}),
        **({"min_samples": args.min_samples} if args.min_samples is not None else {}),
        **({"cluster_selection_epsilon": args.distance_threshold} if args.distance_threshold is not None else {})
    }

    clusterer = hdbscan.HDBSCAN(**clustering_parameters)

    clusterer.fit(embedding_df)
    embedding_df[args.label_attribute] = clusterer.labels_.astype(str)

    if args.output_figure is not None:

        plot_data = {
            "x": embedding_df.to_numpy()[:, 0],
            "y": embedding_df.to_numpy()[:, 1],
        }

        plot_data["cluster"] = clusterer.labels_.astype(str)

        plot_df = pd.DataFrame(plot_data)
        ax = sns.scatterplot(
            data=plot_df,
            x="x",
            y="y",
            hue="cluster",
            alpha=0.5,
        )
        plt.savefig(args.output_figure)
        plt.close()

    if args.output_dataframe is not None:
        embedding_df.to_csv(args.output_dataframe, index_label="strain")

def distance(args):
    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(args.alignment, "fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = list(sequences_by_name.keys())

    hamming_distances = get_hamming_distances(
        list(sequences_by_name.values()),
        args.indel_distance,
    )
    distance_matrix = pd.DataFrame(squareform(hamming_distances))
    distance_matrix.index = sequence_names
    distance_matrix.columns = distance_matrix.index

    distance_matrix.to_csv(args.output)

# Ignore warnings from Numba deprecation:
# https://numba.readthedocs.io/en/stable/reference/deprecation.html#deprecation-of-object-mode-fall-back-behaviour-when-using-jit
# Numba is required by UMAP.
from numba.core.errors import NumbaDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)

import matplotlib; matplotlib.set_loglevel("critical")
import argparse
import Bio.AlignIO
import Bio.SeqIO
from collections import OrderedDict
import hdbscan
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
import re
from scipy.spatial.distance import squareform, pdist
from scipy.stats import linregress
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS
import sys
from umap import UMAP

def distance_tick_formatter(tick_val, tick_pos):
    return str(int(tick_val))

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


def encode_alignment_for_pca_by_integer(alignment_path):
    """Return a data frame representing a N x L matrix encoding of the alignment
    (with N alignment sequences of length L) in the given file path with each
    nucleotide character represented as an integer from 1 to 5 with all non-ACTG
    characters encoded as 5.

    The index of the data frame contains the sequence name for each input
    sequence.

    """
    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(alignment_path, "fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = sorted(list(sequences_by_name.keys()))
    numbers = []
    for sequence_name in sequence_names:
        sequence = sequences_by_name[sequence_name]
        sequence = re.sub(r'[^AGCT]', '5', sequence)
        sequence = list(sequence.replace('A','1').replace('G','2').replace('C', '3').replace('T','4'))
        sequence = [int(j) for j in sequence]
        numbers.append(sequence)

    genomes_df = pd.DataFrame(
        numbers,
        index=sequence_names,
    )

    return genomes_df


def encode_alignment_for_pca_by_genotype(alignment_path):
    """Return a data frame representing a N x 4L matrix encoding of the
    alignment (with N alignment sequences of length L) in the given file path
    with each position and allele combination (genotype) represented as a binary
    integer value. For example, the presence of an "A" allele in the first
    position of an alignment would cause the first four columns of the
    corresponding matrix row to be [1, 0, 0, 0].

    The index of the data frame contains the sequence name for each input
    sequence.

    """
    nucleotide_to_binary = {
        "A": np.array([1, 0, 0, 0], dtype=bool),
        "C": np.array([0, 1, 0, 0], dtype=bool),
        "G": np.array([0, 0, 1, 0], dtype=bool),
        "T": np.array([0, 0, 0, 1], dtype=bool),
    }
    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(alignment_path, "fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = sorted(list(sequences_by_name.keys()))
    numbers = []
    for sequence_name in sequence_names:
        sequence = sequences_by_name[sequence_name]
        sequence_numbers = []
        for position in range(len(sequence)):
            sequence_numbers.append(
                nucleotide_to_binary.get(
                    sequence[position].upper(),
                    np.array([0, 0, 0, 0], dtype=bool),
                )
            )

        numbers.append(np.array(sequence_numbers, dtype=bool).flatten())

    genomes_df = pd.DataFrame(
        numbers,
        index=sequence_names,
    ).astype(np.int8)

    return genomes_df


def encode_alignment_for_pca_by_simplex(alignment_path):
    """Return a data frame representing a N x 3L matrix encoding of the
    alignment (with N alignment sequences of length L) in the given file path
    with each position by three binary values corresponding to a simplex
    representation of the allele at that position as described in Stormo 2011
    (https://academic.oup.com/genetics/article/187/4/1219/6063441). For example,
    the presence of an "A" allele in the first position of an alignment would
    cause the first three columns of the corresponding matrix row to be [1, -1, -1].
    All non-ACGT characters get represented by a tuple of zeros.

    The index of the data frame contains the sequence name for each input
    sequence.

    """
    nucleotide_to_binary = {
        "A": np.array([1, -1, -1], dtype=np.int8),
        "C": np.array([-1, 1, -1], dtype=np.int8),
        "G": np.array([-1, -1, 1], dtype=np.int8),
        "T": np.array([1, 1, 1], dtype=np.int8),
    }
    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(alignment_path, "fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = sorted(list(sequences_by_name.keys()))
    numbers = []
    for sequence_name in sequence_names:
        sequence = sequences_by_name[sequence_name]
        sequence_numbers = []
        for position in range(len(sequence)):
            sequence_numbers.append(
                nucleotide_to_binary.get(
                    sequence[position].upper(),
                    np.array([0, 0, 0], dtype=np.int8),
                )
            )

        numbers.append(np.array(sequence_numbers, dtype=np.int8).flatten())

    genomes_df = pd.DataFrame(
        numbers,
        index=sequence_names,
    )

    return genomes_df


def encode_alignment_for_pca_by_biallelic(alignment_path):
    """Return a data frame representing a N x L matrix encoding of the alignment
    (with N alignment sequences of length L(b) where b identicates biallelic
    positions) in the given file path with each biallelic position represented
    by a binary value where 0 is the reference allele and 1 is the alternate
    allele. The function determines the "reference" allele based on the first
    sequence in the alignment. The alternate allele is the only other ACGT
    allele present in the remaining sequences at that position.

    The index of the data frame contains the sequence name for each input
    sequence.

    """
    valid_nucleotides = {"A", "C", "G", "T"}
    alignment = Bio.AlignIO.read(alignment_path, "fasta")
    reference_sequence = str(alignment[0].seq)

    # Find biallelic positions in the alignment.
    biallelic_positions = []
    for position in range(alignment.get_alignment_length()):
        distinct_nucleotides = set(alignment[:, position].upper()) & valid_nucleotides
        if len(distinct_nucleotides) == 2:
            biallelic_positions.append(position)

    # Create a binary matrix to represent alternate alleles at known biallelic
    # positions.
    numbers = np.zeros((len(alignment), len(biallelic_positions)), dtype=bool)
    sequence_names = []
    for sequence_index, sequence in enumerate(alignment):
        sequence_names.append(sequence.name)

        for position_index, position in enumerate(biallelic_positions):
            if (sequence[position].upper() in valid_nucleotides and
                sequence[position].upper() != reference_sequence[position].upper()):
                numbers[sequence_index, position_index] = True

    genomes_df = pd.DataFrame(
        numbers,
        index=sequence_names,
    ).astype(np.int8).sort_index()

    return genomes_df


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
        print("ERROR: PCA requires an alignment input to create the embedding.", file=sys.stderr)
        sys.exit(1)

    if args.alignment is None and args.command == "t-sne":
        print("ERROR: t-SNE requires an alignment input to initialize the embedding.", file=sys.stderr)
        sys.exit(1)

    if args.distance_matrix is not None and not all((distance_matrix.endswith(".csv") for distance_matrix in args.distance_matrix)):
        print("ERROR: The distance matrix input(s) must be in comma-separate value (CSV) format.", file=sys.stderr)
        sys.exit(1)

    if args.alignment is not None and args.distance_matrix is not None and len(args.alignment) != len(args.distance_matrix):
        print("ERROR: If giving multiple alignments and distance matrices the number of both must match.", file=sys.stderr)
        sys.exit(1)

    # Load distance matrices, if they have been provided, and sum them across
    # all inputs.
    distance_matrix = None
    if args.distance_matrix is not None and args.command != "pca":
        distance_path = args.distance_matrix[0]
        distance_df = pd.read_csv(distance_path, index_col=0).sort_index(axis=0).sort_index(axis=1)
        distance_array = distance_df.values.astype(float)

        # Add the numpy arrays element-wise
        for distance_path in args.distance_matrix[1:]:
            other_distance_df = pd.read_csv(distance_path, index_col=0).sort_index(axis=0).sort_index(axis=1)
            if not np.array_equal(distance_df.index.values, other_distance_df.index.values):
                print("ERROR: The given distance matrices do not have the same sequence names.", file=sys.stderr)
                sys.exit(1)

            other_distance_array = other_distance_df.values.astype(float)
            distance_array = distance_array + other_distance_array

        # Confirm that the distance matrix is square.
        if distance_array.shape[0] != distance_array.shape[1]:
            print("ERROR: Distance matrices must be square (with the same number of rows and columns).")
            sys.exit(1)

        distance_matrix = pd.DataFrame(distance_array, index=distance_df.index)

    # If we have alignments but no distance matrices and we need distances for
    # the method, calculate distances on the fly.
    if args.alignment is not None and distance_matrix is None and (args.command != "pca" or args.output_pairwise_distance_figure):
        for alignment in args.alignment:
            # Calculate a distance matrix on the fly and sort the matrix
            # alphabetically by sequence name, so we can safely sum values
            # across all inputs.
            new_distance_matrix = calculate_distances_from_alignment(
                alignment,
                args.indel_distance,
            ).sort_index(
                axis=0,
            ).sort_index(
                axis=1,
            )

            if distance_matrix is None:
                distance_matrix = new_distance_matrix
            else:
                # Confirm that we have the same strain names in the same order
                # for each matrix.
                if not np.array_equal(
                    new_distance_matrix.index.values,
                    distance_matrix.index.values,
                ):
                    print("ERROR: The given alignments do not have the same sequence names.", file=sys.stderr)
                    sys.exit(1)

                distance_matrix += new_distance_matrix

    # Load embedding parameters from an external CSV file, if possible.
    external_embedding_parameters = None
    if args.embedding_parameters is not None:
        if not args.embedding_parameters.endswith('.csv'):
            print("You must supply a CSV file for embedding parameters.", file=sys.stderr)
            sys.exit(1)
        else:
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
        if args.command == "pca":
            encoding = args.encoding
        else:
            encoding = args.pca_encoding

        genomes_dfs = []
        column_offset = 0
        for alignment in args.alignment:
            if encoding == "integer":
                genomes_df = encode_alignment_for_pca_by_integer(alignment)
            elif encoding == "genotype":
                genomes_df = encode_alignment_for_pca_by_genotype(alignment)
            elif encoding == "simplex":
                genomes_df = encode_alignment_for_pca_by_simplex(alignment)
            elif encoding == "biallelic":
                genomes_df = encode_alignment_for_pca_by_biallelic(alignment)
            else:
                raise NotImplementedError(f"PCA encoding by '{encoding}' is not supported.")

            genomes_df.columns = [
                "Site " + str(k)
                for k in range(column_offset, column_offset + genomes_df.shape[1])
            ]
            column_offset += genomes_df.shape[1]
            genomes_dfs.append(genomes_df)

        all_alignments_have_same_sequence_names = all(
            (
                np.array_equal(alignment.index.values, genomes_dfs[0].index.values)
                for alignment in genomes_dfs
            )
        )
        if not all_alignments_have_same_sequence_names:
            print("ERROR: The given alignments do not have the same sequence names.", file=sys.stderr)
            sys.exit(1)

        # Combine matrix encodings for all alignments, after confirming that
        # they represent the same sequence names.
        genomes_df = pd.concat(genomes_dfs, axis=1)
        sequence_names = genomes_df.index.values

        # If we're using PCA to initialize the t-SNE embedding, confirm that the
        # input alignments used for PCA have the same sequence names as the
        # input distance matrices.
        if (
                args.command == "t-sne" and
                not np.array_equal(distance_matrix.index.values, sequence_names)
        ):
            print("ERROR: The sequence names for the distance matrix inputs do not match the names in the alignment inputs.", file=sys.stderr)
            sys.exit(1)

        #performing PCA on my pandas dataframe
        pca = PCA(
            n_components=n_components,
            svd_solver='full',
            random_state=args.random_seed,
        )
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
        }
    elif args.command == "umap":
        embedding_class = UMAP
        embedding_parameters = {
            "n_neighbors": args.nearest_neighbors,
            "min_dist": args.min_dist,
            "n_components": n_components,
            "init": "spectral",
            "random_state" : args.random_seed,
            "n_jobs": 1,
        }
    elif args.command == "mds":
        embedding_class = MDS
        embedding_parameters = {
            "dissimilarity": "precomputed",
            "n_components": n_components,
            "n_jobs": 1,
            "n_init": 2,
            "random_state" : args.random_seed,
            "normalized_stress": False,
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
        #TODO: distance matrices are no longer symmetrics/not square? Check into this
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
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), dpi=300)
        ax.scatter(plot_df["x"], plot_df["y"], alpha=0.5)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        plt.tight_layout()
        plt.savefig(args.output_figure)
        plt.close()

    if args.output_pairwise_distance_figure:
        # Calculate pairwise Euclidean distances in the embedding.
        euclidean_distances = pdist(embedding_df.values)

        # Transform the genetic distance matrix to the condensed matrix format
        # such that each entry in the output corresponds to the entry at the
        # same position in the Euclidean distances.
        genetic_distances = squareform(distance_matrix)

        regression = linregress(genetic_distances, euclidean_distances)
        slope, intercept, r_value, p_value, std_err = regression
        intercept_sign = "+" if intercept >= 0 else "-"

        # Group Euclidean distances across the range of observed genetic
        # distances, so we can make a separate boxplot of Euclidean distances
        # per genetic distance value.
        genetic_distance_range = range(
            int(genetic_distances.min()),
            int(genetic_distances.max()) + 1,
        )
        grouped_euclidean_distances = []
        for genetic_distance in genetic_distance_range:
            grouped_euclidean_distances.append(
                euclidean_distances[genetic_distances == genetic_distance]
            )

        fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=300)
        boxplot = ax.boxplot(
            grouped_euclidean_distances,
            labels=list(genetic_distance_range),
            positions=list(genetic_distance_range),
            boxprops={
                "linewidth": 0.25,
            },
            medianprops={
                "color": "#999999",
            },
            whiskerprops={
                "linewidth": 0.25,
            },
            flierprops={
                "markersize": 1,
            },
            patch_artist=True,
        )

        # Set colors for boxes.
        for patch in boxplot["boxes"]:
            patch.set_facecolor("lightblue")

        # Plot linear fit behind the boxplots.
        ax.plot(
            genetic_distance_range,
            [slope * distance + intercept for distance in genetic_distance_range],
            color="#999999",
            zorder=-10,
        )

        ax.set_xlabel("Genetic distance")
        ax.set_ylabel("Euclidean distance")

        ax.xaxis.set_major_formatter(distance_tick_formatter)
        ax.xaxis.set_major_locator(MultipleLocator(5))

        ax.set_xlim(left=-1)
        ax.set_ylim(bottom=-1)

        ax.set_title(f"{args.command} (Pearson's $R^2={r_value:.2f}, y = {slope:.2f}x {intercept_sign} {np.abs(intercept):.2f}$)")
        plt.tight_layout()
        plt.savefig(args.output_pairwise_distance_figure)

def cluster(args):

    if not args.embedding.endswith('.csv'):
        print("You must supply a CSV file for the embedding.", file=sys.stderr)
        sys.exit(1)
    else:
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
        clusters = plot_df['cluster'].unique()

        fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=300)

        for i, cluster in enumerate(clusters):
            cluster_data = plot_df[plot_df['cluster'] == cluster]
            ax.scatter(cluster_data["x"], cluster_data["y"], label=f'Cluster {cluster}', alpha=0.5)

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(
            frameon=False,
            bbox_to_anchor=(1.05, 1),
            loc='upper left',
            borderaxespad=0.
        )
        plt.tight_layout()
        plt.savefig(args.output_figure)
        plt.close()

    if args.output_dataframe is not None:
        embedding_df.to_csv(args.output_dataframe, index_label="strain")

def calculate_distances_from_alignment(alignment_path, indel_distance):
    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(alignment_path, "fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = list(sequences_by_name.keys())

    hamming_distances = get_hamming_distances(
        list(sequences_by_name.values()),
        indel_distance,
    )
    distance_matrix = pd.DataFrame(squareform(hamming_distances))
    distance_matrix.index = sequence_names
    distance_matrix.columns = distance_matrix.index

    return distance_matrix

def distance(args):
    distance_matrix = calculate_distances_from_alignment(
        args.alignment,
        args.indel_distance,
    )
    distance_matrix.to_csv(args.output)

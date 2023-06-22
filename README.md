# pathogen-embed
Create reduced dimension embeddings for pathogen sequences


pathogen-embed is an open-source software for scientists, epidemiologists, etc. to run reduced dimension embeddings (PCA, MDS, t-SNE, and UMAP) on viral populations. This is the source code from the paper Cartography written by Sravani Nanduri and John Huddleston.

[Documentation](https://blab.github.io/pathogen-embed/)\
[Source Code](https://github.com/blab/pathogen-embed/tree/main)\
[Bug reports](https://github.com/blab/pathogen-embed/issues)

# Build documentation  
Build the [Documentation](https://blab.github.io/pathogen-embed/):

``` make -C /docs html ```

Clean the docs.

``` make -C /docs clean ```

# Releasing a new version

### Information about each file

#### README.md

contains the description of the package pathogen-embed.

#### setup.py

Gives PyPi the instructions about where to find dependent packages, the authors and relevant links, etc. Also gives the entry points for the console script, which tells Pypi to call the main function of __main__.py. 

#### __init__.py

Initializes the package, creates the parser to parse the command line arguments and pass them into the embed.py function.

#### __main__.py

Calls the "run" function in __init__.py, which calls embed.py. 

#### embed.py

The main code for the package.

# To create new version 

Run 

``` python3 -m build ``` 

This creates the dist folder that gets uploaded to pypi.

``` python3 -m twine upload dist/* ```

Input the username and password, upload new dist files to pypi. Make sure the version of the dist folders does not already exist within pypi. 


## Installing the package

Simply install the package using pip.

```
pip install pathogen-embed
```

# src.embed module

## Command line interface

The full [Documentation](https://blab.github.io/pathogen-embed/). 

The below documentation does not detail the named and positional arguments. 

Reduced dimension embeddings for pathogen sequences


```
usage: embed [-h] [--distance-matrix DISTANCE_MATRIX] [--separator SEPARATOR]
             [--alignment ALIGNMENT] [--cluster-data CLUSTER_DATA]
             [--cluster-threshold CLUSTER_THRESHOLD]
             [--random-seed RANDOM_SEED] [--output-dataframe OUTPUT_DATAFRAME]
             [--output-figure OUTPUT_FIGURE]
             {pca,t-sne,umap,mds} ...
```

### Sub-commands:

#### pca

Principal Component Analysis

```
embed pca [-h] [--components COMPONENTS]
          [--explained-variance EXPLAINED_VARIANCE]
```

#### t-sne

t-distributed Stochastic Neighborhood Embedding

```
embed t-sne [-h] [--perplexity PERPLEXITY] [--learning-rate LEARNING_RATE]
```

#### umap

Uniform Manifold Approximation and Projection

```
embed umap [-h] [--nearest-neighbors NEAREST_NEIGHBORS] [--min-dist MIN_DIST]
```

#### mds

Multidimensional Scaling

```
embed mds [-h] [--components COMPONENTS]
```

## API


### src.embed.get_hamming_distances(genomes)
Calculate pairwise Hamming distances between the given list of genomes
and return the nonredundant array of values for use with scipy’s squareform function.
Bases other than standard nucleotides (A, T, C, G) are ignored.


* **Parameters**

    **genomes** (*list*) – a list of strings corresponding to genomes that should be compared



* **Returns**

    a list of distinct Hamming distances as a vector-form distance vector



* **Return type**

    list


```python
>>> genomes = ["ATGCT", "ATGCT", "ACGCT"]
>>> get_hamming_distances(genomes)
[0, 1, 1]
>>> genomes = ["AT-GCT", "AT--CT", "AC--CT"]
>>> get_hamming_distances(genomes)
[0, 1, 1]
```

# Issues and fixes:

Issue/Fix: Used sphinx-book-theme version 0.3.3 for backwards compatibility (wouldn't render otherwise)

Issue: index.rst: Module "src" has no attribute "make_parser"
Incorrect argparse :module: or :func: values?

Fix: changed module from src to src.embed





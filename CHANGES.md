# CHANGELOG

## 3.1.0

### Features

* pathogen-cluster-mutations: Add this new top-level command to allow users to create a table of mutations that appear in clusters or other previously defined genetic groups ([#36][])

[#36]: https://github.com/blab/pathogen-embed/pull/36

## 3.0.0

### Major Changes

* pathogen-embed: Use "simplex" encoding by default for PCA ([#35][])
* pathogen-cluster: Change default minimum number of samples for a cluster from 5 to 10 ([#35][])

[#35]: https://github.com/blab/pathogen-embed/pull/35

## 2.3.0

### Features

* pathogen-cluster: Add `--distance-matrix` input argument to support HDBSCAN clustering of genetic distances from `pathogen-distance` ([#33][])

[#33]: https://github.com/blab/pathogen-embed/pull/33

## 2.2.1

### Bug Fixes

* Do not internally sort embedding inputs by sequence name ([#32][])
* Let scikit-learn automatically pick SVD algorithm to use for PCA instead of hardcoding the "full" solver ([#31][])

[#31]: https://github.com/blab/pathogen-embed/pull/31
[#32]: https://github.com/blab/pathogen-embed/pull/32

## 2.2.0

### Features

* admin: Publish to PyPI with GitHub Actions ([#30][])

### Bug Fixes

* Use inferred types for external embedding parameters ([#29][])

[#29]: https://github.com/blab/pathogen-embed/pull/29
[#30]: https://github.com/blab/pathogen-embed/pull/30

## 2.1.0

### Features

* Add alternate encodings of nucleotide sequences for PCA in `pathogen-embed` ([#23][])

[#23]: https://github.com/blab/pathogen-embed/pull/23

## 2.0.0

### Major Changes

* Removed seaborn as a dependency in favor of base matplotlib ([#13][])
* Set default learning rate for t-SNE to "auto" such that the learning rate scales with the sample size ([#12][])

### Features

* Add support for multiple alignment and/or distance matrix inputs to `pathogen-embed` ([#19][])
* Add optional output from `pathogen-embed` that produces the boxplot figure of Euclidean by genetic distance ([#14][])

### Bug Fixes

* Display default parameters for subcommands of pathogen-embed ([#12][])

[#12]: https://github.com/blab/pathogen-embed/pull/12
[#13]: https://github.com/blab/pathogen-embed/pull/13
[#14]: https://github.com/blab/pathogen-embed/pull/14
[#19]: https://github.com/blab/pathogen-embed/pull/19

## 1.1.2

### Bug Fixes

* Fix t-SNE keyword argument error associated with recent versions of scikit-learn ([#6][])
* Pass random seed argument from the command line to PCA and MDS implementations ([#6][])

[#6]: https://github.com/blab/pathogen-embed/pull/6

## 1.1.1

### Bug Fixes

* Fix MDS stress output

## 1.1.0

### Features
* Add stress value dataframe to MDS arguments to relay fitness of the embedding

## 1.0.0

### Major Changes

* Created separate commands for embedding, clustering, and distance matrix creation (pathogen-embed, pathogen-cluster, pathogen-distance) ([#2](https://github.com/blab/pathogen-embed/pull/2))

## 0.1.1

### Bug Fixes

* Migrated source code and documentation from [cartography](https://github.com/blab/cartography) to [pathogen-embed](https://github.com/blab/pathogen-embed)

## 0.1.0

### Features

* Only calculate the distance matrix for an alignment if it isn't available already ([194fd74](https://github.com/blab/cartography/commit/194fd746c458d51bb73c962728da6c242a2d00f0))
* Source embedding params from cluster data ([8c26898](https://github.com/blab/cartography/commit/8c268981fa20d59888a92c0f38eedf8b42065db8))
* Initialize t-SNE embeddings with PCA instead of a random initialization ([89fc458](https://github.com/blab/cartography/commit/89fc4583e9af3caab332405dbb3b1f9e10f06c29))

### Bug Fixes

* Avoid re-reading alignment for PCA and t-SNE ([5fb2bbb](https://github.com/blab/cartography/commit/5fb2bbb13d686660dae48d1744f42a69c18fea57))

## 0.0.2

### Bug Fixes

* Issue [on github](https://github.com/blab/cartography/issues/20). The hamming distance calculation
now also works on lowercase fasta files as well as uppercase.


## 0.0.1

* First version of embed-pathogen.

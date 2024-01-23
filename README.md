# pathogen-embed

Create reduced dimension embeddings for pathogen sequences

pathogen-embed is an open-source software for scientists, epidemiologists, etc. to run reduced dimension embeddings (PCA, MDS, t-SNE, and UMAP) on viral populations. This is the source code from the paper Cartography written by Sravani Nanduri and John Huddleston.

 - [Change log](CHANGES.md)
 - [Documentation](https://blab.github.io/pathogen-embed/)
 - [Source Code](https://github.com/blab/pathogen-embed/tree/main)
 - [Bug reports](https://github.com/blab/pathogen-embed/issues)

## Installation

Simply install the package using pip.

```
pip install pathogen-embed
```

## Build documentation

Build the [Documentation](https://blab.github.io/pathogen-embed/):

```
make -C /docs html
```

Clean the docs.

```
make -C /docs clean
```

## Releasing a new version

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

### To create new version

Run

``` python3 -m build ```

This creates the dist folder that gets uploaded to pypi.

``` python3 -m twine upload dist/* ```

Input the username and password, upload new dist files to pypi. Make sure the version of the dist folders does not already exist within pypi.

## Run tests

Run tests with cram.

``` bash
cram tests
```

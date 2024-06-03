Try to run pathogen-embed with PCA and no alignment input.
This should fail with an error.

  $ pathogen-embed \
  >   --output-dataframe embed_pca.csv \
  >   pca
  ERROR: PCA requires an alignment input to create the embedding.
  [1]

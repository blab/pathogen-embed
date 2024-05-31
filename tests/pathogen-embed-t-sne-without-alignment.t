Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output distances.csv

Try to run t-SNE with a distance matrix but no alignment input.
This should fail with an error.

  $ pathogen-embed \
  >   --distance-matrix distances.csv \
  >   --output-dataframe embed_t-sne.csv \
  >   t-sne
  ERROR: t-SNE requires an alignment input to initialize the embedding.
  [1]

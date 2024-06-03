Always use the simplest Numba threading layer, to avoid lower-level OMP warnings.

  $ export NUMBA_THREADING_LAYER=workqueue

Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output distances.csv

Run pathogen-embed with UMAP on a H3N2 HA alignment with the distance matrix.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --distance-matrix distances.csv \
  >   --output-dataframe embed_umap.csv \
  >   umap \
  >   --min-dist 0.1 \
  >   --nearest-neighbors 25

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_umap.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

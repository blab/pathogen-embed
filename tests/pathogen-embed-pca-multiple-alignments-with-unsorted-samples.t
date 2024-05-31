Run pathogen-embed with PCA on a H3N2 HA and H3N2 NA alignments with sequence names in a different order in each file.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output-dataframe embed_pca_unsorted.csv \
  >   pca \
  >    --components 2

Run PCA again but with alignments with sequence names in the same order in each file.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --output-dataframe embed_pca_sorted.csv \
  >   pca \
  >    --components 2

The two PCA embeddings should be identical regardless of the sorting of sequences in the input alignments.

  $ diff -u embed_pca_unsorted.csv embed_pca_sorted.csv

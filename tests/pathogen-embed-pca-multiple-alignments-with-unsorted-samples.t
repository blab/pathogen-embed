Run pathogen-embed with PCA on a H3N2 HA and H3N2 NA alignments with sequence names in a different order in each file.
This should cause an error explaining how to properly sort alignments prior to the embedding.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output-dataframe embed_pca_unsorted.csv \
  >   pca \
  >    --components 2
  ERROR: The given alignments do not have the same sequence names in the same order. Confirm your alignments have the same sequence names and sort your alignments (e.g., `seqkit sort -n alignment.fasta > sorted_alignment.fasta`) so they have the same order.
  [1]

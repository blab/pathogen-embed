Create copies of two input alignments with different number of samples in each.

  $ cp $TESTDIR/data/h3n2_ha_alignment.fasta .
  $ echo ">sample\nACGT" > h3n2_na_alignment.fasta

Run pathogen-embed with PCA on a H3N2 HA and H3N2 NA alignments.

  $ pathogen-embed \
  >   --alignment h3n2_ha_alignment.fasta h3n2_na_alignment.fasta \
  >   --output-dataframe embed_pca.csv \
  >   pca
  ERROR: The given alignments do not have the same sequence names in the same order. Confirm your alignments have the same sequence names and sort your alignments (e.g., `seqkit sort -n alignment.fasta > sorted_alignment.fasta`) so they have the same order.
  [1]

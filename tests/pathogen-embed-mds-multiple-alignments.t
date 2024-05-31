Run pathogen-embed with MDS on a H3N2 HA and H3N2 NA alignments such that distance matrices get calculated on the fly.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output-dataframe embed_mds.csv \
  >   mds \
  >   --components 2

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_mds.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

Repeat the embedding with sorted alignment inputs.
The resulting embeddings should be identical, since the script should internally sort distance matrices by sequence names.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --output-dataframe embed_mds_sorted.csv \
  >   mds \
  >   --components 2
  $ diff -u embed_mds.csv embed_mds_sorted.csv

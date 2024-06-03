Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output ha_distances.csv

Get a distance matrix from a H3N2 NA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output na_distances.csv

Run pathogen-embed with t-SNE on distances from H3N2 HA and H3N2 NA alignments.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --distance-matrix ha_distances.csv na_distances.csv \
  >   --output-dataframe embed_t-sne.csv \
  >   t-sne

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_t-sne.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

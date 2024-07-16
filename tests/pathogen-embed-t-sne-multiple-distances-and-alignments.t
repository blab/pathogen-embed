Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta \
  >   --output ha_distances.csv

Get a distance matrix from a H3N2 NA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --output na_distances.csv

Run pathogen-embed with t-SNE on distances from H3N2 HA and H3N2 NA alignments.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --distance-matrix ha_distances.csv na_distances.csv \
  >   --output-dataframe embed_t-sne.csv \
  >   --output-figure embed_t-sne.pdf \
  >   --output-pairwise-distance-figure t-sne_pairwise_distances.pdf \
  >   t-sne

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_t-sne.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.sorted.fasta | wc -l) ]]

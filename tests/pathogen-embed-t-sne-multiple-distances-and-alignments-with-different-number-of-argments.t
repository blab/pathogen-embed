Get a distance matrix from a H3N2 HA alignment that has been sorted by sequence name.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta \
  >   --output ha_distances.csv

Get a distance matrix from a H3N2 NA alignment that has been sorted by sequence name.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --output na_distances.csv

Run pathogen-embed with t-SNE on distances from H3N2 HA and H3N2 NA alignments.
Only provide one alignment input, though, which should cause an error.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta \
  >   --distance-matrix ha_distances.csv na_distances.csv \
  >   --output-dataframe embed_t-sne.csv \
  >   t-sne
  ERROR: If giving multiple alignments and distance matrices the number of both must match.
  [1]

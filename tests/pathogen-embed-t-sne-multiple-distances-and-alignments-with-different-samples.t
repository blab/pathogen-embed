Get a distance matrix from a H3N2 HA alignment that has been sorted by sequence name.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta \
  >   --output ha_distances_complete.csv

Get a distance matrix from a H3N2 NA alignment that has been sorted by sequence name.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --output na_distances_complete.csv

Remove the second record from the HA and NA distance matrices.
This should produce mismatched records between the alignments and distances, but the pairs of alignments and distances on their own are matched.

  $ cut -f 1,3- -d "," ha_distances_complete.csv | sed 2d > ha_distances.csv
  $ cut -f 1,3- -d "," na_distances_complete.csv | sed 2d > na_distances.csv

Run pathogen-embed with t-SNE on distances from H3N2 HA and H3N2 NA alignments.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --distance-matrix ha_distances.csv na_distances.csv \
  >   --output-dataframe embed_t-sne.csv \
  >   t-sne
  ERROR: The sequence names for the distance matrix inputs do not match the names in the alignment inputs.
  [1]

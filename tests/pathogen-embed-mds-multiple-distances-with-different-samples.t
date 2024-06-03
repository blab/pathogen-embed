Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output ha_distances.csv

Get a distance matrix from a H3N2 NA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output na_distances_complete.csv

Remove the second record from the NA distance matrix, so we have mismatched records between HA and NA.

  $ cut -f 1,3- -d "," na_distances_complete.csv | sed 2d > na_distances.csv

Run pathogen-embed with MDS on mismatched distances with different samples.

  $ pathogen-embed \
  >   --distance-matrix ha_distances.csv na_distances.csv \
  >   --output-dataframe embed_mds.csv \
  >   mds \
  >   --components 2
  ERROR: The given distance matrices do not have the same sequence names.
  [1]

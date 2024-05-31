Get a distance matrix from a H3N2 NA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output na_distances_complete.csv

Remove the second row from the NA distance matrix, so we have a non-square distance matrix.

  $ cut -f 1,3- -d "," na_distances_complete.csv > na_distances.csv

Try to run pathogen-embed with MDS on a non-square distance matrix.
This should fail.

  $ pathogen-embed \
  >   --distance-matrix na_distances.csv \
  >   --output-dataframe embed_mds.csv \
  >   mds \
  >   --components 2
  ERROR: Distance matrices must be square (with the same number of rows and columns).
  [1]

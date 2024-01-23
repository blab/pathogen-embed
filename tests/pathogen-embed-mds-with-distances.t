Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output distances.csv

Run pathogen-embed with MDS on a H3N2 HA alignment with the distance matrix.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --distance-matrix distances.csv \
  >   --output-dataframe embed_mds.csv \
  >   mds \
  >   --components 2 \
  >   --stress embed_mds_stress.csv

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_mds.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output ha_distances.csv

Get a distance matrix from a H3N2 NA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output na_distances.csv

Run pathogen-embed with MDS on distances from H3N2 HA and H3N2 NA alignments.

  $ pathogen-embed \
  >   --distance-matrix ha_distances.csv na_distances.csv \
  >   --output-dataframe embed_mds.csv \
  >   mds \
  >   --components 2

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_mds.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

The order of records in the embedding should be alphabetically sorted and contain the same sequence names as the input distances.

  $ cut -f 1 -d "," ha_distances.csv | sed 1d | sort -k 1,1 > ha_distances_records.txt
  $ cut -f 1 -d "," embed_mds.csv | sed 1d > embed_records.txt
  $ diff -u ha_distances_records.txt embed_records.txt

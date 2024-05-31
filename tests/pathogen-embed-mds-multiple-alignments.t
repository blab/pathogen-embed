Run pathogen-embed with MDS on a H3N2 HA and H3N2 NA alignments.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta $TESTDIR/data/h3n2_na_alignment.fasta \
  >   --output-dataframe embed_mds.csv \
  >   mds \
  >   --components 2 \
  >   --stress embed_mds_stress.csv

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_mds.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

There should be 1 entry for the stress of the MDS embedding.

  $ wc -l < embed_mds_stress.csv
  \s*1 (re)


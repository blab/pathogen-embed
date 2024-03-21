Run pathogen-embed with MDS on a H3N2 HA alignment.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_mds.csv \
  >   --output-figure embed_mds_fig.png \
  >   mds \
  >   --components 2 \
  >   --stress embed_mds_stress.csv
  $ [[ -e embed_mds_fig.png ]]
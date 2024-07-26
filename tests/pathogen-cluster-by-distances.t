Get a distance matrix from a H3N2 HA alignment.

  $ pathogen-distance \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta \
  >   --output ha_distances.csv

Find clusters from the genetic distances.

  $ pathogen-cluster \
  >   --distance-matrix ha_distances.csv \
  >   --label-attribute genetic_label \
  >   --distance-threshold 0.5 \
  >   --output-dataframe cluster_distances.csv

There should be one record in the cluster output for each record in the input distances.

  $ [[ $(wc -l < cluster_distances.csv) == $(wc -l < ha_distances.csv) ]]

The header should include the index strain, each strain from the distance matrix, and the requested cluster label.

  $ head -n 1 cluster_distances.csv
  strain,.*,genetic_label (re)

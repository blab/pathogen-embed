Run pathogen-embed with PCA on a H3N2 HA alignment.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta \
  >   --output-dataframe embed.csv \
  >   pca \
  >   --components 2

Find clusters from the embedding.

  $ pathogen-cluster \
  >   --embedding embed.csv \
  >   --label-attribute cluster_label \
  >   --distance-threshold 0.5 \
  >   --output-dataframe cluster_embed.csv

Find mutations per cluster.

  $ pathogen-cluster-mutations \
  >   --reference-sequence $TESTDIR/data/h3n2_ha_reference.fasta \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta \
  >   --clusters cluster_embed.csv \
  >   --cluster-column cluster_label \
  >   --min-allele-count 10 \
  >   --min-allele-frequency 0.5 \
  >   --output mutations_cluster_embed.csv

Confirm that the mutation table output has the correct structure and more than one row.

  $ head -n 1 mutations_cluster_embed.csv
  mutation,cluster_count,distinct_clusters,cluster_column

  $ [[ $(sed 1d mutations_cluster_embed.csv | wc -l | sed 's/ //g') > 0 ]]

Run pathogen-embed with PCA on a H3N2 HA alignment.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_pca.csv \
  >   pca \
  >   --components 2

Find clusters in the PCA embedding.

  $ pathogen-cluster \
  >   --embedding embed_pca.csv \
  >   --label-attribute pca_label \
  >   --distance-threshold 0.5 \
  >   --output-dataframe cluster_embed_pca.csv

There should be one record in the cluster output for each record in the input embedding.

  $ [[ $(wc -l < cluster_embed_pca.csv) == $(wc -l < embed_pca.csv) ]]

The header should include the strain, the PCA components, and the requested cluster label.

  $ head -n 1 cluster_embed_pca.csv
  strain,pca1,pca2,pca_label

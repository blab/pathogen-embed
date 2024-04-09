Run pathogen cluster and check that error message is thrown when non csv input is given for embedding

  $ pathogen-cluster \
  >   --embedding embed_pca.tsv \
  >   --label-attribute pca_label \
  >   --distance-threshold 0.5 \
  >   --output-dataframe cluster_embed_pca.csv

== stderr
You must supply a CSV file for the embedding.
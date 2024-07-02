Run pathogen-embed with t-SNE on a H3N2 HA alignment.
Use the simplex encoding of the alignment for the PCA initialization.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_t-sne.csv \
  >   t-sne \
  >   --perplexity 25 \
  >   --learning-rate 100 \
  >   --pca-encoding simplex

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_t-sne.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

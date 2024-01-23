Run pathogen-embed with t-SNE on a H3N2 HA alignment.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_t-sne.csv \
  >   t-sne \
  >   --perplexity 100 \
  >   --learning-rate 100

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_t-sne.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

Run pathogen-embed with t-SNE on a H3N2 HA alignment with parameters loaded from a file.
This should print changes from the defaults to the parameters in the file to the screen.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --embedding-parameters $TESTDIR/data/h3n2_ha_t-sne_parameters.csv \
  >   --output-dataframe embed_t-sne.csv \
  >   t-sne
  INFO: Replacing embedding parameter perplexity value of '30.0' with '45.0' provided by '(.*)'. (re)
  INFO: Replacing embedding parameter learning_rate value of 'auto' with '100.0' provided by '(.*)'. (re)

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_t-sne.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

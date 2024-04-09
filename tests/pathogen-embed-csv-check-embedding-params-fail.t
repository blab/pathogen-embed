Run pathogen embed and check that embedding parameters with not csv input fails

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_t-sne.tsv \
  >   --output-figure figure.png \
  >   --embedding-parameters value.tsv \
  >   t-sne \
  >   --perplexity 25

== stderr
You must supply a CSV file for embedding parameters.
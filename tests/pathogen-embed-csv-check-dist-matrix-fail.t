Run pathogen-embed with a tsv instead of a csv to see if it throws an error.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_t-sne.tsv \
  >   t-sne \
  >   --perplexity 25
You must supply a CSV file for distance_matrix.
[1]

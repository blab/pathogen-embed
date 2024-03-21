Run pathogen-embed with t-SNE on a H3N2 HA alignment, specifically looking for the pairwise distance PNG

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_t-sne.csv \
  >   --output-pairwise-distance-figure t-sne_pairwise_output.png \
  >   t-sne \
  >   --perplexity 25 \
  >   --learning-rate 100

There should be one record in the embedding per input sequence in the alignment.

  $ [[ -e t-sne_pairwise_output.png ]]

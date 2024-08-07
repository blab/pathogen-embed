Run pathogen-embed with PCA on a H3N2 HA and H3N2 NA alignments.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.sorted.fasta $TESTDIR/data/h3n2_na_alignment.sorted.fasta \
  >   --output-dataframe embed_pca.csv \
  >   pca \
  >   --components 2 \
  >   --explained-variance embed_pca_variance.csv

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_pca.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.sorted.fasta | wc -l) ]]

There should be 2 components of variance explained.

  $ sed 1d embed_pca_variance.csv | wc -l
  \s*2 (re)

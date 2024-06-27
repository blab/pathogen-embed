Run pathogen-embed with PCA on a H3N2 HA alignment using the default "integer" encoding.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_pca_integer.csv \
  >   pca \
  >   --components 2 \
  >   --encoding integer

Run pathogen-embed with PCA on a H3N2 HA alignment using the "genotype" encoding.

  $ pathogen-embed \
  >   --alignment $TESTDIR/data/h3n2_ha_alignment.fasta \
  >   --output-dataframe embed_pca_genotype.csv \
  >   pca \
  >   --components 2 \
  >   --encoding genotype

There should be one record in the embedding per input sequence in the alignment.

  $ [[ $(sed 1d embed_pca_integer.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]
  $ [[ $(sed 1d embed_pca_genotype.csv | wc -l) == $(grep "^>" $TESTDIR/data/h3n2_ha_alignment.fasta | wc -l) ]]

The two PCA outputs should differ because of their different input encoding.

  $ diff -q embed_pca_integer.csv embed_pca_genotype.csv
  Files embed_pca_integer.csv and embed_pca_genotype.csv differ
  [1]

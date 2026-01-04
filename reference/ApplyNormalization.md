# Applies normalization on sequencing data

Applies normalization on sequencing data

## Usage

``` r
ApplyNormalization(
  data,
  scaling.method = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
  ruv.norm = TRUE,
  ruv.k = 1,
  ruv.drop = 0,
  control.idx = NULL,
  sc.idx = NULL,
  enrich.idx = NULL,
  spike.in.prefix = NULL,
  synthetic.id = NULL
)
```

## Arguments

- data:

  A un-normalized count data matrix of shape n x p, where n is the
  number of samples and p is the number of features.

- scaling.method:

  Vector of normalization methods that are applied to the data.
  Available methods are: `c("TC", "UQ", "TMM", "DESeq", "PossionSeq")`.
  Select one or multiple methods. By default all normalization methods
  will be applied.

- ruv.norm:

  Whether to perform RUV normalization.

- ruv.k:

  The number of factors of unwanted variation to be estimated from the
  data.

- ruv.drop:

  The number of singular values to drop in the estimation of unwanted
  variation, default drop the first singular value that represent the
  difference between enrichment and input.

- control.idx:

  Vector of the negative control genes for RUV normalization.

- sc.idx:

  A numeric matrix specifying the replicate samples for which to compute
  the count differences used to estimate the factors of unwanted
  variation.

- enrich.idx:

  Matrix with two rows indicating the column index of enrichment and
  input samples in the raw/normalized count data matrix. The first row
  is the column index of input and the second row is the column index of
  enrichment samples.

- spike.in.prefix:

  A character specify the prefix of spike-in id.

- synthetic.id:

  Character or vector of string specifying the name of synthetic RNAs.

## Value

List of objects containing normalized data and associated normalization
factors.

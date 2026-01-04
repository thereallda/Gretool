# Perform RUV normalization

Perform RUV normalization

## Usage

``` r
normRUV(
  data,
  control.idx = NULL,
  sc.idx = NULL,
  method = c("RUVg", "RUVs", "RUVse"),
  k = 1,
  drop = 0,
  log = TRUE
)
```

## Arguments

- data:

  A un-normalized count data matrix of shape n x p, where n is the
  number of samples and p is the number of features.

- control.idx:

  Vector of control genes' id.

- sc.idx:

  A numeric matrix specifying the replicate samples for which to compute
  the count differences used to estimate the factors of unwanted
  variation.

- method:

  Perform RUVg or RUVs normalization.

- k:

  The number of factors of unwanted variation to be estimated from the
  data.

- drop:

  The number of singular values to drop in the estimation of unwanted
  variation, default drop the first singular value that represent the
  difference between enrichment and input.

- log:

  Whether to perform log2-transformation with 1 offset on data matrix
  (default: TRUE), while normalized counts are returned in non-log
  format.

## Value

List containing normalized counts and adjust factors for adjusting
unwanted variation.

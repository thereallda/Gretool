# Remove unwanted variation using control genes

Adapted from `RUVSeq` (https://doi.org/10.1038%2Fnbt.2931)

## Usage

``` r
enRUVg(object, log = TRUE, k = 1, tolerance = 1e-08, control.idx, drop = 1)
```

## Arguments

- object:

  A counts matrix.

- log:

  Whether to perform log2-transformation with 1 offset on data matrix,
  default: TRUE.

- k:

  The number of factors of unwanted variation to be estimated from the
  data.

- tolerance:

  Tolerance in the selection of the number of positive singular values,
  i.e., a singular value must be larger than tolerance to be considered
  positive.

- control.idx:

  ID of control genes.

- drop:

  The number of singular values to drop in the estimation of unwanted
  variation, default drop the first singular value that represent the
  difference between enrichment and input

## Value

list contain a matrix of unwanted factors (W) and corrected counts
matrix (normalizedCounts).

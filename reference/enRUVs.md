# Remove unwanted variation using control genes within replicates

Adapted from `RUVSeq` (https://doi.org/10.1038%2Fnbt.2931)

## Usage

``` r
enRUVs(
  object,
  log = TRUE,
  k = 2,
  tolerance = 1e-08,
  control.idx,
  sc.idx,
  drop = 1
)
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

- sc.idx:

  A numeric matrix specifying the replicate samples for which to compute
  the count differences used to estimate the factors of unwanted
  variation (see details).

- drop:

  The number of singular values to drop in the estimation of unwanted
  variation, default drop the first singular value that represent the
  difference between enrichment and input

## Value

list contain a matrix of unwanted factors (W) and corrected counts
matrix (normalizedCounts).

## Details

Each row of sc.idx should correspond to a set of replicate samples. The
number of columns is the size of the largest set of replicates; rows for
smaller sets are padded with -1 values.

For example, if the sets of replicate samples are (1,2,3),(4,5),(6,7,8),
then scIdx should be 1 2 3 4 5 -1 6 7 8

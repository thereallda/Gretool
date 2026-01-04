# Perform PossionSeq normalization

Perform PossionSeq normalization

## Usage

``` r
normPossionSeq(data, ...)
```

## Arguments

- data:

  A un-normalized count data matrix of shape n x p, where n is the
  number of samples and p is the number of features.

- ...:

  Additional parameters can be passed to `PS.Est.Depth`.

## Value

List containing normalized counts and normalized factors for library
size.

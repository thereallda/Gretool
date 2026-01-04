# Esimate sequencing depth

Functions from PossionSeq (https://doi.org/10.1093/biostatistics/kxr031)

## Usage

``` r
PS.Est.Depth(n, iter = 5, ct.sum = 5, ct.mean = 0.5)
```

## Arguments

- n:

  The data matrix. The rows are counts for a gene, and the columns are
  counts from an experiment.

- iter:

  Number of iterations used. Default value: 5. The default value is
  usually a good choice.

- ct.sum:

  if the total number of reads of a gene across all experiments \<=
  ct.sum, this gene will not be considered for estimating sequencing
  depth. Default value: 5.

- ct.mean:

  if the mean number of reads of a gene across all experiments \<=
  ct.mean, this gene will not be considered for estimating sequencing
  depth. Default value: 0.5.

## Value

estimated sequencing depth. a vector. their product is 1.

## Examples

``` r
mat <- matrix(rnbinom(300, mu=100, size=1), ncol=5)
PS.Est.Depth(mat)
#> [1] 1.3508303 1.0316168 0.9758812 0.8248881 0.8914332
```

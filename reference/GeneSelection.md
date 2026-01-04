# Gene set selection

Select negative control genes for RUV normalization, positive and
negative evaluation genes for assessment.

## Usage

``` r
GeneSelection(object, n.neg.control = 1000, n.pos.eval = 500, n.neg.eval = 500)
```

## Arguments

- object:

  Gretool object

- n.neg.control:

  Number of negative control genes for RUV normalization, default: 1000.

- n.pos.eval:

  Number of positive evaluation genes for wanted variation assessment,
  default: 500.

- n.neg.eval:

  Number of negative evaluation genes for unwanted variation assessment,
  default: 500.

## Value

list of genes

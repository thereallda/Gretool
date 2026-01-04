# Filter Low Expressed Genes

Filter Low Expressed Genes

## Usage

``` r
FilterLowExprGene(x, group = NULL, min.count = 10)
```

## Arguments

- x:

  Input data, can be counts matrix or object.

- group:

  Vector or factor giving group membership for a oneway layout, if
  appropriate, default: NULL.

- min.count:

  Minimum count required for at least some samples, default: 10.

## Value

Filtered count matrix or updated object.

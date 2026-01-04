# Accessor of gene set

Get gene set.

## Usage

``` r
getGeneSet(object, name = c("NegControl", "NegEvaluation", "PosEvaluation"))

# S4 method for class 'Gretool,character'
getGeneSet(object, name = c("NegControl", "NegEvaluation", "PosEvaluation"))
```

## Arguments

- object:

  Gretool.

- name:

  Name of the gene set, one of `NegControl`, `NegEvaluation`, or
  `PosEvaluation`

## Value

Vector of gene ID

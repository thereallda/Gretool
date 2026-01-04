# Combine list of DE results

Combine list of DE results

## Usage

``` r
reduceRes(res.ls, logfc.col, levels = names(res.ls))
```

## Arguments

- res.ls:

  Named list of differential analysis results tables. Each elements in
  the list correspond to a table of differential analysis results
  between two groups of samples.

- logfc.col:

  Column name of the log fold-change.

- levels:

  Factor levels of the groups, default order by the element order of
  `res.ls`.

## Value

data.frame

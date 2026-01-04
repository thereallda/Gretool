# Wrapper of edgeR procedure

Wrapper of edgeR procedure

## Usage

``` r
edgeRDE(
  counts,
  group,
  norm.factors = NULL,
  adjust.factors = NULL,
  design.formula = NULL,
  contrast.df = NULL,
  coef = NULL,
  logfc.cutoff = 1,
  p.cutoff = 0.05,
  prior.count = 0.125,
  only.pos = TRUE
)
```

## Arguments

- counts:

  A un-normalized counts data matrix.

- group:

  Vector of length p mapping the columns of `counts` to corresponding
  samples group.

- norm.factors:

  Vector of normalization factors with p length.

- adjust.factors:

  Matrix with each column indicates the adjusting factors that estimated
  from RUV.

- design.formula:

  Formula

- contrast.df:

  Data frame of contrast, where extracting results as first column vs.
  second column.

- coef:

  Integer or character vector indicating which coefficients of the
  linear model are to be tested equal to zero. Values must be columns or
  column names of design.

- logfc.cutoff:

  Filter genes by at least X-fold difference (log2-scale) between the
  two groups of samples, default: 1.

- p.cutoff:

  Filter genes by no more than Y adjusted p-value, default: 0.05.

- prior.count:

  Average prior count to be added to observation to shrink the estimated
  log-fold-changes towards zero, default: 0.125.

- only.pos:

  Only return positive genes in filtered results `res.sig.ls`, default:
  TRUE.

## Value

List containing differential analysis object, result table and filtered
result table.

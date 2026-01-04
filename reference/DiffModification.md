# Differential modification by linear model test

Differential modification by linear model test

## Usage

``` r
DiffModification(
  object,
  col.data,
  contrast.df,
  norm.method = NULL,
  test.gene.id = NULL,
  assay.colname = "assay",
  biology.colname = "bio",
  input.id = "Input",
  enrich.id = "Enrich",
  norm = FALSE,
  adjust = TRUE,
  p.cutoff = 0.05,
  pseudo.count = 0.125,
  eps = 1e-04,
  shrink_lambda = 0
)
```

## Arguments

- object:

- col.data:

  `data.frame` with at least three columns (indicate condition, enrich
  and biological groups). Rows of `col.data` correspond to columns of
  `data`.

- contrast.df:

  Data frame of contrast, where extracting results as first column vs.
  second column.

- norm.method:

  Normalization methods to use

- test.gene.id:

  List of gene ids to be tested

- assay.colname:

  Name of the assay column in `col.data`, e.g., "enrich".

- biology.colname:

  Name of the biological column in `col.data`, e.g., "bio".

- input.id:

  Input library id, must be consistent with the enrich column of
  `col.data`, e.g., "Input".

- enrich.id:

  Enrich library id, must be consistent with the enrich column of
  `col.data`, e.g., "Enrich".

- adjust:

  Whether to perform adjustment based on input expression levels.

- p.cutoff:

  P values Cutoff for significant differential modification genes.

- pseudo.count:

  Minimal offset to prevent negative values.

- eps:

  Threshold to prevent negative values, default: 1e-4

- shrink_lambda:

  Factor to control the shrinkage of expression size factor from input,
  default: 0.

## Value

List of differential modification results

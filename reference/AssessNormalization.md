# Assessments of normalization performance

Assessments of normalization performance

## Usage

``` r
AssessNormalization(
  data.ls,
  bio.group = NULL,
  enrich.group = NULL,
  batch.group = NULL,
  eval.pam.k = 2:6,
  eval.pc.n = 3,
  log = TRUE,
  pos.eval.set = NULL,
  neg.eval.set = NULL
)
```

## Arguments

- data.ls:

  List containing normalized counts and adjust factors for adjusting
  unwanted variation. Output of `ApplyNormalization`.

- bio.group:

  Vector of index indicating the column index of samples of each
  biological groups in the raw/normalized count data matrix.

- enrich.group:

  Vector of index indicating the column index of enrichment and input
  samples in the raw/normalized count data matrix.

- batch.group:

  Vector of index indicating the column index of each batch groups in
  the raw/normalized count data matrix.

- eval.pam.k:

  Integer or vector of integers indicates the number of clusters for PAM
  clustering in performance evaluation, default: 2:6.

- eval.pc.n:

  Integer indicates the evaluation metrics will be calculated in the
  first nth PCs, default: 3.

- log:

  Whether to perform log2-transformation with 1 offset on data matrix,
  default: TRUE.

- pos.eval.set:

  Vector of genes id.

- neg.eval.set:

  Vector of genes id.

## Value

List containing the metrics matrix and the ranking matrix, both sorted
by the score of methods from top to bottom.

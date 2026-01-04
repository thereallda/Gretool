# Gretool object and constructor

`Gretool` object extends the `SummarizedExperiment` class. The
`createGretool` is a easy constructor of `Gretool` object

## Usage

``` r
createGretool(
  data,
  col.data,
  spike.in.prefix = NULL,
  input.id = "Input",
  enrich.id = "Enrich",
  synthetic.id = NULL
)
```

## Arguments

- data:

  A un-normalized count data matrix of shape n x p, where n is the
  number of samples and p is the number of features.

- col.data:

  `data.frame` with at least two columns (indicate condition and enrich
  groups). Rows of `col.data` correspond to columns of `data`.

- spike.in.prefix:

  A character specify the prefix of spike-in id, e.g., "^FB" stands for
  fly spike-in id, default: NULL.

- input.id:

  Input library id, must be consistent with the enrich column of
  `col.data`, e.g., "Input".

- enrich.id:

  Enrich library id, must be consistent with the enrich column of
  `col.data`, e.g., "Enrich".

- synthetic.id:

  Vector of synthetic RNA id, e.g. c("Syn1","Syn2"), default: NULL.

## Value

Gretool object

## Details

Description of each slot:  
`assay`
[`SummarizedExperiment::Assays`](https://rdrr.io/pkg/SummarizedExperiment/man/Assays-class.html)
object, contains all counts.  
`counts` list for holding raw/normalized counts of sample and
spike_in.  
`norm_factors` list for holding normalization factors of sample and
spike_in.  
`norm_metrics` data frame with normalization methods in row and metrics
in columns.  
`norm_score` data frame with normalization methods in row and scores in
columns.  
`enrichment` list for holding all differential analysis results.  
`enrichment_filtered` lists for holding filtered differential analysis
results.  
`parameter` list of parameters.  

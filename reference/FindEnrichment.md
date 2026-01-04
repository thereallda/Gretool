# Find enriched genes between enrich and input samples

Find enriched genes between enrich and input samples

## Usage

``` r
FindEnrichment(
  object,
  slot = c("sample", "spike_in"),
  norm.method,
  logfc.cutoff = 1,
  p.cutoff = 0.05,
  only.pos = TRUE,
  ...
)
```

## Arguments

- object:

  object.

- slot:

  Which slot, one of `sample` or `spike_in`.

- norm.method:

  Which normalization method to use for differential analysis, must be
  one of the methods presented in the selected slot.

- logfc.cutoff:

  Filter genes by at least X-fold difference (log2-scale) between the
  two groups of samples, default: 1.

- p.cutoff:

  Filter genes by no more than `p.cutoff` adjusted p-value, default:
  0.05.

- only.pos:

  Whether significant differential genes only consider genes with
  positive fold-change.

- ...:

  Additional parameters can be passed to `edgeRDE`.

## Value

updated object

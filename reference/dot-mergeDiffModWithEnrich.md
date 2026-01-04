# Merge differential modification data with enrichment data

Merge differential modification data with enrichment data

## Usage

``` r
.mergeDiffModWithEnrich(
  diffmod_res,
  enrich_list,
  contrast_name,
  fc.thresh = log2(1.2),
  p.thresh = 0.05
)
```

## Arguments

- diffmod_res:

  Table of differential modification results from `DiffModification`.

- enrich_list:

  Enriched genes from `FindEnrichment`

- contrast_name:

  Contrast name, e.g., "A_B".

- fc.thresh:

  Fold-change threshold for significant differential modification.

- p.thresh:

  P values threshold for significant differential modification.

## Value

list

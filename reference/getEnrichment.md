# Accessor of enrichment results

Get all or filtered enrichment

## Usage

``` r
getEnrichment(object, slot = c("sample", "spike_in"), filter = FALSE)

# S4 method for class 'Gretool,character'
getEnrichment(object, slot = c("sample", "spike_in"), filter = FALSE)
```

## Arguments

- object:

  Gretool.

- slot:

  Which slot to get, one of `sample` or `spike_in`.

- filter:

  Whether to get the filtered enrichment, default FALSE.

## Value

list of enrichment table

# Outlier Test

Rosner’s outlier test on principal component 1 to assess and remove
potential outliers.

## Usage

``` r
OutlierTest(x, return = FALSE, remove = FALSE)
```

## Arguments

- x:

  object

- return:

  Whether to return object or simply perform test, default FALSE

- remove:

  Whether to remove outliers. If TRUE, must be paired with
  "return=TRUE", default: FALSE

## Value

updated object

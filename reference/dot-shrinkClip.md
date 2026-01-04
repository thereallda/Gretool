# Shrink factors toward 1 and clip bounds

Shrink factors toward 1 and clip bounds

## Usage

``` r
.shrinkClip(
  size_factor,
  shrink_lambda = 0,
  min_factor = 1e-06,
  max_factor = 1e+06
)
```

## Arguments

- size_factor:

  from `.expressionFactor`.

- shrink_lambda:

  see description in `DiffModification`.

- min_factor:

  minimum factor to clip, default: 1e-6.

- max_factor:

  maximum factor to clip, default: 1e6.

## Value

shrunk expression factors

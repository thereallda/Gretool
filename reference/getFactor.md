# Accessor of normalization factors

Accessor of normalization factors

## Usage

``` r
getFactor(object, slot = c("sample", "spike_in"), method)

# S4 method for class 'Gretool,character,character'
getFactor(object, slot = c("sample", "spike_in"), method)
```

## Arguments

- object:

  Gretool.

- slot:

  Which slot to get, one of `sample` or `spike_in`.

- method:

  Which normalization methods to get, must be one of the methods
  presented in the selected slot.

## Value

vector or list of factors.

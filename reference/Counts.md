# Accessor for counts

Accessor for counts

## Usage

``` r
Counts(object, slot = c("sample", "spike_in"), method)

Counts(object, slot = c("sample", "spike_in"), method) <- value

# S4 method for class 'Gretool,character,character'
Counts(object, slot = c("sample", "spike_in"), method)

# S4 method for class 'Gretool,character,character,matrix'
Counts(object, slot = c("sample", "spike_in"), method) <- value
```

## Arguments

- object:

  Gretool.

- slot:

  Which slot to get, one of `sample` or `spike_in`.

- method:

  Which counts matrix to get, must be one of the raw or normalized
  counts matrix presented in the selected slot.

- value:

  Raw or normalized counts matrix.

## Value

matrix.

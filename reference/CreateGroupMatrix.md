# Create a matrix for RUVSeq

Create a matrix for RUVSeq

## Usage

``` r
CreateGroupMatrix(group.vec)
```

## Arguments

- group.vec:

  A vector indicating membership in a group.

## Value

A matrix.

## Examples

``` r
CreateGroupMatrix(c('a','b','b','c','c','c','a','d','d'))
#>      [,1] [,2] [,3]
#> [1,]    1    7   -1
#> [2,]    2    3   -1
#> [3,]    4    5    6
#> [4,]    8    9   -1
```

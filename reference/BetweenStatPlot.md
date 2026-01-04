# Box-violin plot comparing values between groups

Box-violin plot comparing values between groups

## Usage

``` r
BetweenStatPlot(
  data,
  x,
  y,
  color,
  palette = NULL,
  test = c("wilcox.test", "t.test", "none"),
  add.p = c("p", "p.adj"),
  comparisons = NULL,
  step.increase = 0.3,
  title = NULL
)
```

## Arguments

- data:

  A data frame (or a tibble).

- x:

  The grouping variable from the `data`.

- y:

  The value variable from the `data`.

- color:

  The color variable from the `data`.

- palette:

  The color palette for different groups.

- test:

  Perform wilcoxon rank sum test or t-test or no test, must be one of
  c("wilcox.test", "t.test", "none").

- add.p:

  Label p-value or adjusted p-value, must be one of c("p", "p.adj").

- comparisons:

  A list of length-2 vectors specifying the groups of interest to be
  compared. For example to compare groups "A" vs "B" and "B" vs "C", the
  argument is as follow: comparisons = list(c("A", "B"), c("B", "C"))

- step.increase:

  Numeric vector with the increase in fraction of total height for every
  additional comparison to minimize overlap.

- title:

  Plot title, default NULL.

## Value

ggplot2 object

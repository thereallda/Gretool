# PCA plot from counts matrix

PCA plot from counts matrix

## Usage

``` r
PCAplot(
  object,
  use.pc = c(1, 2),
  color = NULL,
  label = NULL,
  shape = NULL,
  title = NULL,
  vst.norm = FALSE,
  palette = NULL,
  repel = TRUE
)
```

## Arguments

- object:

  A count matrix.

- use.pc:

  Which two PCs to be used, default PC1 in x-axis and PC2 in y-axis.

- color:

  Vector indicates the color mapping of samples, default NULL.

- label:

  Vector of sample names or labels, default NULL.

- shape:

  Vector indicates the shape mapping of samples, default NULL.

- title:

  Plot title, default NULL.

- vst.norm:

  Whether to perform `vst` transformation, default FALSE.

- palette:

  The color palette for different groups.

- repel:

  Whether to use `ggrepel` to avoid overlapping text labels or not,
  default TRUE.

## Value

ggplot2 object

## Details

Perform PCA based on matrix using `prcomp`, and visualized with scatter
plot.

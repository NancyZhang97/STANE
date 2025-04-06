
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STANE

<!-- badges: start -->
<!-- badges: end -->

The goal of STANE (Spatial Transcriptomics Analysis considering Niche
Effects) is to estimate cell type-specific proportions in spot-level
spatial transcriptomics data as well as gene expression alterations
caused by niche effects through cell-cell interactions. The STANE
applies an iterated computation process which use estimated niche
effects to adjust altered gene expression in different spots, in order
to achieve a more accurate results.

## Installation

You can install the development version of STANE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("NancyZhang97/STANE")
```

## Example

In the following sections, we will use example dataset to

``` r
library(STANE)
## basic example code
```

## Reference

\[1\] Heiland D, Kueckelhaus J (2024). *SPATAData: The MILOLab Spatial
Transriptomic Database*. R package version 0.0.0.9000.

\[2\] Hao et al. Dictionary learning for integrative, multimodal and
scalable single-cell analysis. Nature Biotechnology (2023) \[Seurat V5\]

\[3\] Cilluffo G, Sottile G, La Grutta S, Muggeo V (2020). “The Induced
Smoothed lasso: A practical framework for hypothesis testing in high
dimensional regression.” *Statistical Methods in Medical Research*,
*29*(3), 765-777. <doi:10.1177/0962280219842890>
<https://doi.org/10.1177/0962280219842890>.

\[4\] Friedman J, Tibshirani R, Hastie T (2010). “Regularization Paths
for Generalized Linear Models via Coordinate Descent.” *Journal of
Statistical Software*, *33*(1), 1-22. <doi:10.18637/jss.v033.i01>
<https://doi.org/10.18637/jss.v033.i01>.

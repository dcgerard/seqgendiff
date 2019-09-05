
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RNA-Seq Generation/Modification for Simulation

[![Travis-CI Build
Status](https://travis-ci.org/dcgerard/seqgendiff.svg?branch=master)](https://travis-ci.org/dcgerard/seqgendiff)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/dcgerard/seqgendiff?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/seqgendiff)
[![Coverage
Status](https://img.shields.io/codecov/c/github/dcgerard/seqgendiff/master.svg)](https://codecov.io/github/dcgerard/seqgendiff?branch=master)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![CRAN
status](https://www.r-pkg.org/badges/version/seqgendiff)](https://cran.r-project.org/package=seqgendiff)
[![](http://cranlogs.r-pkg.org/badges/grand-total/seqgendiff)](https://cran.r-project.org/package=seqgendiff)

This package will take real RNA-seq data (either single-cell or bulk)
and alter it by adding signal to it. This signal is in the form of a
generalized linear model with a log (base-2) link function under a
Poisson / negative binomial / mixture of negative binomials
distribution. The advantage of this way of simulating data is that you
can see how your method behaves when the simulated data exhibit common
(and annoying) features of real data. This is without you having to
specify these features *a priori*. We call the way we add signal
“binomial thinning”.

The main functions are:

  - `select_counts`: Subsample the columns and rows of a real RNA-seq
    count matrix. You would then feed this sub-matrix into one of the
    thinning functions below.
  - `thin_diff`: The function most users should be using for
    general-purpose binomial thinning. For the special applications of
    the two-group model or library/gene thinning, see the functions
    listed below.
  - `thin_2group`: The specific application of thinning in the two-group
    model.
  - `thin_lib`: The specific application of library size thinning.
  - `thin_gene`: The specific application of total gene expression
    thinning.
  - `thin_all`: The specific application of thinning all counts.
  - `effective_cor`: Returns an estimate of the actual correlation
    between the surrogate variables and a user-specified design matrix.
  - `ThinDataToSummarizedExperiment`: Converts a `ThinData` object to a
    `SummarizedExperiment` object.
  - `ThinDataToDESeqDataSet`: Converts a `ThinData` object to a
    `DESeqDataSet` object.

If you find a bug or want a new feature, please submit an
[issue](http://github.com/dcgerard/seqgendiff/issues).

Check out [NEWS](NEWS.md) for updates.

# Installation

To install from CRAN, run the following code in R:

``` r
install.packages("seqgendiff")
```

To install the latest version of seqgendiff, run the following code in
R:

``` r
install.packages("devtools")
devtools::install_github("dcgerard/seqgendiff")
```

To get started, check out the vignettes by running the following in R:

``` r
library(seqgendiff)
browseVignettes(package = "seqgendiff")
```

Or you can check out the vignettes I post online:
<https://dcgerard.github.io/seqgendiff/>.

# Citation

If you use this package, please cite:

> Gerard D (2019). “Data-based RNA-seq Simulations by Binomial
> Thinning.” *bioRxiv*. doi:
> [10.1101/758524](https://doi.org/10.1101/758524).

A BibTeX entry for LaTeX users is

``` tex
@article{,
    author = {Gerard, David},
    title = {Data-based {RNA}-seq Simulations by Binomial Thinning},
    elocation-id = {758524},
    year = {2019},
    doi = {10.1101/758524},
    publisher = {Cold Spring Harbor Laboratory},
    journal = {bioRxiv}
}
```

# Code of Conduct

Please note that the ‘seqgendiff’ project is released with a
[Contributor Code of
Conduct](https://github.com/dcgerard/seqgendiff/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.


<!-- README.md is generated from README.Rmd. Please edit that file -->
Sequence Generation for Differential Expression Analysis
========================================================

[![Travis-CI Build Status](https://travis-ci.org/dcgerard/seqgendiff.svg?branch=master)](https://travis-ci.org/dcgerard/seqgendiff) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dcgerard/seqgendiff?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/seqgendiff) [![Coverage Status](https://img.shields.io/codecov/c/github/dcgerard/seqgendiff/master.svg)](https://codecov.io/github/dcgerard/seqgendiff?branch=master) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

This package is designed to take real RNA-seq data and alter it by adding signal to a known proportion of the genes. The advantage of this way of simulating data is that you can see how your method behaves when the simulated data exhibit common (and annoying) features of real data. For example, in the real world data are not normally distributed and unobserved confounding is a major issue. If your sample size is set small enough (i.e. under realistic settings), then this package will simulate data that exhibit these characteristics.

Installation
============

To install, run the following code in R:

``` r
install.packages("devtools")
devtools::install_github("dcgerard/seqgendiff")
```

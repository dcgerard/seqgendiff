#' seqgendiff: Sequence Generation/Modification for Differential Expression Analysis Simulations.
#'
#' This package is designed to take real RNA-seq data and alter it by
#' adding a known amount of signal. You can then use this modified dataset
#' in simulation studies for differential expression analysis. The
#' advantage of this way of simulating data is that you can see how
#' your method behaves when the simulated data exhibit common
#' (and annoying) features of real data. For example, in the real world
#' data are not normally (or negative binomially) distributed and
#' unobserved confounding is a major issue. This package will simulate
#' data that exhibit these characteristics.
#'
#' @section seqgendiff Functions:
#' \describe{
#'   \item{\code{\link{thin_diff}}}{For the function most users should
#'       be using for general-purpose Poisson thinning. For the special
#'       applications of the two-group model or library/gene thinning, see
#'       the functions listed below.}
#'   \item{\code{\link{thin_2group}}}{For the specific application of
#'       thinning in the two-group model.}
#'   \item{\code{\link{thin_lib}}}{For the specific application of
#'       library size thinning.}
#'   \item{\code{\link{thin_gene}}}{For the specific application of
#'       total gene expression thinning.}
#'   \item{\code{\link{ThinDataToSummarizedExperiment}}}{For converting a
#'       ThinData object to a SummarizedExperiment object.}
#'   \item{\code{\link{ThinDataToDESeqDataSet}}}{For converting a
#'       ThinData object to a DESeqDataSet object.}
#' }
#'
#' @docType package
#' @name seqgendiff
#'
#' @author David Gerard
NULL

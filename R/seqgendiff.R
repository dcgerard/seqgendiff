#' seqgendiff: RNA-Seq Generation/Modification for Simulation
#'
#' This package is designed to take real RNA-seq data and alter it by
#' adding a known amount of signal. You can then use this modified dataset
#' in simulation studies for differential expression analysis, factor
#' analysis, confounder adjustment, or library size adjustment. The
#' advantage of this way of simulating data is that you can see how
#' your method behaves when the simulated data exhibit common
#' (and annoying) features of real data. For example, in the real world
#' data are not normally (or negative binomially) distributed and
#' unobserved confounding is a major issue. This package will simulate
#' data that exhibit these characteristics. The methods used in this
#' package are described in detail in Gerard (2020).
#'
#' @section seqgendiff Functions:
#' \describe{
#'   \item{\code{\link{select_counts}()}}{Subsample the columns and rows
#'       of a real RNA-seq count matrix. You would then feed this sub-matrix
#'       into one of the thinning functions below.}
#'   \item{\code{\link{thin_diff}()}}{The function most users should
#'       be using for general-purpose binomial thinning. For the special
#'       applications of the two-group model or library/gene thinning, see
#'       the functions listed below.}
#'   \item{\code{\link{thin_2group}()}}{The specific application of
#'       thinning in the two-group model.}
#'   \item{\code{\link{thin_lib}()}}{The specific application of
#'       library size thinning.}
#'   \item{\code{\link{thin_gene}()}}{The specific application of
#'       total gene expression thinning.}
#'   \item{\code{\link{thin_all}()}}{The specific application of thinning
#'       all counts.}
#'   \item{\code{\link{effective_cor}()}}{Returns an estimate of the actual
#'       correlation between surrogate variables and a user-specified
#'       design matrix.}
#'   \item{\code{\link{ThinDataToSummarizedExperiment}()}}{Converts a
#'       ThinData object to a SummarizedExperiment object.}
#'   \item{\code{\link{ThinDataToDESeqDataSet}()}}{Converts a
#'       ThinData object to a DESeqDataSet object.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gerard, D (2020). "Data-based RNA-seq simulations by binomial thinning." \emph{BMC Bioinformatics}. 21(1), 206. \doi{10.1186/s12859-020-3450-9}.}
#' }
#'
#' @docType package
#' @name seqgendiff-package
#' @aliases seqgendiff
#'
#' @author David Gerard
"_PACKAGE"

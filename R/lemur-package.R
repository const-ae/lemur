
#'
#' @import methods
#' @import stats
#' @rawNamespace import(SingleCellExperiment, except = weights)
#' @importFrom S4Vectors metadata `metadata<-` DataFrame
#' @importFrom SummarizedExperiment assayNames assay `assay<-` colData `colData<-` rowData `rowData<-`
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom BiocGenerics design
#' @importFrom utils .DollarNames capture.output head txtProgressBar setTxtProgressBar
#' @importFrom Matrix t
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

## usethis namespace: start
#' @useDynLib lemur, .registration = TRUE
## usethis namespace: end
NULL

# Satisfy NOTE from CRAN CMD check
# Both are used in 'find_de_neighborhoods'
utils::globalVariables(c("cntrst", "..did_indicator"))

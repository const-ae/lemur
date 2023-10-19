#' The `lemur_fit` class
#'
#' The `lemur_fit` class extends [`SingleCellExperiment`] and provides
#' additional accessors to get the values of the values produced by [`lemur`].
#'
#' @param object the `lemur_fit` object for the [`BiocGenerics::design`] generic
#' @param x,i,j,...,drop the `lemur_fit` object and indices for the `[` subsetting operator
#'
#' @details
#'
#' To access the values produced by [`lemur`], use the dollar notation (`$`):
#' \describe{
#'  \item{`fit$n_embedding`}{the number of embedding dimensions.}
#'  \item{`fit$design`}{the specification of the design in [`lemur`]. Usually this is a [`stats::formula`].}
#'  \item{`fit$base_point`}{a matrix (`nrow(fit) * fit$n_embedding`) with the base point for the Grassmann exponential map.}
#'  \item{`fit$coefficients`}{a three-dimensional tensor (`nrow(fit) * fit$n_embedding * ncol(fit$design_matrix)`) with the coefficients for
#'    the exponential map.}
#'  \item{`fit$embedding`}{a matrix (`fit$n_embedding * ncol(fit)`) with the low dimensional position for each cell.}
#'  \item{`fit$design_matrix`}{a matrix with covariates for each cell (`ncol(fit) * ncol(fit$design_matrix)`).}
#'  \item{`fit$linear_coefficients`}{a matrix (`nrow(fit) * ncol(fit$design_matrix)`) with the coefficients for the linear regression.}
#'  \item{`fit$alignment_coefficients`}{a 3D tensor with the coefficients for the alignment (`fit$n_embedding * fit$n_embedding * ncol(fit$design_matrix)`)}
#'  \item{`fit$alignment_design`}{an alternative design specification for the alignment. This is typically a [`stats::formula`].}
#'  \item{`fit$alignment_design_matrix`}{an alternative design matrix specification for the alignment.}
#'  \item{`fit$contrast`}{a parsed version of the contrast specification from the `test_de` function or `NULL`.}
#'  \item{`fit$colData`}{the column annotation `DataFrame`.}
#'  \item{`fit$rowData`}{the row annotation `DataFrame`.}
#' }
#'
#' @seealso [`lemur`], [`predict`][predict.lemur_fit], [`residuals`][residuals,lemur_fit-method]
#'
#' @returns An object of class `lemur_fit`.
#'
#' @rdname lemur_fit
#' @aliases lemur_fit
#'
#' @examples
#' # The easiest way to make a lemur_fit object, is to call `lemur`
#' data(glioblastoma_example_data)
#' fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition,
#'              n_emb = 5, verbose = FALSE)
#'
#' fit$n_embedding
#' fit$embedding[,1:10]
#' fit$n_embedding
#' fit$embedding[,1:10]
#' fit$design_matrix[1:10,]
#' fit$coefficients[1:3,,]
#'
#' @export
.lemur_fit <- setClass("lemur_fit", contains = "SingleCellExperiment")

lemur_fit <- function(data, col_data, row_data,
                      n_embedding,
                      design, design_matrix, linear_coefficients,
                      base_point, coefficients, embedding,
                      alignment_coefficients,
                      alignment_design, alignment_design_matrix,
                      use_assay, is_test_data){
  if(is.null(data)){
    data <- SingleCellExperiment(assays = list())
  }else if(! is(data, "SummarizedExperiment")){
    data <- SingleCellExperiment(assays = setNames(list(data), use_assay))
  }else{
    data <- as(data, "SingleCellExperiment")
  }

  n_features <- nrow(data)
  linearFit <- LinearEmbeddingMatrix(design_matrix, linear_coefficients)
  colnames(linearFit) <- colnames(linear_coefficients)

  # Do some more graceful merging of the existing and the new information.
  colData(data) <- (colData(data) %update_values% col_data) |> as("DataFrame")
  SingleCellExperiment::int_colData(data)$is_test_data <- is_test_data
  rowData(data) <- (rowData(data) %update_values% row_data) |> as("DataFrame")
  reducedDims(data) <- reducedDims(data) %update_values% list(linearFit = linearFit, embedding = t(embedding))
  metadata(data) <- metadata(data) %update_values% list(
    n_embedding = n_embedding, design = design, base_point = base_point,
    coefficients = coefficients, alignment_coefficients = alignment_coefficients,
    alignment_design = alignment_design, alignment_design_matrix = alignment_design_matrix,
    use_assay = use_assay, row_mask = rep(TRUE, n_features))

  .lemur_fit(data)
}


S4Vectors::setValidity2("lemur_fit", function(obj){
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  msg <- NULL

  n_features_original <- length(metadata(obj)$row_mask)
  n_features <- nrow(obj)
  n_obs <- ncol(obj)

  n_embedding <-  obj$n_embedding
  if(is.null(n_embedding)) msg <- c(msg, "'n_embedding' must not be NULL")
  base_point <- obj$base_point
  if(is.null(base_point)) msg <- c(msg, "'base_point' must not be NULL")
  coefficients <- obj$coefficients
  if(is.null(coefficients)) msg <- c(msg, "'coefficients' must not be NULL")
  alignment_coefficients <- obj$alignment_coefficients
  if(is.null(alignment_coefficients)) msg <- c(msg, "'alignment_coefficients' must not be NULL")
  alignment_design <- obj$alignment_design
  alignment_design_matrix <- obj$alignment_design_matrix
  if(is.null(alignment_design_matrix)) msg <- c(msg, "'alignment_design_matrix' must not be NULL")
  embedding <- obj$embedding
  if(is.null(embedding)) msg <- c(msg, "'embedding' must not be NULL")
  design <- obj$design
  design_matrix <- sampleFactors(reducedDim(obj, "linearFit"))
  if(is.null(design_matrix)) msg <- c(msg, "'design_matrix' must not be NULL")
  linear_coefficients <- featureLoadings(reducedDim(obj, "linearFit"))
  if(is.null(linear_coefficients)) msg <- c(msg, "'linear_coefficients' must not be NULL")
  use_assay <- obj$use_assay
  if(is.null(use_assay)) msg <- c(msg, "'use_assay' must not be NULL")
  is_test_data <- obj$is_test_data
  if(is.null(is_test_data)) msg <- c(msg, "'is_test_data' must not be NULL")
  col_names <- colnames(obj)
  row_names <- rownames(obj)

  if(! is.null(design_matrix) && nrow(design_matrix) != n_obs) msg <- c(msg, "`nrow(design_matrix)` does not match number of observations")
  if(! is.null(linear_coefficients) && nrow(linear_coefficients) != n_features_original) msg <- c(msg, "`nrow(linear_coefficients)` does not match `n_features_original`")
  if(! is.null(linear_coefficients) && ncol(linear_coefficients) != ncol(design_matrix)) msg <- c(msg, "`ncol(linear_coefficients)` does not match `ncol(design_matrix)`")
  if(! is.null(base_point) && nrow(base_point) != n_features_original) msg <- c(msg, "`nrow(base_point)` does not match `n_features_original`")
  if(! is.null(base_point) && ncol(base_point) != n_embedding) msg <- c(msg, "`ncol(base_point)` does not match `n_embedding`")
  if(! is.null(coefficients) && ! is.array(coefficients) || length(dim(coefficients)) != 3) msg <- c(msg, "`coefficients` must be a three dimensional array")
  if(! is.null(coefficients) && dim(coefficients)[1] != n_features_original) msg <- c(msg, "`dim(coefficients)[1]` does not match `n_features_original`")
  if(! is.null(coefficients) && dim(coefficients)[2] != n_embedding) msg <- c(msg, "`dim(coefficients)[2]` does not match `n_embedding`")
  if(! is.null(coefficients) && dim(coefficients)[3] != ncol(design_matrix)) msg <- c(msg, "`dim(coefficients)[3]` does not match `ncol(design_matrix)`")
  if(! is.null(embedding) && nrow(embedding) != n_embedding) msg <- c(msg, "`nrow(embedding)` does not match `n_embedding`")
  if(! is.null(embedding) && ncol(embedding) != n_obs) msg <- c(msg, "`ncol(embedding)` does not match number of observations")
  if(! is.null(alignment_coefficients) && ! is.array(alignment_coefficients) || length(dim(alignment_coefficients)) != 3) msg <- c(msg, "`alignment_coefficients` must be a three dimensional array")
  if(! is.null(alignment_coefficients) && dim(alignment_coefficients)[1] != n_embedding) msg <- c(msg, "`dim(alignment_coefficients)[1]` does not match `n_embedding`")
  if(! is.null(alignment_coefficients) && dim(alignment_coefficients)[2] != n_embedding + 1) msg <- c(msg, "`dim(alignment_coefficients)[2]` does not match `n_embedding + 1`")
  if(! is.null(alignment_coefficients) && dim(alignment_coefficients)[3] != ncol(alignment_design_matrix)) msg <- c(msg, "`dim(alignment_coefficients)[3]` does not match `ncol(alignment_design_matrix)`")
  if(! is.null(alignment_design_matrix) && nrow(alignment_design_matrix) != n_obs) msg <- c(msg, "`nrow(alignment_design_matrix)` does not match number of observations")
  if(! is.null(alignment_design) &&  ! inherits(alignment_design, "formula")) msg <- c(msg, "`alignment_design` must inherit from formula or be NULL")
  if(! is.null(design) &&  ! inherits(design, "formula")) msg <- c(msg, "`design` must inherit from formula or be NULL")
  if(! is.null(use_assay) &&  ! use_assay %in% assayNames(obj)) msg <- c(msg, "`use_assay` must be one of the assays")
  if(! is.null(is_test_data) &&  ! is.logical(is_test_data)) msg <- c(msg, "`is_test_data` must be a logical vector")
  if(! is.null(is_test_data) &&  length(is_test_data) != n_obs ) msg <- c(msg, "length `is_test_data` must match the number of observations")
  if(! is.null(col_names) && length(col_names) != length(unique(col_names))) msg <- c(msg, "`colnames` are not unique")
  if(! is.null(row_names) && length(row_names) != length(unique(row_names))) msg <- c(msg, "`rownames` are not unique")

  if(is.null(msg)){
    TRUE
  }else{
    msg
  }
})

# Subsetting

#' @rdname lemur_fit
#' @export
setMethod("[", c("lemur_fit", "ANY", "ANY"), function(x, i, j, ...) {
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  i_missing <- missing(i)
  j_missing <- missing(j)

  if (! i_missing) {
    # Update metadata
    ii <- convert_subset_to_index(i, rownames(x))
    old_mask <- metadata(x)$row_mask
    metadata(x)$row_mask[] <- FALSE
    metadata(x)$row_mask[old_mask][ii] <- TRUE
  }
  if(! j_missing){
    jj <- convert_subset_to_index(j, colnames(x))
    metadata(x)[["alignment_design_matrix"]] <- metadata(x)[["alignment_design_matrix"]][jj,,drop=FALSE]
  }

  callNextMethod()
})

convert_subset_to_index <- function(subset, names){
  if (is.character(subset)) {
    orig <- subset
    subset <- match(subset, names)
    if (any(bad <- is.na(subset))) {
      bad_examples <- toString(orig[bad], width = 100)
      stop("index out of bounds: ", bad_examples)
    }
  }
  return(as.vector(subset))
}

.methods_to_suggest <- c("n_embedding", "embedding",
                         "design", "design_matrix", "base_point",
                         "coefficients", "linear_coefficients", "alignment_coefficients",
                         "alignment_design", "alignment_design_matrix",
                         "contrast", "use_assay", "colData", "rowData",
                         "test_data", "training_data", "is_test_data")


#' @rdname lemur_fit
#' @export
setMethod("design", signature = "lemur_fit", function(object){
  metadata(object)[["design"]]
})



#' @rdname cash-lemur_fit-method
#' @export
.DollarNames.lemur_fit <- function(x, pattern = ""){
  grep(pattern, .methods_to_suggest, value = TRUE)
}

#' Access values from a `lemur_fit`
#'
#' @param x the `lemur_fit`
#' @param pattern the pattern from looking up potential values interactively
#' @param name the name of the value behind the dollar
#' @param value the replacement value. This only works for `colData` and
#'   `rowData`.
#'
#' @returns The respective value stored in the `lemur_fit` object.
#'
#' @seealso [`lemur_fit-class`] for more documentation on the
#'   accessor functions.
#' @aliases dollar_methods
setMethod("$", "lemur_fit",
          function(x, name){
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  if(! name %in% .methods_to_suggest){
    stop("Illegal name after '$' sign: ", name)
  }
  switch(name,
    n_embedding =             metadata(x)[["n_embedding"]],
    design =                  design(x),
    base_point =              metadata(x)[["base_point"]],
    coefficients =            metadata(x)[["coefficients"]],
    embedding =               t(reducedDim(x, "embedding")),
    design_matrix =           sampleFactors(reducedDim(x, "linearFit")),
    linear_coefficients =     featureLoadings(reducedDim(x, "linearFit")),
    alignment_design =        metadata(x)[["alignment_design"]],
    alignment_design_matrix = metadata(x)[["alignment_design_matrix"]],
    alignment_coefficients =  metadata(x)[["alignment_coefficients"]],
    contrast =                metadata(x)[["contrast"]],
    use_assay =               metadata(x)[["use_assay"]],
    colData =                 colData(x),
    rowData =                 rowData(x),
    test_data =               get_test_data(x),
    training_data =           get_training_data(x),
    is_test_data =            int_colData(x)[["is_test_data"]],
    stop("Invalid `name` value.")
  )
})

get_test_data <- function(fit){
  fit[,fit$is_test_data]
}

get_training_data <- function(fit){
  fit[,!fit$is_test_data]
}

#' @rdname cash-lemur_fit-method
setReplaceMethod("$", "lemur_fit",
                 function(x, name, value){
  if(! name %in% .methods_to_suggest){
   stop("Illegal name after '$' sign: ", name)
  }
  switch(name,
        colData = {SummarizedExperiment::colData(x) <- value},
        rowData = {SummarizedExperiment::rowData(x) <- value},
        stop("It is illegal to modify the content of lemur_fit object")
  )
  x
})



















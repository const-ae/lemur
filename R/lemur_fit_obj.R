
#' @export
.lemur_fit_obj <- setClass("lemur_fit_obj", contains = "SingleCellExperiment")

lemur_fit_obj <- function(data_mat, col_data, row_data,
                       n_ambient, n_embedding,
                       ambient_coordsystem, ambient_offset,
                       design, design_matrix, linear_coefficients,
                       base_point, coefficients, embedding,
                       alignment_method, alignment_rotation, alignment_stretching,
                       alignment_design, alignment_design_matrix){


  if(is.null(data_mat)){
    n_features <- nrow(row_data)
    assays <- list()
  }else{
    n_features <- nrow(data_mat)
    assays <- list(expr = data_mat)
  }
  linearFit <- LinearEmbeddingMatrix(design_matrix, linear_coefficients)
  colnames(linearFit) <- colnames(linear_coefficients)
  sce <- SingleCellExperiment(assays, colData = col_data, rowData = row_data,
                              reducedDims = list(linearFit = linearFit,
                                                 embedding = t(embedding)),
                              metadata = list(n_ambient = n_ambient, n_embedding = n_embedding,
                                              ambient_coordsystem = ambient_coordsystem, ambient_offset = ambient_offset,
                                              design = design,
                                              base_point = base_point, coefficients = coefficients,
                                              alignment_method = alignment_method, alignment_rotation = alignment_rotation, alignment_stretching = alignment_stretching,
                                              alignment_design = alignment_design, alignment_design_matrix = alignment_design_matrix))
  .lemur_fit_obj(sce)
}


S4Vectors::setValidity2("lemur_fit_obj", function(obj){
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  msg <- NULL

  if(length(assayNames(obj)) > 0 && assayNames(obj)[1] != "expr"){
    msg <- c(msg, "'expr' must be the first assay")
  }

  # n_features <- sum(int_metadata(obj)$row_mask)
  n_features <- nrow(obj)
  n_obs <- ncol(obj)

  n_ambient <-  obj$n_ambient
  eff_n_ambient <- if(is.infinite(n_ambient)){
    n_features
  }else{
    n_ambient
  }
  if(is.null(n_ambient)) msg <- c(msg, "'n_ambient' must not be NULL")
  n_embedding <-  obj$n_embedding
  if(is.null(n_embedding)) msg <- c(msg, "'n_embedding' must not be NULL")
  ambient_coordsystem <- obj$ambient_coordsystem
  if(is.null(ambient_coordsystem)) msg <- c(msg, "'ambient_coordsystem' must not be NULL")
  ambient_offset <- obj$ambient_offset
  if(is.null(ambient_offset)) msg <- c(msg, "'ambient_offset' must not be NULL")
  base_point <- obj$base_point
  if(is.null(base_point)) msg <- c(msg, "'base_point' must not be NULL")
  coefficients <- obj$coefficients
  if(is.null(coefficients)) msg <- c(msg, "'coefficients' must not be NULL")
  alignment_method <- obj$alignment_method
  alignment_rotation <- obj$alignment_rotation
  if(is.null(alignment_rotation)) msg <- c(msg, "'alignment_rotation' must not be NULL")
  alignment_stretching <- obj$alignment_stretching
  if(is.null(alignment_stretching)) msg <- c(msg, "'alignment_stretching' must not be NULL")
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

  if(! is.null(ambient_coordsystem) && nrow(ambient_coordsystem) != n_features) msg <- c(msg, "`nrow(ambient_coordsystem)` does not match the number of features")
  if(! is.null(ambient_coordsystem) && ncol(ambient_coordsystem) != eff_n_ambient) msg <- c(msg, "`ncol(ambient_coordsystem)` does not match `n_ambient`")
  if(! is.null(ambient_offset) && length(ambient_offset) != n_features) msg <- c(msg, "`length(ambient_offset)` does not match number of features")
  if(! is.null(design_matrix) && nrow(design_matrix) != n_obs) msg <- c(msg, "`nrow(design_matrix)` does not match number of observations")
  if(! is.null(linear_coefficients) && nrow(linear_coefficients) != eff_n_ambient) msg <- c(msg, "`nrow(linear_coefficients)` does not match `n_ambient`")
  if(! is.null(linear_coefficients) && ncol(linear_coefficients) != ncol(design_matrix)) msg <- c(msg, "`ncol(linear_coefficients)` does not match `ncol(design_matrix)`")
  if(! is.null(base_point) && nrow(base_point) != eff_n_ambient) msg <- c(msg, "`nrow(base_point)` does not match `n_ambient`")
  if(! is.null(base_point) && ncol(base_point) != n_embedding) msg <- c(msg, "`ncol(base_point)` does not match `n_embedding`")
  if(! is.null(coefficients) && ! is.array(coefficients) || length(dim(coefficients)) != 3) msg <- c(msg, "`coefficients` must be a three dimensional array")
  if(! is.null(coefficients) && dim(coefficients)[1] != eff_n_ambient) msg <- c(msg, "`dim(coefficients)[1]` does not match `n_ambient`")
  if(! is.null(coefficients) && dim(coefficients)[2] != n_embedding) msg <- c(msg, "`dim(coefficients)[2]` does not match `n_embedding`")
  if(! is.null(coefficients) && dim(coefficients)[3] != ncol(design_matrix)) msg <- c(msg, "`dim(coefficients)[3]` does not match `ncol(design_matrix)`")
  if(! is.null(embedding) && nrow(embedding) != n_embedding) msg <- c(msg, "`nrow(embedding)` does not match `n_embedding`")
  if(! is.null(embedding) && ncol(embedding) != n_obs) msg <- c(msg, "`ncol(embedding)` does not match number of observations")
  if(! is.null(alignment_method) && length(alignment_method) != 1 && length(alignment_method) != n_obs) msg <- c(msg, "`length(alignment_method)` must either be a single value or a vector with the same length as the number of observation")
  if(! is.null(alignment_rotation) && ! is.array(alignment_rotation) || length(dim(alignment_rotation)) != 3) msg <- c(msg, "`alignment_rotation` must be a three dimensional array")
  if(! is.null(alignment_rotation) && dim(alignment_rotation)[1] != n_embedding) msg <- c(msg, "`dim(alignment_rotation)[1]` does not match `n_embedding`")
  if(! is.null(alignment_rotation) && dim(alignment_rotation)[2] != n_embedding) msg <- c(msg, "`dim(alignment_rotation)[2]` does not match `n_embedding`")
  if(! is.null(alignment_rotation) && dim(alignment_rotation)[3] != ncol(alignment_design_matrix)) msg <- c(msg, "`dim(alignment_rotation)[3]` does not match `ncol(alignment_design_matrix)`")
  if(! is.null(alignment_stretching) && ! is.array(alignment_stretching) || length(dim(alignment_stretching)) != 3) msg <- c(msg, "`alignment_stretching` must be a three dimensional array")
  if(! is.null(alignment_stretching) && dim(alignment_stretching)[1] != n_embedding) msg <- c(msg, "`dim(alignment_stretching)[1]` does not match `n_embedding`")
  if(! is.null(alignment_stretching) && dim(alignment_stretching)[2] != n_embedding) msg <- c(msg, "`dim(alignment_stretching)[2]` does not match `n_embedding`")
  if(! is.null(alignment_stretching) && dim(alignment_stretching)[3] != ncol(alignment_design_matrix)) msg <- c(msg, "`dim(alignment_stretching)[3]` does not match `ncol(alignment_design_matrix)`")
  if(! is.null(alignment_design_matrix) && nrow(alignment_design_matrix) != n_obs) msg <- c(msg, "`nrow(alignment_design_matrix)` does not match number of observations")
  if(! is.null(alignment_design) &&  ! inherits(alignment_design, "formula")) msg <- c(msg, "`alignment_design` must inherit from formula or be NULL")
  if(! is.null(design) &&  ! inherits(design, "formula")) msg <- c(msg, "`design` must inherit from formula or be NULL")


  if(is.null(msg)){
    TRUE
  }else{
    msg
  }
})

# Subsetting

#' @export
setMethod("[", c("lemur_fit_obj", "ANY", "ANY"), function(x, i, j, ...) {
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }


  i_missing <- missing(i)
  j_missing <- missing(j)

  if (!missing(i)) {
    # Update metadata
    ii <- SingleCellExperiment:::.convert_subset_index(i, rownames(x))
    metadata(x)[["ambient_offset"]] <- metadata(x)[["ambient_offset"]][ii]
    metadata(x)[["ambient_coordsystem"]] <- metadata(x)[["ambient_coordsystem"]][ii,,drop=FALSE]
  }
  if(! missing(j)){
    jj <- SingleCellExperiment:::.convert_subset_index(j, rownames(x))
    metadata(x)[["alignment_design_matrix"]] <- metadata(x)[["alignment_design_matrix"]][jj,,drop=FALSE]
    am <- metadata(x)[["alignment_method"]]
    if(! isTRUE(am) && ! isFALSE(am)){
      metadata(x)[["alignment_method"]] <- metadata(x)[["alignment_method"]][jj]
    }
  }

  callNextMethod()
})



.methods_to_suggest <- c("n_ambient", "n_embedding",
                         "ambient_coordsystem", "ambient_offset",
                         "design", "base_point",
                         "coefficients", "embedding", "design_matrix", "linear_coefficients",
                         "alignment_method", "alignment_rotation", "alignment_stretching", "alignment_design", "alignment_design_matrix",
                         "colData", "rowData")

#' Get different features and elements of the 'lemur_fit_obj' object
#'
#' The following elements can all be accessed using the
#' fluent dollar notation (ie. \code{fit$n_ambient}).
#'
#' @usage
#'
#' fit$n_ambient
#' fit$n_embedding
#' fit$ambient_coordsystem
#' fit$ambient_offset
#' fit$design
#' fit$base_point
#' fit$coefficients
#' fit$embedding
#' fit$design_matrix
#' fit$linear_coefficients
#' fit$alignment_method
#' fit$alignment_rotation
#' fit$alignment_stretching
#' fit$alignment_design
#' fit$alignment_design_matrix
#' fit$colData
#' fit$rowData
#'
#' @param object the 'lemur_fit_obj' object
#'
#' @return The return value depends on the value
#' \describe{
#'  \item{n_ambient}{the number of ambient dimensions. `Inf` indicates that no dimensionality reduction was performed.}
#'  \item{n_embedding}{the number of embedding dimensions.}
#'  \item{ambient_coordsystem}{a matrix from the ambient PCA.}
#'  \item{ambient_offset}{a vector with the gene means.}
#'  \item{design}{the specification of the design in [`lemur`]. Usually this is a [`stats::formula`].}
#'  \item{base_point}{a matrix (`nrow(fit) * fit$n_embedding`) with the base point for the Grassmann exponential map.}
#'  \item{coefficients}{a three-dimensional tensor (`nrow(fit) * fit$n_embedding * ncol(fit$design_matrix)`) with the coefficients for
#'    the exponential map.}
#'  \item{embedding}{a matrix (`fit$n_embedding * ncol(fit)`) with the low dimensional position for each cell.}
#'  \item{design_matrix}{a matrix with covariates for each cell (`ncol(fit) * ncol(fit$design_matrix)`).}
#'  \item{linear_coefficients}{a matrix (`nrow(fit) * ncol(fit$design_matrix)`) with the coefficients for the linear regression.}
#'  \item{alignment_method}{a boolean. *Might be deleted or changed in a future version*.}
#'  \item{alignment_rotation}{a 3D tensor with the coefficients for the alignment rotation (`fit$n_embedding * fit$n_embedding * ncol(fit$design_matrix)`)}
#'  \item{alignment_stretching}{a 3D tensor with the coefficients for the alignment stretching (`fit$n_embedding * fit$n_embedding * ncol(fit$design_matrix)`)}
#'  \item{alignment_design}{an alternative design specification for the alignment. This is typically a [`stats::formula`].}
#'  \item{alignment_design_matrix}{an alternative design matrix specification for the alignment.}
#'  \item{colData}{the column annotation `DataFrame`.}
#'  \item{rowData}{the row annotation `DataFrame`.}
#' }
#'
#' @name accessor_methods
NULL


# The generic of design is in BiocGenerics

#' @rdname accessor_methods
#' @export
setMethod("design", signature = "lemur_fit_obj", function(object){
  metadata(object)[["design"]]
})



#' @rdname cash-lemur_fit_obj-method
#' @export
.DollarNames.lemur_fit_obj <- function(x, pattern = ""){
  grep(pattern, .methods_to_suggest, value = TRUE)
}

#' Fluent use of accessor methods
#'
#'
#' @seealso \link{accessor_methods} for more documentation on the
#'   accessor functions.
#' @aliases dollar_methods
setMethod("$", "lemur_fit_obj",
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
    n_ambient =               metadata(x)[["n_ambient"]],
    n_embedding =             metadata(x)[["n_embedding"]],
    design =                  design(x),
    base_point =              metadata(x)[["base_point"]],
    coefficients =            metadata(x)[["coefficients"]],
    embedding =               t(reducedDim(x, "embedding")),
    design_matrix =           sampleFactors(reducedDim(x, "linearFit")),
    linear_coefficients =     featureLoadings(reducedDim(x, "linearFit")),
    alignment_design =        metadata(x)[["alignment_design"]],
    alignment_design_matrix = metadata(x)[["alignment_design_matrix"]],
    alignment_stretching =    metadata(x)[["alignment_stretching"]],
    alignment_rotation =      metadata(x)[["alignment_rotation"]],
    alignment_method =        metadata(x)[["alignment_method"]],
    ambient_coordsystem =     metadata(x)[["ambient_coordsystem"]],
    ambient_offset =          metadata(x)[["ambient_offset"]],
    colData =                 colData(x),
    rowData =                 rowData(x),
    stop("Invalid `name` value.")
  )
})

#' @rdname cash-lemur_fit_obj-method
setReplaceMethod("$", "lemur_fit_obj",
                 function(x, name, value){
 if(! name %in% .methods_to_suggest){
   stop("Illegal name after '$' sign: ", name)
 }
 switch(name,
        colData = `colData<-`(x, value),
        rowData = `rowData<-`(x, value),
        stop("It is illegal to modify the content of lemur_fit_obj object")
 )
})



















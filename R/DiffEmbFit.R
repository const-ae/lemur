
#' @export
.DiffEmbFit <- setClass("DiffEmbFit", contains = "SingleCellExperiment")

DiffEmbFit <- function(data_mat, col_data, row_data,
                       n_ambient, n_embedding,
                       ambient_coordsystem, ambient_offset,
                       design, design_matrix, linear_coefficients,
                       diffemb_basepoint, diffemb_coefficients, diffemb_embedding,
                       alignment_method, alignment_rotation, alignment_stretching,
                       alignment_design_matrix, bootstrap_samples = NULL, knn_graph = NULL){


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
                                                 diffemb_embedding = t(diffemb_embedding)),
                              metadata = list(n_ambient = n_ambient, n_embedding = n_embedding,
                                              ambient_coordsystem = ambient_coordsystem, ambient_offset = ambient_offset,
                                              design = design,
                                              diffemb_basepoint = diffemb_basepoint, diffemb_coefficients = diffemb_coefficients,
                                              alignment_method = alignment_method, alignment_rotation = alignment_rotation, alignment_stretching = alignment_stretching,
                                              alignment_design_matrix = alignment_design_matrix,
                                              bootstrap_samples = bootstrap_samples, knn_graph = knn_graph))
  .DiffEmbFit(sce)
}


S4Vectors::setValidity2("DiffEmbFit", function(obj){
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
  diffemb_basepoint <- obj$diffemb_basepoint
  if(is.null(diffemb_basepoint)) msg <- c(msg, "'diffemb_basepoint' must not be NULL")
  diffemb_coefficients <- obj$diffemb_coefficients
  if(is.null(diffemb_coefficients)) msg <- c(msg, "'diffemb_coefficients' must not be NULL")
  alignment_method <- obj$alignment_method
  alignment_rotation <- obj$alignment_rotation
  if(is.null(alignment_rotation)) msg <- c(msg, "'alignment_rotation' must not be NULL")
  alignment_stretching <- obj$alignment_stretching
  if(is.null(alignment_stretching)) msg <- c(msg, "'alignment_stretching' must not be NULL")
  alignment_design_matrix <- obj$alignment_design_matrix
  if(is.null(alignment_design_matrix)) msg <- c(msg, "'alignment_stretching' must not be NULL")
  diffemb_embedding <- obj$diffemb_embedding
  if(is.null(diffemb_embedding)) msg <- c(msg, "'diffemb_embedding' must not be NULL")
  design <- obj$design
  design_matrix <- sampleFactors(reducedDim(obj, "linearFit"))
  if(is.null(design_matrix)) msg <- c(msg, "'design_matrix' must not be NULL")
  linear_coefficients <- featureLoadings(reducedDim(obj, "linearFit"))
  if(is.null(linear_coefficients)) msg <- c(msg, "'linear_coefficients' must not be NULL")
  bootstrap_samples <- obj$bootstrap_samples
  knn_graph <- obj$knn_graph


  if(! is.null(ambient_coordsystem) && nrow(ambient_coordsystem) != n_features) msg <- c(msg, "`nrow(ambient_coordsystem)` does not match the number of features")
  if(! is.null(ambient_coordsystem) && ncol(ambient_coordsystem) != eff_n_ambient) msg <- c(msg, "`ncol(ambient_coordsystem)` does not match `n_ambient`")
  if(! is.null(ambient_offset) && length(ambient_offset) != n_features) msg <- c(msg, "`length(ambient_offset)` does not match number of features")
  if(! is.null(design_matrix) && nrow(design_matrix) != n_obs) msg <- c(msg, "`nrow(design_matrix)` does not match number of observations")
  if(! is.null(linear_coefficients) && nrow(linear_coefficients) != eff_n_ambient) msg <- c(msg, "`nrow(linear_coefficients)` does not match `n_ambient`")
  if(! is.null(linear_coefficients) && ncol(linear_coefficients) != ncol(design_matrix)) msg <- c(msg, "`ncol(linear_coefficients)` does not match `ncol(design_matrix)`")
  if(! is.null(diffemb_basepoint) && nrow(diffemb_basepoint) != eff_n_ambient) msg <- c(msg, "`nrow(diffemb_basepoint)` does not match `n_ambient`")
  if(! is.null(diffemb_basepoint) && ncol(diffemb_basepoint) != n_embedding) msg <- c(msg, "`ncol(diffemb_basepoint)` does not match `n_embedding`")
  if(! is.null(diffemb_coefficients) && ! is.array(diffemb_coefficients) || length(dim(diffemb_coefficients)) != 3) msg <- c(msg, "`diffemb_coefficients` must be a three dimensional array")
  if(! is.null(diffemb_coefficients) && dim(diffemb_coefficients)[1] != eff_n_ambient) msg <- c(msg, "`dim(diffemb_coefficients)[1]` does not match `n_ambient`")
  if(! is.null(diffemb_coefficients) && dim(diffemb_coefficients)[2] != n_embedding) msg <- c(msg, "`dim(diffemb_coefficients)[2]` does not match `n_embedding`")
  if(! is.null(diffemb_coefficients) && dim(diffemb_coefficients)[3] != ncol(design_matrix)) msg <- c(msg, "`dim(diffemb_coefficients)[3]` does not match `ncol(design_matrix)`")
  if(! is.null(diffemb_embedding) && nrow(diffemb_embedding) != n_embedding) msg <- c(msg, "`nrow(diffemb_embedding)` does not match `n_embedding`")
  if(! is.null(diffemb_embedding) && ncol(diffemb_embedding) != n_obs) msg <- c(msg, "`ncol(diffemb_embedding)` does not match number of observations")
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
  if(! is.null(design) &&  ! inherits(design, "formula")) msg <- c(msg, "`design` must inherit from formula or be NULL")
  if(! is.null(bootstrap_samples) && any(vapply(bootstrap_samples, \(samp) ! is(samp, "DiffEmbFit"), FUN.VALUE = logical(1L)))) msg <- c(msg, "all bootstrap samples must be valid 'DiffEmbFit' objects.")
  if(! is.null(knn_graph) && ! inherits(knn_graph, "igraph")) msg <- c(msg, "knn_graph must be an object of class 'igraph'.")
  if(! is.null(knn_graph) && igraph::vcount(knn_graph) != n_obs) msg <- c(msg, "knn_graph must have one vertex for each observation.")


  if(is.null(msg)){
    TRUE
  }else{
    msg
  }
})

# Subsetting

#' @export
setMethod("[", c("DiffEmbFit", "ANY", "ANY"), function(x, i, j, ...) {
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }


  i_missing <- missing(i)
  j_missing <- missing(j)
  if(! is.null(x$bootstrap_samples)){
    metadata(x)$bootstrap_samples <- lapply(x$bootstrap_samples, \(bs) {
      if(!i_missing && !j_missing) bs[i,j,...]
      else if(!i_missing && j_missing) bs[i,,...]
      else if(i_missing && !j_missing) bs[,j,...]
      else bs[]
    })
  }

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
    if(! is.null(metadata(x)[["knn_graph"]])) metadata(x)[["knn_graph"]] <- igraph::induced_subgraph(metadata(x)[["knn_graph"]], jj)
  }

  callNextMethod()
})



.methods_to_suggest <- c("n_ambient", "n_embedding",
                         "ambient_coordsystem", "ambient_offset",
                         "design", "diffemb_basepoint",
                         "diffemb_coefficients", "diffemb_embedding", "design_matrix", "linear_coefficients",
                         "alignment_method", "alignment_rotation", "alignment_stretching", "alignment_design_matrix", 
                         "bootstrap_samples", "knn_graph", "colData", "rowData")

#' Get different features and elements of the 'DiffEmbFit' object
#'
#' The functions listed below can all be accessed using the
#' fluent dollar notation (ie. \code{fit$n_ambient}) without
#' any additional parentheses.
#'
#'
#' @param object the 'DiffEmbFit' object
#'
#' @return See the documentation of the generics to find out what each method returns
#'
#' @name accessor_methods
NULL



setGeneric("n_ambient", function(object, ...) standardGeneric("n_ambient"))

#' @rdname accessor_methods
#' @export
setMethod("n_ambient", signature = "DiffEmbFit", function(object){
  metadata(object)[["n_ambient"]]
})

setGeneric("n_embedding", function(object, ...) standardGeneric("n_embedding"))

#' @rdname accessor_methods
#' @export
setMethod("n_embedding", signature = "DiffEmbFit", function(object){
  metadata(object)[["n_embedding"]]
})

setGeneric("ambient_coordsystem", function(object, ...) standardGeneric("ambient_coordsystem"))

#' @rdname accessor_methods
#' @export
setMethod("ambient_coordsystem", signature = "DiffEmbFit", function(object){
  metadata(object)[["ambient_coordsystem"]]
})

setGeneric("ambient_offset", function(object, ...) standardGeneric("ambient_offset"))

#' @rdname accessor_methods
#' @export
setMethod("ambient_offset", signature = "DiffEmbFit", function(object){
  metadata(object)[["ambient_offset"]]
})

# The generic of design is in BiocGenerics

#' @rdname accessor_methods
#' @export
setMethod("design", signature = "DiffEmbFit", function(object){
  metadata(object)[["design"]]
})

setGeneric("diffemb_basepoint", function(object, ...) standardGeneric("diffemb_basepoint"))

#' @rdname accessor_methods
#' @export
setMethod("diffemb_basepoint", signature = "DiffEmbFit", function(object){
  metadata(object)[["diffemb_basepoint"]]
})

setGeneric("diffemb_coefficients", function(object, ...) standardGeneric("diffemb_coefficients"))

#' @rdname accessor_methods
#' @export
setMethod("diffemb_coefficients", signature = "DiffEmbFit", function(object){
  metadata(object)[["diffemb_coefficients"]]
})

setGeneric("diffemb_embedding", function(object, ...) standardGeneric("diffemb_embedding"))

#' @rdname accessor_methods
#' @export
setMethod("diffemb_embedding", signature = "DiffEmbFit", function(object){
  t(reducedDim(object, "diffemb_embedding"))
})

setGeneric("alignment_method", function(object, ...) standardGeneric("alignment_method"))

#' @rdname accessor_methods
#' @export
setMethod("alignment_method", signature = "DiffEmbFit", function(object){
  metadata(object)[["alignment_method"]]
})



setGeneric("alignment_rotation", function(object, ...) standardGeneric("alignment_rotation"))

#' @rdname accessor_methods
#' @export
setMethod("alignment_rotation", signature = "DiffEmbFit", function(object){
  metadata(object)[["alignment_rotation"]]
})

setGeneric("alignment_stretching", function(object, ...) standardGeneric("alignment_stretching"))

#' @rdname accessor_methods
#' @export
setMethod("alignment_stretching", signature = "DiffEmbFit", function(object){
  metadata(object)[["alignment_stretching"]]
})

setGeneric("alignment_design_matrix", function(object, ...) standardGeneric("alignment_design_matrix"))

#' @rdname accessor_methods
#' @export
setMethod("alignment_design_matrix", signature = "DiffEmbFit", function(object){
  metadata(object)[["alignment_design_matrix"]]
})


setGeneric("design_matrix", function(object, ...) standardGeneric("design_matrix"))

#' @rdname accessor_methods
#' @export
setMethod("design_matrix", signature = "DiffEmbFit", function(object){
  sampleFactors(reducedDim(object, "linearFit"))
})

setGeneric("linear_coefficients", function(object, ...) standardGeneric("linear_coefficients"))

#' @rdname accessor_methods
#' @export
setMethod("linear_coefficients", signature = "DiffEmbFit", function(object){
  featureLoadings(reducedDim(object, "linearFit"))
})

setGeneric("bootstrap_samples", function(object, ...) standardGeneric("bootstrap_samples"))

#' @rdname accessor_methods
#' @export
setMethod("bootstrap_samples", signature = "DiffEmbFit", function(object){
  metadata(object)[["bootstrap_samples"]]
})

setGeneric("knn_graph", function(object, ...) standardGeneric("knn_graph"))

#' @rdname accessor_methods
#' @export
setMethod("knn_graph", signature = "DiffEmbFit", function(object){
  metadata(object)[["knn_graph"]]
})


#' @rdname cash-DiffEmbFit-method
#' @export
.DollarNames.DiffEmbFit <- function(x, pattern = ""){
  grep(pattern, .methods_to_suggest, value = TRUE)
}

#' Fluent use of accessor methods
#'
#'
#' @seealso \link{accessor_methods} for more documentation on the
#'   accessor functions.
#' @aliases dollar_methods
setMethod("$", "DiffEmbFit",
          function(x, name){
            old <- S4Vectors:::disableValidity()
            if (!isTRUE(old)) {
              S4Vectors:::disableValidity(TRUE)
              on.exit(S4Vectors:::disableValidity(old))
            }

            if(! name %in% .methods_to_suggest){
              stop("Illegal name after '$' sign: ", name)
            }
            getGeneric(name)(x)
          })

#' @rdname cash-DiffEmbFit-method
setReplaceMethod("$", "DiffEmbFit",
                 function(x, name, value){
                   # if(name %in% c("rowData", "feature_parameters")){
                   #   getGeneric(paste0(name, "<-"))(x, value = value)
                   # }else{
                   #   stop("It is illegal to modify the content of ", name)
                   # }
                   stop("It is illegal to modify the content of DiffEmbFit object")
                 })



















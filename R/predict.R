

#' Predict values from `lemur_fit` object
#'
#' @param object an `lemur_fit` object
#' @param newdata a data.frame which passed to [`model.matrix`] with
#'   `design` to make the `newdesign` matrix
#' @param newdesign a matrix with the covariates for which the output
#'   is predicted. If `NULL`, the `object$design_matrix` is used. If
#'   it is a vector it is repeated `ncol(embedding)` times to create
#'   a design matrix with the same entry for each cell.
#' @param embedding the low-dimensional cell position for which the
#'   output is predicted.
#' @param with_linear_model a boolean to indicate if the linear regression
#'   offset is included in the prediction.
#' @param with_embedding a boolean to indicate if the embedding contributes
#'   to the output.
#' @param with_alignment a boolean to indicate if the alignment effect
#'   is removed from the output.
#' @param ... additional parameters passed to `predict_impl`.
#'
#' @return a matrix with the same dimension `nrow(object) * nrow(newdesign)`.
#'
#'
#' @seealso [`residuals`][residuals,lemur_fit-method]
#'
#' @export
predict.lemur_fit <- function(object, newdata = NULL, newdesign = NULL,
                               embedding = object$embedding,
                               with_linear_model = TRUE,
                               with_embedding = TRUE,
                               with_alignment = TRUE,
                               ...){
  predict_impl(object, newdata = newdata, newdesign = newdesign, embedding = embedding,
               with_linear_model = with_linear_model,
               with_embedding = with_embedding, with_alignment = with_alignment, ...)

}

predict_impl <- function(object, newdata = NULL, newdesign = NULL,
                         embedding = object$embedding,
                         with_linear_model = TRUE,
                         with_embedding = TRUE,
                         with_alignment = TRUE,
                         n_embedding = object$n_embedding,
                         design_matrix = object$design_matrix, design = object$design,
                         linear_coefficients = object$linear_coefficients,
                         coefficients = object$coefficients,
                         base_point = object$base_point,
                         alignment_coefficients = object$alignment_coefficients,
                         alignment_design_matrix = object$alignment_design_matrix,
                         row_mask = metadata(object)$row_mask,
                         ...){
  if(is.null(newdesign) && is.null(newdata)){
    newdesign <- design_matrix
  }else if(! is.null(newdata)){
    if(is.null(design)) stop("'newdata' is provided, but 'object' does not contain a design formula.")
    newdesign <- model.matrix(design, newdata)
  }else if(! is.matrix(newdesign)){
    newdesign <- matrix(newdesign, nrow = ncol(embedding), ncol = length(newdesign), byrow = TRUE)
  }
  if(! is.matrix(alignment_design_matrix)){
    alignment_design_matrix <- matrix(alignment_design_matrix, nrow = ncol(embedding), ncol = length(alignment_design_matrix), byrow = TRUE)
  }

  if(all(dim(design_matrix) == dim(alignment_design_matrix)) && all(design_matrix == alignment_design_matrix)){
    # The design matrices were identical, presume that the newdesign should also be identical
    alignment_design_matrix <- newdesign
  }

  if(nrow(newdesign) != ncol(embedding)){
    stop("The number of rows in 'newdesign' must match the number of columns in 'embedding'")
  }
  if(nrow(newdesign) != nrow(alignment_design_matrix)){
    stop("The number of rows in 'newdesign' (", nrow(newdesign) ,") and 'alignment_design_matrix'(", nrow(alignment_design_matrix) ,")  must be the same")
  }
  approx <- if(with_linear_model){
    linear_coefficients %*% t(newdesign)
  }else{
    matrix(0, nrow = min(nrow(linear_coefficients)), ncol = nrow(newdesign))
  }

  if(with_embedding){
    mm_groups <- get_groups(newdesign, n_groups = 100)
    mm_al_groups <- get_groups(alignment_design_matrix, n_groups = 100)
    stopifnot(length(mm_groups) == length(mm_al_groups))
    mmg <- unique(cbind(mm_groups, mm_al_groups))
    for(idx in seq_len(nrow(mmg))){
      gr1 <- mmg[idx,1]
      gr2 <- mmg[idx,2]
      covar1 <- newdesign[which(mm_groups == gr1)[1],]
      diffemb <- grassmann_map(sum_tangent_vectors(coefficients, covar1), base_point)
      alignment <- if(with_alignment){
        covar2 <- alignment_design_matrix[which(mm_al_groups == gr2)[1],]
        reverse_linear_transformation(alignment_coefficients, covar2)
      }else{
        diag(nrow = n_embedding)
      }
      sel <- gr1 == mm_groups & gr2 == mm_al_groups
      approx[,sel] <- approx[,sel] + diffemb %*% alignment %*% embedding[,sel]
    }
  }

  approx <- approx[row_mask,,drop=FALSE]
  colnames(approx) <- rownames(newdesign)
  rownames(approx) <- rownames(object)
  approx
}


#' Predict values from `lemur_fit` object
#'
#' @inheritParams predict.lemur_fit
#' @param ... ignored.
#'
#' @export
setMethod("residuals", signature = "lemur_fit", function(object,
                                                          with_linear_model = TRUE,
                                                          with_embedding = TRUE, ...){
  residuals_impl(object, with_linear_model = with_linear_model, with_embedding = with_embedding)
})


residuals_impl <- function(object,
                           with_linear_model = TRUE,
                           with_embedding = TRUE){
  assay(object, "expr") - predict(object, with_linear_model = with_linear_model, with_embedding = with_embedding)
}



get_residuals_for_alt_fit <- function(fit, Y = assay(fit, "expr"), reduced_design_mat, with_linear_model = TRUE, with_embedding = TRUE){
  if(with_embedding){
    fit_alt <- lemur_impl(Y, design_matrix = reduced_design_mat,  n_embedding = fit$n_embedding,
                          base_point = fit$base_point,  verbose = FALSE)
    Y - predict_impl(object = fit_alt, embedding = fit_alt$embedding,
                 with_linear_model = TRUE, with_embedding = TRUE,
                 n_embedding = fit_alt$n_embedding,
                 design_matrix = fit_alt$design_matrix, design = fit_alt$design,
                 linear_coefficients = fit_alt$linear_coefficients, coefficients = fit_alt$coefficients,
                 base_point = fit_alt$base_point, alignment_design_matrix = fit_alt$alignment_design_matrix,
                 alignment_coefficients = fit_alt$alignment_coefficients,
                 row_mask = rep(TRUE, nrow(Y)))
  }else{
    fit_alt <- lemur_impl(Y, design_matrix = reduced_design_mat,
                          n_embedding = 0,
                          base_point = matrix(nrow = nrow(fit), ncol = 0),
                          verbose = FALSE)
    Y - predict_impl(object = fit_alt, embedding = fit_alt$embedding,
                 with_linear_model = TRUE, with_embedding = FALSE,
                 n_embedding = fit_alt$n_embedding,
                 design_matrix = fit_alt$design_matrix, design = fit_alt$design,
                 linear_coefficients = fit_alt$linear_coefficients, coefficients = fit_alt$coefficients,
                 base_point = fit_alt$base_point, alignment_design_matrix = fit$alignment_design_matrix,
                 alignment_coefficients = fit_alt$alignment_coefficients,
                 row_mask = rep(TRUE, nrow(Y)))
  }
}



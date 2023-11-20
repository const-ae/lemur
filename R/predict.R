

#' Predict values from `lemur_fit` object
#'
#' @param object an `lemur_fit` object
#' @param newdata a data.frame which passed to [`model.matrix`] with
#'   `design` to make the `newdesign` matrix
#' @param newdesign a matrix with the covariates for which the output
#'   is predicted. If `NULL`, the `object$design_matrix` is used. If
#'   it is a vector it is repeated `ncol(embedding)` times to create
#'   a design matrix with the same entry for each cell.
#' @param newcondition an unquoted expression with a call to `cond()` specifying
#'   the covariates of the prediction. See the `contrast` argument in [test_de]
#'   for more details. Note that combinations of multiple calls to `cond()` are
#'   not allowed (e.g., `cond(a = 1) - cond(a = 2)`). If specified, `newdata`
#'   and `newdesign` are ignored.
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
#' @returns A matrix with the same dimension `nrow(object) * nrow(newdesign)`.
#'
#'
#' @seealso [`residuals`][residuals,lemur_fit-method]
#'
#' @examples
#'
#' data(glioblastoma_example_data)
#' fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition,
#'              n_emb = 5, verbose = FALSE)
#'
#' pred <- predict(fit)
#'
#' pred_ctrl <- predict(fit, newdesign = c(1, 0, 0, 0, 0, 0))
#' pred_trt <-  predict(fit, newdesign = c(1, 0, 0, 0, 0, 1))
#' # This is the same as the test_de result
#' fit <- test_de(fit, cond(condition = "panobinostat") - cond(condition = "ctrl"))
#' all.equal(SummarizedExperiment::assay(fit, "DE"), pred_trt - pred_ctrl,
#'           check.attributes = FALSE)
#'
#' @export
predict.lemur_fit <- function(object, newdata = NULL, newdesign = NULL,
                              newcondition = NULL,
                               embedding = object$embedding,
                               with_linear_model = TRUE,
                               with_embedding = TRUE,
                               with_alignment = TRUE,
                               ...){
  predict_impl(object, newdata = newdata, newdesign = newdesign, newcondition = {{newcondition}},
               embedding = embedding, with_linear_model = with_linear_model,
               with_embedding = with_embedding, with_alignment = with_alignment, ...)

}

predict_impl <- function(object, newdata = NULL, newdesign = NULL,
                         newcondition = NULL,
                         embedding = object$embedding,
                         with_linear_model = TRUE,
                         with_embedding = TRUE,
                         with_alignment = TRUE,
                         n_embedding = object$n_embedding,
                         design_matrix = object$design_matrix,
                         design = object$design,
                         linear_coefficients = object$linear_coefficients,
                         coefficients = object$coefficients,
                         base_point = object$base_point,
                         alignment_coefficients = object$alignment_coefficients,
                         alignment_design = object$alignment_design,
                         alignment_design_matrix = object$alignment_design_matrix,
                         row_mask = metadata(object)$row_mask,
                         ...){
  if(! rlang::quo_is_null(rlang::quo({{newcondition}}))){
    if(! is.null(newdesign) || !is.null(newdata)) warning("If 'newcondition' is used, 'newdesign' and 'newdata' are ignored.")
    newdesign <- parse_contrast({{newcondition}}, design)
    alignment_design_matrix <- parse_contrast({{newcondition}}, alignment_design)
    if(inherits(newdesign, "contrast_relation")) stop("Contrast relations using + or -  are not allowed")
    newdesign <- matrix(newdesign, nrow = ncol(embedding), ncol = length(newdesign), byrow = TRUE)
    alignment_design_matrix <- matrix(alignment_design_matrix, nrow = ncol(embedding), ncol = length(alignment_design_matrix), byrow = TRUE)
  }else if(is.null(newdesign) && is.null(newdata)){
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
    linear_coefficients[row_mask,,drop=FALSE] %*% t(newdesign)
  }else{
    matrix(0, nrow = length(row_mask), ncol = nrow(newdesign))
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
      if(with_alignment){
        covar2 <- alignment_design_matrix[which(mm_al_groups == gr2)[1],]
        alignment <- reverse_linear_transformation(alignment_coefficients, covar2)
        offset <- c(matrix(alignment_coefficients[,1,], ncol = length(covar2)) %*% covar2)
      }else{
        alignment <- diag(nrow = n_embedding)
        offset <- 0
      }
      sel <- gr1 == mm_groups & gr2 == mm_al_groups
      approx[,sel] <- approx[,sel] + diffemb[row_mask,,drop=FALSE] %*% (alignment %*% (embedding[,sel] - offset))
    }
  }

  colnames(approx) <- rownames(newdesign)
  rownames(approx) <- rownames(object)
  approx
}


#' Predict values from `lemur_fit` object
#'
#' @inheritParams predict.lemur_fit
#' @param ... ignored.
#'
#' @returns A matrix with the same dimension `dim(object)`.
#'
#' @seealso [predict.lemur_fit]
#'
#' @examples
#' data(glioblastoma_example_data)
#' fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition,
#'              n_emb = 5, verbose = FALSE)
#'
#' resid <- residuals(fit)
#' dim(resid)
#'
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
  assay(object, object$use_assay) - predict(object, with_linear_model = with_linear_model, with_embedding = with_embedding)
}



get_residuals_for_alt_fit <- function(fit, Y = assay(fit, fit$use_assay), reduced_design_mat, with_linear_model = TRUE, with_embedding = TRUE){
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



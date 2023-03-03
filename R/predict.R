

#' Predict values from `lemur_fit_obj` object
#'
#' @export
predict.lemur_fit_obj <- function(object, newdata = NULL, newdesign = NULL,
                               embedding = object$embedding,
                               with_ambient_pca = TRUE,
                               with_linear_model = TRUE,
                               with_differential_embedding = TRUE,
                               with_alignment = TRUE,
                               ...){
  predict_impl(object, newdata = newdata, newdesign = newdesign, embedding = embedding,
               with_ambient_pca = with_ambient_pca, with_linear_model = with_linear_model,
               with_differential_embedding = with_differential_embedding, with_alignment = with_alignment, ...)

}

predict_impl <- function(object, newdata = NULL, newdesign = NULL,
                         embedding = object$embedding,
                         with_ambient_pca = TRUE,
                         with_linear_model = TRUE,
                         with_differential_embedding = TRUE,
                         with_alignment = TRUE,
                         n_ambient = object$n_ambient, n_embedding = object$n_embedding,
                         design_matrix = object$design_matrix, design = object$design,
                         ambient_coordsystem = object$ambient_coordsystem, ambient_offset = object$ambient_offset,
                         linear_coefficients = object$linear_coefficients,
                         coefficients = object$coefficients,
                         basepoint = object$basepoint,
                         alignment_rotation = object$alignment_rotation,
                         alignment_stretching = object$alignment_stretching,
                         alignment_design_matrix = object$alignment_design_matrix,
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
    matrix(0, nrow = min(n_ambient, nrow(linear_coefficients)), ncol = nrow(newdesign))
  }

  if(with_differential_embedding){
    mm_groups <- get_groups(newdesign, n_groups = 100)
    mm_al_groups <- get_groups(alignment_design_matrix, n_groups = 100)
    stopifnot(length(mm_groups) == length(mm_al_groups))
    mmg <- unique(cbind(mm_groups, mm_al_groups))
    for(idx in seq_len(nrow(mmg))){
      gr1 <- mmg[idx,1]
      gr2 <- mmg[idx,2]
      covar1 <- newdesign[which(mm_groups == gr1)[1],]
      diffemb <- grassmann_map(sum_tangent_vectors(coefficients, covar1), basepoint)
      alignment <- if(with_alignment){
        covar2 <- alignment_design_matrix[which(mm_al_groups == gr2)[1],]
        spd_map(sum_tangent_vectors(alignment_stretching, covar2), diag(nrow = n_embedding)) %*%
          rotation_map(sum_tangent_vectors(alignment_rotation, covar2), diag(nrow = n_embedding))
      }else{
        diag(nrow = n_embedding)
      }
      sel <- gr1 == mm_groups & gr2 == mm_al_groups
      approx[,sel] <- approx[,sel] + diffemb %*% alignment %*% embedding[,sel]
    }
  }

  if(with_ambient_pca){
    res <- as.matrix(ambient_coordsystem %*% approx + ambient_offset)
    colnames(res) <- rownames(newdesign)
    rownames(res) <- rownames(object)
    res
  }else{
    colnames(approx) <- rownames(newdesign)
    rownames(approx) <- paste0("latent_", seq_len(nrow(approx)))
    approx
  }
}


#' Predict values from `lemur_fit_obj` object
#'
#' @export
setMethod("residuals", signature = "lemur_fit_obj", function(object,
                                                          with_linear_model = TRUE,
                                                          with_differential_embedding = TRUE, ...){
  residuals_impl(object, with_linear_model = with_linear_model, with_differential_embedding = with_differential_embedding)
})


residuals_impl <- function(object,
                           with_linear_model = TRUE,
                           with_differential_embedding = TRUE){
  assay(object, "expr") - predict(object, with_linear_model = with_linear_model, with_differential_embedding = with_differential_embedding)
}



get_residuals_for_alt_fit <- function(fit, Y = assay(fit, "expr"), reduced_design_mat, with_linear_model = TRUE, with_differential_embedding = TRUE){
  if(with_differential_embedding){
    fit_alt <- lemur_impl(matrix(nrow = nrow(fit), ncol = 0), design_matrix = reduced_design_mat,
                                           n_ambient = fit$n_ambient, n_embedding = fit$n_embedding,
                                           alignment = fit$alignment_method, base_point = fit$basepoint,
                                           amb_pca = list(coordsystem = fit$ambient_coordsystem,
                                                          embedding = as.matrix(t(fit$ambient_coordsystem) %*% (Y - fit$ambient_offset)),
                                                          offset = fit$ambient_offset),
                                           verbose = FALSE)
    Y - predict_impl(object = NULL, embedding = fit_alt$embedding,
                 with_linear_model = TRUE, with_differential_embedding = TRUE,
                 n_ambient = fit_alt$n_ambient, n_embedding = fit_alt$n_embedding,
                 design_matrix = fit_alt$design_matrix, design = fit_alt$design,
                 ambient_coordsystem = fit_alt$ambient_coordsystem, ambient_offset = fit_alt$ambient_offset,
                 linear_coefficients = fit_alt$linear_coefficients, coefficients = fit_alt$coefficients,
                 basepoint = fit_alt$basepoint, alignment_design_matrix = fit_alt$alignment_design_matrix,
                 alignment_rotation = fit_alt$alignment_rotation, alignment_stretching = fit_alt$alignment_stretching)
  }else{
    fit_alt <- lemur_impl(matrix(nrow = nrow(fit), ncol = 0), design_matrix = reduced_design_mat,
                                           n_ambient = fit$n_ambient, n_embedding = 0,
                                           alignment = fit$alignment_method, base_point = matrix(nrow = fit$n_ambient, ncol = 0),
                                           amb_pca = list(coordsystem = fit$ambient_coordsystem,
                                                          embedding = as.matrix(t(fit$ambient_coordsystem) %*% (Y - fit$ambient_offset)),
                                                          offset = fit$ambient_offset),
                                           verbose = FALSE)
    Y - predict_impl(object = NULL, embedding = fit_alt$embedding,
                 with_linear_model = TRUE, with_differential_embedding = FALSE,
                 n_ambient = fit_alt$n_ambient, n_embedding = fit_alt$n_embedding,
                 design_matrix = fit_alt$design_matrix, design = fit_alt$design,
                 ambient_coordsystem = fit_alt$ambient_coordsystem, ambient_offset = fit_alt$ambient_offset,
                 linear_coefficients = fit_alt$linear_coefficients, coefficients = fit_alt$coefficients,
                 basepoint = fit_alt$basepoint, alignment_design_matrix = fit$alignment_design_matrix,
                 alignment_rotation = fit_alt$alignment_rotation, alignment_stretching = fit_alt$alignment_stretching)
  }
}



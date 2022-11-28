

#' Predict values from `DiffEmbFit` object
#'
#' @export
setMethod("predict", signature = "DiffEmbFit", function(object, newdata = NULL, newdesign = NULL,
                                                        diffemb_embedding = object$diffemb_embedding,
                                                        with_ambient_pca = TRUE,
                                                        with_linear_model = TRUE,
                                                        with_differential_embedding = TRUE,
                                                        with_alignment = TRUE,
                                                        ...){
  predict_impl(object, newdata = newdata, newdesign = newdesign, diffemb_embedding = diffemb_embedding,
               with_ambient_pca = with_ambient_pca, with_linear_model = with_linear_model,
               with_differential_embedding = with_differential_embedding, with_alignment = with_alignment)

})

predict_impl <- function(object, newdata = NULL, newdesign = NULL,
                         diffemb_embedding = object$diffemb_embedding,
                         with_ambient_pca = TRUE,
                         with_linear_model = TRUE,
                         with_differential_embedding = TRUE,
                         with_alignment = TRUE,
                         n_ambient = object$n_ambient, n_embedding = object$n_embedding,
                         design_matrix = object$design_matrix, design = object$design,
                         ambient_coordsystem = object$ambient_coordsystem, ambient_offset = object$ambient_offset,
                         linear_coefficients = object$linear_coefficients,
                         diffemb_coefficients = object$diffemb_coefficients,
                         diffemb_basepoint = object$diffemb_basepoint,  alignment_coefficients = object$alignment_coefficients,
                         ...){
  if(is.null(newdesign) && is.null(newdata)){
    newdesign <- design_matrix
  }else if(! is.null(newdata)){
    if(is.null(design)) stop("'newdata' is provided, but 'object' does not contain a design formula.")
    newdesign <- model.matrix(design, newdata)
  }else if(! is.matrix(newdesign)){
    newdesign <- matrix(newdesign, nrow = ncol(diffemb_embedding), ncol = length(newdesign), byrow = TRUE)
  }
  if(nrow(newdesign) != ncol(diffemb_embedding)){
    stop("The number of rows in 'newdesign' must match the number of columns in 'diffemb_embedding'")
  }
  approx <- if(with_linear_model){
    linear_coefficients %*% t(newdesign)
  }else{
    matrix(0, nrow =min(n_ambient, nrow(object)), ncol = nrow(newdesign))
  }

  if(with_differential_embedding){
    mm_groups <- get_groups(newdesign, n_groups = 100)
    for(gr in unique(mm_groups)){
      covar <- newdesign[which(mm_groups == gr)[1],]
      diffemb <- grassmann_map(sum_tangent_vectors(diffemb_coefficients, covar), diffemb_basepoint)
      alignment <- if(with_alignment){
        rotation_map(sum_tangent_vectors(alignment_coefficients, covar), diag(nrow = n_embedding))
      }else{
        diag(nrow = n_embedding)
      }
      approx[,gr == mm_groups] <- approx[,gr == mm_groups] + diffemb %*% alignment %*% diffemb_embedding[,gr == mm_groups]
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


#' Predict values from `DiffEmbFit` object
#'
#' @export
setMethod("residuals", signature = "DiffEmbFit", function(object,
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
    fit_alt <- differential_embedding_impl(matrix(nrow = nrow(fit), ncol = 0), design_matrix = reduced_design_mat,
                                           n_ambient = fit$n_ambient, n_embedding = fit$n_embedding,
                                           alignment = fit$alignment_method, base_point = fit$diffemb_basepoint,
                                           amb_pca = list(coordsystem = fit$ambient_coordsystem,
                                                          embedding = as.matrix(t(fit$ambient_coordsystem) %*% (Y - fit$ambient_offset)),
                                                          offset = fit$ambient_offset),
                                           verbose = FALSE)
    Y - predict_impl(object = NULL, diffemb_embedding = fit_alt$diffemb_embedding,
                 with_linear_model = TRUE, with_differential_embedding = TRUE,
                 n_ambient = fit_alt$n_ambient, n_embedding = fit_alt$n_embedding,
                 design_matrix = fit_alt$design_matrix, design = fit_alt$design,
                 ambient_coordsystem = fit_alt$ambient_coordsystem, ambient_offset = fit_alt$ambient_offset,
                 linear_coefficients = fit_alt$linear_coefficients, diffemb_coefficients = fit_alt$diffemb_coefficients,
                 diffemb_basepoint = fit_alt$diffemb_basepoint, alignment_coefficients = fit_alt$alignment_coefficients)
  }else{
    fit_alt <- differential_embedding_impl(matrix(nrow = nrow(fit), ncol = 0), design_matrix = reduced_design_mat,
                                           n_ambient = fit$n_ambient, n_embedding = 0,
                                           alignment = fit$alignment_method, base_point = matrix(nrow = fit$n_ambient, ncol = 0),
                                           amb_pca = list(coordsystem = fit$ambient_coordsystem,
                                                          embedding = as.matrix(t(fit$ambient_coordsystem) %*% (Y - fit$ambient_offset)),
                                                          offset = fit$ambient_offset),
                                           verbose = FALSE)
    Y - predict_impl(object = NULL, diffemb_embedding = fit_alt$diffemb_embedding,
                 with_linear_model = TRUE, with_differential_embedding = FALSE,
                 n_ambient = fit_alt$n_ambient, n_embedding = fit_alt$n_embedding,
                 design_matrix = fit_alt$design_matrix, design = fit_alt$design,
                 ambient_coordsystem = fit_alt$ambient_coordsystem, ambient_offset = fit_alt$ambient_offset,
                 linear_coefficients = fit_alt$linear_coefficients, diffemb_coefficients = fit_alt$diffemb_coefficients,
                 diffemb_basepoint = fit_alt$diffemb_basepoint, alignment_coefficients = fit_alt$alignment_coefficients)
  }
}



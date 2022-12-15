

#' Main function to fit the differential embedding object
#'
#' @param data a matrix with obsevations in the columns and features in the rows.
#'   Or a `SummarizedExperiment` / `SingleCellExperiment` object
#' @param design a formula referring to global objects or column in the `colData` of `data`
#'   and `col_data` argument
#' @param col_data an optional data frame with `ncol(data)` rows.
#' @param n_ambient the dimension of the ambient PCA. This should be large enough
#'   to capture all relevant variation in the data
#' @param n_embedding the dimension of the $k$-plane that is rotated through space.
#'   Needs to be smaller than `n_ambient`.
#' @param alignment optional specification how or if points should be aligned. This
#'   can also be done in a separate step by calling [`align_embeddings`].
#' @param base_point a string specifying how to find the base point for the geodesic
#'   regression. Alternatively, an orthogonal matrix with `n_ambient` \eqn{\times} `n_embedding`
#'   dimension.
#' @param use_assay if `data` is a `SummarizedExperiment` / `SingleCellExperiment` object,
#'   which assay should be used.
#' @param ... additional parameters that are passed on to the internal function `differential_embedding_impl`.
#'   For example, if you have pre-calculated the ambient PCA, you can provide it as a list with
#'   the `coordsystem`, `embedding` and `offset` fields, to speed-up the fit.
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#'
#' @export
differential_embedding <- function(data, design = ~ 1, col_data = NULL,
                                   n_ambient = 30, n_embedding = 15,
                                   alignment = FALSE,
                                   base_point = c("global_embedding", "mean"),
                                   use_assay = "logcounts",
                                   ...,
                                   verbose = TRUE){

  data_mat <- handle_data_parameter(data, on_disk = FALSE, assay = use_assay)

  col_data <- glmGamPoi:::get_col_data(data, col_data)
  des <- handle_design_parameter(design, data, col_data)


  res <- differential_embedding_impl(data_mat, des$model_matrix,
                                     n_ambient = n_ambient, n_embedding = n_embedding,
                                     alignment = alignment, base_point = base_point,
                                     verbose = verbose, ...)

  DiffEmbFit(data_mat, col_data = col_data, row_data = if(is(data, "SummarizedExperiment")) rowData(data) else NULL,
             n_ambient = res$n_ambient, n_embedding = res$n_embedding,
             ambient_coordsystem = res$ambient_coordsystem, ambient_offset = res$ambient_offset,
             design = des$design_formula, design_matrix = res$design_matrix,
             linear_coefficients = res$linear_coefficients,
             diffemb_basepoint = res$diffemb_basepoint,
             diffemb_coefficients = res$diffemb_coefficients,
             diffemb_embedding = res$diffemb_embedding,
             alignment_method = res$alignment_method,
             alignment_coefficients = res$alignment_coefficients)
}


differential_embedding_impl <- function(Y, design_matrix,
                                        n_ambient = 30, n_embedding = 15,
                                        alignment = FALSE,
                                        base_point = c("global_embedding", "mean"),
                                        amb_pca = NULL,
                                        linear_coefficients = NULL,
                                        diffemb_coefficients = NULL,
                                        diffemb_embedding = NULL,
                                        alignment_coefficients = NULL,
                                        n_iter = 10, tol = 1e-8,
                                        reshuffling_fraction = 0,
                                        verbose = TRUE){
  alignment_fixed_but_embedding_fitted <- ! is.null(alignment_coefficients) && is.null(diffemb_embedding)

  # Set reduced dimensions
  stopifnot(n_ambient >= 0 && n_embedding >= 0)
  n_embedding <- min(n_embedding, nrow(Y), n_ambient)
  linear_coef_fixed <-  ! is.null(linear_coefficients)
  diffemb_coef_fixed <- ! is.null(diffemb_coefficients)
  diffemb_embedding_fixed <- ! is.null(diffemb_embedding)
  alignment_coef_fixed <- ! is.null(alignment_coefficients)

  # Reduce to ambient space
  if(is.null(amb_pca)){
    if(n_ambient > nrow(Y)){
      if(verbose) message("Skip ambient PCA step")
      n_ambient <- Inf
      offset <- MatrixGenerics::rowMeans2(Y)
      amb_pca <- list(
        coordsystem = Matrix::Diagonal(n = nrow(Y)),
        embedding = Y - offset,
        offset = offset
      )
    }else{
      if(verbose) message("Fit ambient PCA")
      amb_pca <- pca(Y, n_ambient)
    }
    n_ambient_eff <- min(nrow(Y), n_ambient)
  }else{
    # Check that amb_pca is correct
    stopifnot(all(names(amb_pca) %in% c("coordsystem", "embedding", "offset")))
    n_ambient_eff <- min(nrow(Y), n_ambient)
    stopifnot(ncol(amb_pca$coordsystem) == n_ambient_eff)
    stopifnot(nrow(amb_pca$embedding) == n_ambient_eff)
    stopifnot(length(amb_pca$offset) == nrow(Y))
    if(ncol(Y) > 0){
      rand_sel <- sample(seq_len(ncol(Y)), min(ncol(Y), 100))
      pred_emb <- t(amb_pca$coordsystem) %*% (Y[,rand_sel,drop=FALSE] - amb_pca$offset)
      if(! all(abs(pred_emb - amb_pca$embedding[,rand_sel,drop=FALSE]) < 1e-8)){
        stop("The provided ambient PCA ('amb_pca') does not match the observed data. Does 'use_assay' match the assay that was used to calculate the PCA?")
      }
    }
  }

  if(reshuffling_fraction != 0){
    stopifnot(reshuffling_fraction > 0 && reshuffling_fraction <= 1)
    warning("Reshuffling elements from the design matrix is a beta feature, that was ",
            "designed to regularize the differential inference for unmatched populations.\n",
            "On the down-side a large 'reshuffle_fraction' can adversely affect the accuracy of the inference.",
            "Please use with care.", call. = FALSE)
    sel <- sample.int(nrow(design_matrix), size = round(nrow(design_matrix) * reshuffling_fraction), replace = FALSE)
    shuf_sel <- sample(sel)
    design_matrix[sel,] <- design_matrix[shuf_sel,]
  }

  # Initialize values
  if(linear_coef_fixed){
    Y_clean <- amb_pca$embedding - linear_coefficients %*% t(design_matrix)
  }else{
    if(verbose) message("Regress out global effects")
    linear_fit <- lm.fit(design_matrix, t(amb_pca$embedding))
    linear_coefficients <- t(linear_fit$coefficients)
    Y_clean <- t(linear_fit$residuals)
  }
  if(!is.matrix(base_point)){
    if(verbose) message("Find base point for differential embedding")
    base_point <- find_base_point(Y_clean, base_point, n_embedding = n_embedding)
  }

  last_round_error <- sum(amb_pca$embedding^2)
  if(verbose) message("Fit differential embedding model")
  if(verbose) message("-Iteration: ", 0, "\terror: ", sprintf("%.3g", last_round_error))
  for(iter in seq_len(n_iter)){
    if(! diffemb_coef_fixed){
      if(verbose) message("---Fit Grassmann linear model")
      diffemb_coefficients <- grassmann_lm(Y_clean, design = design_matrix, base_point = base_point)
    }
    if(! diffemb_embedding_fixed){
      diffemb_embedding <- project_data_on_diffemb(Y_clean, design = design_matrix,
                                                   diffemb_coefficients = diffemb_coefficients, base_point = base_point)
    }
    if(! linear_coef_fixed){
      if(verbose) message("---Update linear regression")
      Y_clean <- amb_pca$embedding - project_diffemb_into_data_space(diffemb_embedding, design = design_matrix,
                                                                     diffemb_coefficients = diffemb_coefficients, base_point = base_point)
      if(any(is.na(Y_clean))){
        linear_fit <- list(coefficients = matrix(NA, nrow = ncol(design_matrix), ncol = nrow(Y_clean)),
                           residuals = NA)
      }else{
        linear_fit <- lm.fit(design_matrix, t(Y_clean))
      }
      linear_coefficients <- t(linear_fit$coefficients)
    }
    Y_clean <- amb_pca$embedding - linear_coefficients %*% t(design_matrix)
    error <- sum(linear_fit$residuals^2)
    if(verbose) message("-Iteration: ", iter, "\terror: ", sprintf("%.3g", error))
    if(is.na(error) || abs(last_round_error - error) / (error + 0.5) < tol){
      if(verbose) message("Converged")
      break
    }else{
      last_round_error <- error
    }
  }

  if(alignment_fixed_but_embedding_fitted){
    # Rotate the diffemb_embedding if it wasn't provided
    stop("Fixing 'alignment_coefficients' without fixing 'diffemb_embedding' is not implemented")
  }else if(! alignment_coef_fixed && ! isFALSE(alignment)){
    if(verbose) message("Align points")
    stop("Cannot handle 'alignment='", alignment)
  }else{
    alignment_coefficients <- array(0, c(n_embedding, n_embedding, ncol(design_matrix)))
  }



  list(n_ambient = n_ambient, n_embedding = n_embedding,
       ambient_coordsystem = amb_pca$coordsystem, ambient_offset = amb_pca$offset,
       design_matrix = design_matrix, data = Y,
       linear_coefficients = linear_coefficients,
       diffemb_basepoint = base_point,
       diffemb_coefficients = diffemb_coefficients,
       diffemb_embedding = diffemb_embedding,
       alignment_method = alignment,
       alignment_coefficients = alignment_coefficients)
}


find_base_point <- function(Y_clean, base_point, n_embedding){
  n_ambient <- nrow(Y_clean)
  if(is.matrix(base_point)){
    stopifnot(nrow(base_point) == n_ambient)
    stopifnot(ncol(base_point) == n_embedding)

    # Check if it is orthogonal
    orth <- t(base_point) %*% base_point
    if(sum((orth - diag(nrow = n_embedding))^2) > 1e-8){
      stop("The provided 'base_point'  is not orthogonal")
    }
    base_point
  }else{
    base_point_meth <- match.arg(base_point, c("global_embedding", "mean"))
    if(base_point_meth == "global_embedding"){
      pca(Y_clean, n = n_embedding)$coordsystem
    }else if(base_point_meth == "mean"){
      stop("'base_point = \"mean\"' is not implemented. Please use 'global_embedding'.")
    }
  }
}


project_diffemb_into_data_space <- function(diffemb_embedding, design, diffemb_coefficients, base_point){
  n_amb <- nrow(base_point)
  res <- matrix(NA, nrow = n_amb, ncol = ncol(diffemb_embedding))
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  for(gr in unique(mm_groups)){
    covars <- design[which(mm_groups == gr)[1], ]
    res[,mm_groups == gr] <- grassmann_map(sum_tangent_vectors(diffemb_coefficients, covars), base_point) %*% diffemb_embedding[,mm_groups == gr,drop=FALSE]
  }
  res
}

project_data_on_diffemb <- function(Y_clean, design, diffemb_coefficients, base_point){
  n_emb <- ncol(base_point)
  res <- matrix(NA, nrow = n_emb, ncol = ncol(Y_clean))
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  for(gr in unique(mm_groups)){
    covars <- design[which(mm_groups == gr)[1], ]
    res[,mm_groups == gr] <- t(grassmann_map(sum_tangent_vectors(diffemb_coefficients, covars), base_point)) %*% Y_clean[,mm_groups == gr,drop=FALSE]
  }
  res
}

sum_tangent_vectors <- function(tangent_block, covariates){
  stopifnot(length(covariates) == dim(tangent_block)[3])
  res <- matrix(0, nrow = dim(tangent_block)[1], ncol = dim(tangent_block)[2])
  for(idx in seq_len(length(covariates))){
    res <- res + tangent_block[,,idx] * covariates[idx]
  }
  res
}

correct_fit <- function(fit, correction){
  reducedDim(fit, "diffemb_embedding") <- t(correction$diffemb_embedding)
  metadata(fit)[["alignment_coefficients"]] <-  metadata(fit)[["alignment_coefficients"]] + correction$rotation_coefficients
  metadata(fit)[["bootstrap_samples"]] <- lapply(metadata(fit)[["bootstrap_samples"]], \(samp){
    reducedDim(samp, "diffemb_embedding") <- t(apply_rotation(
      apply_stretching(samp$diffemb_embedding, correction$stretch_coefficients, design_matrix, correction$stretch_base_point),
      correction$rotation_coefficients, design_matrix, correction$rotation_base_point))
    metadata(samp)[["alignment_coefficients"]] <- metadata(samp)[["alignment_coefficients"]] + correction$rotation_coefficients
    samp
  })
  fit
}






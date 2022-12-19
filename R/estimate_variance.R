
#' Bootstrap the fit to get an estimate of the parameter variances
#'
#' @param fit the result produced by [`differential_embedding`]
#' @param n_bootstrap_samples how many bootstrap samples to produce. Note that
#'   the function has a high startup cost, but each additional bootstrap is fast.
#' @param refit_ambient_pca Fitting the ambient PCA is often the slowest step when
#'   fitting the differential embedding model. Thus, by default, it is set to `FALSE`
#'   for the bootstrap samples.
#' @param refit_linear_model,refit_differential_embedding,refit_alignment additional
#'   flags to set which part of the model are copied from the original fit and which are
#'   bootstrapped. By default, they are all `TRUE`
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#'
#' @export
estimate_variance <- function(fit, n_bootstrap_samples = 100,
                              replicates_by = NULL,
                              refit_ambient_pca = FALSE,
                              refit_linear_model = TRUE,
                              refit_differential_embedding = TRUE,
                              refit_alignment = TRUE,
                              verbose = TRUE){
  if(! refit_ambient_pca){
    # This is the slowest step of the bootstrapping
    # (I don't think I can do much about it though)
    original_embedding <- as.matrix(t(fit$ambient_coordsystem) %*% (assay(fit, "expr") - fit$ambient_offset))
  }

  replicates <- rlang::eval_tidy(rlang::enquo(replicates_by), data = as.list(fit$colData))
  if(! is.null(replicates)){
    if(length(unique(replicates)) <= 5) warning("The number of replicates is small. The variance estimates might be unreliable.")
  }

  bootstraps <- lapply(seq_len(n_bootstrap_samples), function(idx){
    if(verbose) message("Start bootstrap iteration ", idx)
    # resampling <- sample(seq_len(ncol(fit)), size = ncol(fit), replace = TRUE)
    resampling <- resample(ncol(fit), cluster = replicates)

    design_matrix <- fit$design_matrix[resampling,,drop=FALSE]
    alignment_design_matrix <- fit$alignment_design_matrix[resampling,,drop=FALSE]
    base_point <- fit$diffemb_basepoint

    if(refit_ambient_pca){
      data_mat <- assay(fit, "expr")[,resampling,drop=FALSE]
      amb_pca <- NULL
    }else{
      # Project data on coordinate system
      embedding <- original_embedding[,resampling,drop=FALSE]
      # If amb_pca is provided, 'differential_embedding_impl' doesn't need 'data_mat'
      # except to check that 'n_ambient' is correct
      data_mat <- matrix(nrow = nrow(fit), ncol = 0)
      amb_pca <- list(coordsystem = fit$ambient_coordsystem, embedding = embedding, offset = fit$ambient_offset)
    }

    if(refit_linear_model){
      linear_coefficients <- NULL
    }else{
      linear_coefficients <- fit$linear_coefficients
    }

    if(refit_differential_embedding){
      diffemb_coefficients <- NULL
      diffemb_embedding <- NULL
    }else{
      diffemb_coefficients <- fit$diffemb_coefficients
      diffemb_embedding <- fit$diffemb_embedding[,resampling,drop=FALSE]
    }

    alignment_method <- fit$alignment_method
    if(length(alignment_method) == ncol(fit)){
      alignment_method <- alignment_method[resampling]
    }
    if(refit_alignment){
      alignment_rotation <- NULL
      alignment_stretching <- NULL
    }else{
      alignment_rotation <- fit$alignment_rotation
      alignment_stretching <- fit$alignment_stretching
    }
    res <- differential_embedding_impl(data_mat, design_matrix = design_matrix,
                                n_ambient = fit$n_ambient, n_embedding = fit$n_embedding,
                                alignment = alignment_method, base_point = base_point,
                                amb_pca = amb_pca,
                                linear_coefficients = linear_coefficients,
                                diffemb_coefficients = diffemb_coefficients,
                                diffemb_embedding = diffemb_embedding,
                                alignment_rotation = alignment_rotation,
                                alignment_stretching = alignment_stretching,
                                alignment_design_matrix = alignment_design_matrix,
                                verbose = FALSE)

    # Use the results of the resampled fit, to make 'diffemb_embedding' correspond to the original data
    refitted_Y_clean <- if(refit_ambient_pca){
      as.matrix(t(res$ambient_coordsystem) %*% (assay(fit, "expr") - res$ambient_offset) - res$linear_coefficients %*% t(fit$design_matrix))
    }else{
      original_embedding - res$linear_coefficients %*% t(fit$design_matrix)
    }
    refitted_diffemb_embedding <- project_data_on_diffemb(refitted_Y_clean, design = fit$design_matrix,
                                                          diffemb_coefficients = res$diffemb_coefficients,
                                                          base_point = base_point)

    # Create DiffEmbFit object
    samp <- DiffEmbFit(NULL, col_data = as.data.frame(matrix(nrow = ncol(fit), ncol = 0)),
               row_data = as.data.frame(matrix(nrow = nrow(fit), ncol = 0)),
               n_ambient = res$n_ambient, n_embedding = res$n_embedding,
               ambient_coordsystem = res$ambient_coordsystem, ambient_offset = res$ambient_offset,
               design = fit$design, design_matrix = fit$design_matrix,
               linear_coefficients = res$linear_coefficients,
               diffemb_basepoint = res$diffemb_basepoint,
               diffemb_coefficients = res$diffemb_coefficients,
               diffemb_embedding = refitted_diffemb_embedding,
               alignment_method = alignment_method,
               alignment_rotation = res$alignment_rotation,
               alignment_stretching = res$alignment_stretching,
               alignment_design_matrix = fit$alignment_design_matrix)
    rownames(samp) <- rownames(fit)
    colnames(samp) <- colnames(fit)
    samp
  })
  metadata(fit)[["bootstrap_samples"]] <- bootstraps
  fit
}





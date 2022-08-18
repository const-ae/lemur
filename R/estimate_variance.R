

estimate_variance <- function(fit, n_bootstrap_samples = 100,
                              refit_ambient_pca = TRUE,
                              refit_linear_model = TRUE,
                              refit_differential_embedding = TRUE,
                              refit_alignment = TRUE){

  bootstraps <- lapply(seq_len(n_bootstrap_samples), function(idx){

    resampling <- sample(seq_len(ncol(fit)), size = ncol(fit), replace = TRUE)
    design_matrix <- fit$design_matrix[resampling,,drop=FALSE]
    data_mat <- assay(fit, "expr")[,resampling,drop=FALSE]
    base_point <- fit$diffemb_basepoint

    if(refit_ambient_pca){
      amb_pca <- NULL
    }else{
      # Project data on coordinate system
      embedding <- t(fit$ambient_coordsystem) %*% data_mat
      # If amb_pca is provided, 'differential_embedding_impl' doesn't need 'data_mat'
      # except to check that 'n_ambient' is correct
      data_mat <- matrix(nrow = nrow(data_mat), ncol = 0)
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
      alignment_coefficients <- NULL
    }else{
      alignment_coefficients <- fit$alignment_coefficients
    }

    res <- differential_embedding_impl(data_mat, design_matrix = design_matrix,
                                n_ambient = fit$n_ambient, n_embedding = fit$n_embedding,
                                alignment = alignment_method, base_point = base_point,
                                amb_pca = amb_pca,
                                linear_coefficients = linear_coefficients,
                                diffemb_coefficients = diffemb_coefficients,
                                diffemb_embedding = diffemb_embedding,
                                alignment_coefficients = alignment_coefficients, verbose = FALSE)

    # Use the results of the resampled fit, to make 'diffemb_embedding' corresponding to the original data
    refitted_Y_clean <- t(res$ambient_coordsystem) %*% (assay(fit, "expr") - res$ambient_offset) - res$linear_coefficients %*% t(fit$design_matrix)
    refitted_diffemb_embedding <-project_data_on_diffemb(refitted_Y_clean, design = fit$design_matrix,
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
               alignment_coefficients = res$alignment_coefficients)
    rownames(samp) <- rownames(fit)
    colnames(samp) <- colnames(fit)
    samp
  })
  metadata(fit)[["bootstrap_samples"]] <- bootstraps
  fit
}





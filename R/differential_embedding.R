


differential_embedding <- function(data, design = ~ 1, col_data = NULL,
                                   n_ambient = 30, n_embedding = 15,
                                   alignment = FALSE,
                                   base_point = c("global_embedding", "mean"),
                                   ...,
                                   verbose = TRUE){

  data_mat <- glmGamPoi:::handle_data_parameter(data, FALSE)

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
                                        verbose = TRUE){
  alignment_fixed_but_embedding_fitted <- ! is.null(alignment_coefficients) && is.null(diffemb_embedding)

  # Set reduced dimensions
  stopifnot(n_embedding > 0 && n_embedding > 0)
  n_ambient <- min(nrow(Y), n_ambient)
  n_embedding <- min(n_embedding, n_ambient)

  # Reduce to ambient space
  if(is.null(amb_pca)){
    if(verbose) message("Fit ambient PCA")
    amb_pca <- pca(Y, n_ambient)
  }

  # Regress out overall effects
  if(is.null(linear_coefficients)){
    if(verbose) message("Regress out global effects")
    linear_fit <- lm.fit(design_matrix, t(amb_pca$embedding))
    linear_coefficients <- t(linear_fit$coefficients)
    Y_clean <- t(linear_fit$residuals)
  }else{
    Y_clean <- amb_pca$embedding - linear_coefficients %*% t(design_matrix)
  }

  if(!is.matrix(base_point)){
    if(verbose) message("Find base point for differential embedding")
    base_point <- find_base_point(Y_clean, base_point, n_embedding = n_embedding)
  }
  if(is.null(diffemb_coefficients)){
    if(verbose) message("Fit Grassmann linear model")
    diffemb_coefficients <- grassmann_lm(Y_clean, design = design_matrix, base_point = base_point)
  }
  if(is.null(diffemb_embedding)){
    diffemb_embedding <- project_data_on_diffemb(Y_clean, design = design_matrix,
                                                 diffemb_coefficients = diffemb_coefficients, base_point = base_point)
  }

  if(alignment_fixed_but_embedding_fitted){
    # Rotate the diffemb_embedding if it wasn't provided
    stop("Fixing 'alignment_coefficients' without fixing 'diffemb_embedding' is not implemented")
  }else if(is.null(alignment_coefficients)){
    if(verbose) message("Align points")
    align_res <- align_points_impl(alignment, diffemb_embedding, design_matrix)
    diffemb_embedding <- align_res$diffemb_embedding
    alignment_coefficients <- align_res$alignment_coefficients
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


project_data_on_diffemb <- function(Y_clean, design, diffemb_coefficients, base_point){
  n_emb <- ncol(base_point)
  res <- matrix(NA, nrow = n_emb, ncol = ncol(Y_clean))
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  for(gr in unique(mm_groups)){
    covars <- design[which(mm_groups == gr)[1], ]
    res[,mm_groups == gr] <- t(grassmann_map(sum_tangent_vectors(diffemb_coefficients, covars), base_point)) %*% Y_clean[,mm_groups ==gr,drop=FALSE]
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


align_embeddings <- function(fit, alignment = TRUE, verbose = FALSE){
  if(isTRUE(alignment)){
    # Cluster each condition
    exp_group <- get_groups(design_matrix, 10)
    exp_group_levels <- unique(exp_group)
    Y_clean <- t(fit$ambient_coordsystem) %*% assay(fit, "expr") - fit$linear_coefficients %*% t(fit$design_matrix)

    clusters <- lapply(exp_group_levels, \(gr) kmeans(t(Y_clean[,exp_group == gr,drop=FALSE]), nstart = 20, centers = 10))
    stop("Not implemented")
  }


  align_res <- align_points_impl(alignment, fit$diffemb_embedding, fit$design_matrix, verbose = verbose)
  metadata(fit)[["alignment_method"]] <- alignment
  reducedDim(fit, "diffemb_embedding") <- t(align_res$diffemb_embedding)
  metadata(fit)[["alignment_coefficients"]] <- align_res$alignment_coefficients
  exp_group <- get_groups(fit$design_matrix, 10)
  exp_group_levels <- unique(exp_group)
  base_point <- diag(nrow = fit$n_embedding)
  metadata(fit)[["bootstrap_samples"]] <- lapply(metadata(fit)[["bootstrap_samples"]], \(samp){
    # Rotate all datapoints
    for(gr in exp_group_levels){
      dir <- sum_tangent_vectors(-align_res$alignment_coefficients, fit$design_matrix[which(exp_group == gr)[1],])
      reducedDim(fit, "diffemb_embedding")[gr == exp_group, ] <- t(rotation_map(dir, base_point) %*% samp$diffemb_embedding[,gr == exp_group])
    }
    metadata(samp)[["alignment_method"]] <- alignment
    metadata(samp)[["alignment_coefficients"]] <- align_res$alignment_coefficients
    samp
  })

  fit
}


align_points_impl <- function(alignment, diffemb_embedding, design_matrix,
                              n_iter = 5, verbose = FALSE){
  n_embedding <- nrow(diffemb_embedding)
  n_obs <- nrow(design_matrix)
  base_point <- diag(nrow = n_embedding)

   if(isFALSE(alignment)){
    coef <- array(0, dim = c(n_embedding, n_embedding, ncol(design_matrix)))
  }else{
    stopifnot(length(alignment) == n_obs)
    exp_group <- get_groups(design_matrix, 10)
    exp_group_levels <- unique(exp_group)
    n_exp_group_levels <- length(exp_group_levels)

    if(is.null(exp_group)) stop("Too many different experimental groups to do alignment.")
    alignment_levels <- unique(alignment)
    alignment_levels <- alignment_levels[! is.na(alignment_levels)]
    n_alignment_levels <- length(alignment_levels)
    alignment_red <- rep(alignment_levels, each = n_exp_group_levels)
    exp_group_red <- rep(exp_group_levels, times = n_alignment_levels)

    # Take the first row of each experimental design group level and repeated for each alignment level
    new_design <- mply_dbl(exp_group_levels, \(gr) design_matrix[which(exp_group == gr)[1],], ncol = ncol(design_matrix)) |>
      duplicate_rows(times = n_alignment_levels)
    # Take the mean of observations per alignment and experimental design group
    new_embedding <- matrix(NA, nrow = n_embedding, ncol = n_exp_group_levels * n_alignment_levels)
    idx <- 1
    for(al in alignment_levels){
      for(gr in exp_group_levels){
        new_embedding[,idx] <- rowMeans(diffemb_embedding[,!is.na(alignment) & al == alignment & gr == exp_group,drop=FALSE])
        idx <- idx + 1
      }
    }
    coef <- array(0, dim  = c(n_embedding, n_embedding, ncol(design_matrix)))
    new_embedding_cp <- new_embedding
    for(iter in seq_len(n_iter)){
      # Rotate centers
      for(gr in exp_group_levels){
        new_embedding[,gr == exp_group_red] <- rotation_map(sum_tangent_vectors(coef, design_matrix[which(exp_group == gr)[1],]), base_point) %*% new_embedding_cp[,gr == exp_group_red]
      }
      # Calculate the mean per alignment level and replicate for each experimental design group
      mean_embeddings <- t(mply_dbl(alignment_levels, \(al) rowMeans(new_embedding[,al == alignment_red,drop=FALSE], na.rm = TRUE), ncol = n_embedding)) |>
        duplicate_cols(each = n_exp_group_levels)
      if(verbose) message("Iter: ", iter, "\tError: ", sum((new_embedding - mean_embeddings)^2, na.rm = TRUE))
      # Align data to centers
      coef <- coef + rotation_lm(data = mean_embeddings,  design = new_design, obs_embedding = new_embedding,
                                 base_point = base_point, tangent_regression = TRUE)
    }


    # Rotate all datapoints
    for(gr in exp_group_levels){
      diffemb_embedding[,gr == exp_group] <- rotation_map(sum_tangent_vectors(coef, design_matrix[which(exp_group == gr)[1],]), base_point) %*% diffemb_embedding[,gr == exp_group]
    }

  }

  list(alignment_coefficients = -coef, diffemb_embedding = diffemb_embedding)
}


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
  if(is.null(newdesign) && is.null(newdata)){
    newdesign <- object$design_matrix
  }else if(! is.null(newdata)){
    if(is.null(fit$design)) stop("'newdata' is provided, but 'object' does not contain a design formula.")
    newdesign <- model.matrix(fit$design, newdata)
  }else if(! is.matrix(newdesign)){
    newdesign <- matrix(newdesign, nrow = ncol(object), ncol = length(newdesign), byrow = TRUE)
  }
  if(nrow(newdesign) != ncol(diffemb_embedding)){
    stop("The number of rows in 'newdesign' must match the number of columns in 'diffemb_embedding'")
  }
  approx <- if(with_linear_model){
    object$linear_coefficients %*% t(newdesign)
  }else{
    matrix(0, nrow = object$n_embedding, ncol = nrow(newdesign))
  }

  mm_groups <- get_groups(newdesign, n_groups = 100)
  for(gr in unique(mm_groups)){
    covar <- newdesign[which(mm_groups == gr)[1],]
    diffemb <- if(with_differential_embedding){
      grassmann_map(sum_tangent_vectors(object$diffemb_coefficients, covar), object$diffemb_basepoint)
    }else{
      diag(nrow = object$n_ambient, ncol = object$n_embedding)
    }
    alignment <- if(with_alignment){
      rotation_map(sum_tangent_vectors(object$alignment_coefficients, covar), diag(nrow = object$n_embedding))
    }else{
      diag(nrow = object$n_embedding)
    }
    approx[,gr == mm_groups] <- diffemb %*% alignment %*% diffemb_embedding[,gr == mm_groups]
  }

  if(with_ambient_pca){
    object$ambient_coordsystem %*% approx + object$ambient_offset
  }else{
    approx
  }
})





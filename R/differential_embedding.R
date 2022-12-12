

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
    Y_clean <- amb_pca$embedding
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
  }else if(! alignment_coef_fixed){
    if(verbose && ! isFALSE(alignment)) message("Align points")
    align_res <- align_points_impl(alignment, diffemb_embedding, design_matrix, verbose = verbose)
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



align_neighbors <- function(fit, neighbors = NULL, method = c("rotation", "stretching", "rotation+stretching"), ...,
                            design_matrix = fit$design_matrix, verbose = TRUE){
  method <- match.arg(method)
  mnn_list <- if(is.null(neighbors)){
    # Over-cluster the logcounts and find mutual nearest neighbors across conditions
    get_mutual_neighbors(assay(fit), design_matrix, ...)
  }else if(is.matrix(neighbors)){
    # Assume neighbors is the corrected embedding from MNN or Harmony
    get_mutual_neighbors(neighbors, design_matrix, ...)
  }else if(is.list(neighbors)){
    if(length(neighbors) != 2) stop("Neighbor list must have two entries")
    if(! is.null(names(neighbors))){
      if(! all(names(neighbors) %in% c("first", "second"))) stop("The names of the neighbors list must be 'first' and 'second'.")
      list(first = as.list(neighbors$first), second = as.list(neighbors$second))
    }else{
      list(first = as.list(neighbors[[1]]), second = as.list(neighbors[[2]]))
    }
  }

  align_res <- align_neighbors_impl(mnn_list, fit$diffemb_embedding, design_matrix, method = method, verbose = verbose, ...)


  reducedDim(fit, "diffemb_embedding") <- t(align_res$diffemb_embedding)
  metadata(fit)[["alignment_coefficients"]] <-  metadata(fit)[["alignment_coefficients"]] + align_res$rotation_coefficients
  metadata(fit)[["bootstrap_samples"]] <- lapply(metadata(fit)[["bootstrap_samples"]], \(samp){
    reducedDim(samp, "diffemb_embedding") <- t(apply_rotation(
      apply_stretching(samp$diffemb_embedding, align_res$stretch_coefficients, design_matrix, align_res$stretch_base_point),
      align_res$rotation_coefficients, design_matrix, align_res$rotation_base_point))
    metadata(samp)[["alignment_coefficients"]] <- metadata(samp)[["alignment_coefficients"]] + align_res$rotation_coefficients
    samp
  })

  fit
}

align_neighbors_impl <- function(mnn_list, diffemb_embedding, design_matrix, method = c("rotation", "stretching", "rotation+stretching"),
                                n_iter = 1, tolerance = 1e-8,  verbose = TRUE, ...){
  method <- match.arg(method)

  n_embedding <- nrow(diffemb_embedding)
  base_point <- diag(nrow = n_embedding)

  d1_per_mnn <- do.call(rbind, lapply(seq_along(mnn_list$first), \(idx){
    design_matrix[mnn_list$first[[idx]][1],]
  }))
  d2_per_mnn <- do.call(rbind, lapply(seq_along(mnn_list$second), \(idx){
    design_matrix[mnn_list$second[[idx]][1],]
  }))
  p1_per_mnn <-  do.call(cbind, lapply(seq_along(mnn_list$first), \(idx){
    rowMeans(diffemb_embedding[, mnn_list$first[[idx]],drop=FALSE])
  }))
  p2_per_mnn <-  do.call(cbind, lapply(seq_along(mnn_list$second), \(idx){
    rowMeans(diffemb_embedding[, mnn_list$second[[idx]],drop=FALSE])
  }))
  w1 <- lengths(mnn_list$first) / (lengths(mnn_list$first) + lengths(mnn_list$second))
  w2 <- lengths(mnn_list$first) / (lengths(mnn_list$first) + lengths(mnn_list$second))
  centers_per_mnn <- t(t(p1_per_mnn) * w2 + t(p2_per_mnn) * w1)
  M <- cbind(centers_per_mnn, centers_per_mnn)
  Y <- cbind(p1_per_mnn, p2_per_mnn)
  D <- rbind(d1_per_mnn, d2_per_mnn)

  if(method == "rotation"){
    rotation_coef <- rotation_lm(M, design = D, obs_embedding = Y, base_point = base_point)
    stretch_coef <- array(0, dim(rotation_coef))
  }else if(method == "stretching"){
    stretch_coef <- spd_lm(M, design = D, obs_embedding = Y, base_point = base_point)
    rotation_coef <- array(0, dim(stretch_coef))
  }else if(method == "rotation+stretching"){
    rotation_coef <- rotation_lm(M, design = D, obs_embedding = Y, base_point = base_point)
    # rotation_coef <- array(0, dim(stretch_coef))
    error <- error_last_round <- mean((Y - M)^2)
    if(verbose) message("Initial error: ", error)
    for(idx in seq_len(10)){
      # Apply **inverse** of rotation to means before fitting stretching
      Mprime <- apply_rotation(M, -rotation_coef, D, base_point)
      stretch_coef <- spd_lm(Mprime, design = D, obs_embedding = Y, base_point = base_point)
      # Stretch the observations before fitting the rotation
      Yprime <- apply_stretching(Y, stretch_coef, D, base_point)
      rotation_coef <- rotation_lm(M, design = D, obs_embedding = Yprime, base_point = base_point)
      # Calculate error
      error <- mean((apply_rotation(apply_stretching(Y, stretch_coef, D, base_point), rotation_coef, D, base_point) - M)^2)
      if(verbose) message("Error: ", error)
      if(abs(error_last_round - error) / (error + 0.5) < tolerance){
        break
      }
      error_last_round <- error
    }
  }

  if(verbose){
    error <- mean((apply_rotation(apply_stretching(Y, stretch_coef, D, base_point), rotation_coef, D, base_point) - M)^2)
    message("Final error: ", error)
  }

  diffemb_embedding <- apply_rotation(
      apply_stretching(diffemb_embedding, stretch_coef, design_matrix, base_point),
    rotation_coef, design_matrix, base_point)


  list(rotation_coefficients = rotation_coef, stretch_coefficients = stretch_coef,
       diffemb_embedding = diffemb_embedding, rotation_base_point = base_point, stretch_base_point= base_point)
}

apply_rotation <- function(A, rotation_coef, design, base_point){
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  groups <- unique(mm_groups)
  for(gr in groups){
    A[,mm_groups == gr] <- rotation_map(sum_tangent_vectors(rotation_coef, design[which(mm_groups == gr)[1],]),
                                        base_point) %*% A[,mm_groups == gr]
  }
  A
}

apply_stretching <- function(A, stretch_coef, design, base_point){
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  groups <- unique(mm_groups)
  for(gr in groups){
    A[,mm_groups == gr] <- spd_map(sum_tangent_vectors(stretch_coef, design[which(mm_groups == gr)[1],]),
                                   base_point) %*% A[,mm_groups == gr]
  }
  A
}




get_mutual_neighbors <- function(data, design_matrix, cell_per_cluster = 20, n_mnn = 10){
  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)
  if(is.null(mm_groups)){
    stop("The model matrix contains too many groups. Is maybe one of the covariates continuous?\n",
         "This error could be removed, but this feature hasn't been implemented yet.")
  }
  groups <- unique(mm_groups)

  mnn_list <- list(first = list(), second = list())
  if(cell_per_cluster <= 1){
    for(idx1 in seq_along(groups)){
      for(idx2 in seq_along(groups)){
        if(idx1 < idx2){
          sel1 <- mm_groups == idx1
          sel2 <- mm_groups == idx2
          k1 <- min(n_mnn, ceiling(sum(sel1) / 4), ceiling(sum(sel2) / 4))
          mnn <- BiocNeighbors::findMutualNN(t(data[,sel1]), t(data[,sel2]), k1 = k1)
          mnn_list$first <- c(mnn_list$first, as.list(which(sel1)[mnn$first]))
          mnn_list$second <- c(mnn_list$second, as.list(which(sel2)[mnn$second]))
        }
      }
    }
  }else{
    clusters <- lapply(groups, \(gr){
      sel <- mm_groups == gr
      k <- max(1, round(sum(sel) / cell_per_cluster))
      cl <- kmeans(t(data[,sel,drop=FALSE]), centers = k, nstart = 10, iter.max = 100)
      list(Y = t(cl$centers), indices = which(sel), cluster = unname(cl$cluster))
    })
    for(idx1 in seq_along(groups)){
      for(idx2 in seq_along(groups)){
        if(idx1 < idx2){
          cl1 <- clusters[[idx1]]
          cl2 <- clusters[[idx2]]
          k1 <- min(n_mnn, ceiling(ncol(cl1$Y) / 4), ceiling(ncol(cl2$Y) / 4))
          mnn <- BiocNeighbors::findMutualNN(t(cl1$Y), t(cl2$Y), k1 = k1)
          mnn_list$first <- c(mnn_list$first, lapply(seq_along(mnn$first), \(idx3){
            cl1$indices[cl1$cluster == mnn$first[idx3]]
          }))
          mnn_list$second <- c(mnn_list$second, lapply(seq_along(mnn$second), \(idx3){
            cl2$indices[cl2$cluster == mnn$second[idx3]]
          }))
        }
      }
    }
  }
  mnn_list
}


#' Enforce additional alignment of cell clusters beyond the direct differential embedding
#'
#' @param alignment a factor of length `ncol(data)`. The method tries to put elements with the
#'   the same factor level close to each other. `NA` entries are ignored.
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#' @param ... additional parameters that are passed on to the internal function `align_points_impl`.
#'   For example, `n_iter`, `tolerance`, and `target_layout`
#'
#'
#' @export
align_embeddings <- function(fit, alignment = TRUE, verbose = TRUE, ...){
  if(isTRUE(alignment)){
    # Cluster each condition
    exp_group <- get_groups(design_matrix, 10)
    exp_group_levels <- unique(exp_group)
    Y_clean <- as.matrix(t(fit$ambient_coordsystem) %*% assay(fit, "expr") - fit$linear_coefficients %*% t(fit$design_matrix))

    clusters <- lapply(exp_group_levels, \(gr) kmeans(t(Y_clean[,exp_group == gr,drop=FALSE]), nstart = 20, centers = 10))
    stop("Not implemented")
  }


  align_res <- align_points_impl(alignment, fit$diffemb_embedding, fit$design_matrix, verbose = verbose, ...)
  metadata(fit)[["alignment_method"]] <- alignment
  reducedDim(fit, "diffemb_embedding") <- t(align_res$diffemb_embedding)
  metadata(fit)[["alignment_coefficients"]] <-  metadata(fit)[["alignment_coefficients"]] + align_res$alignment_coefficients
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
    metadata(samp)[["alignment_coefficients"]] <- metadata(samp)[["alignment_coefficients"]] + align_res$alignment_coefficients
    samp
  })

  fit
}

# Make this function work by solving |M - Beta %*% (P * X)| where M is the mean per k-NN clique,
# P is the current center of the clusters and X is the corresponding line from the design matrix
align_points_impl <- function(alignment, diffemb_embedding, design_matrix,
                              n_iter = 1, tolerance = 1e-8, target_layout = c("mean", "spherical_mds"), verbose = TRUE){

  target_layout <- if(target_layout[1] == "mean"){
    function(embedding, alignment_red, exp_group_red){
      t(mply_dbl(unique(alignment_red), \(al) rowMeans(embedding[,al == alignment_red,drop=FALSE], na.rm = TRUE), ncol = nrow(embedding)))
    }
  }else if(target_layout[1] == "spherical_mds"){
    function(embedding, alignment_red, exp_group_red){
      dists <- lapply(unique(exp_group_red), \(gr) as.matrix(dist_sphere(t(embedding[,exp_group_red == gr,drop=FALSE]))))
      median_dists <- apply(stack_slice(dists), 1:2, median, na.rm = TRUE)
      sphere_mds <- smacof::smacofSphere(median_dists, ndim = min(nrow(embedding), ncol(median_dists)-1), eps = 1e-10, itmax = 1e4)
      mds_layout <- normalize_vec_length(rbind(t(sphere_mds$conf), matrix(0, nrow = nrow(embedding) - ncol(sphere_mds$conf), ncol = ncol(median_dists))))
      n_emb <-  nrow(mds_layout)
      reflec_mat <- diag(c(rep(1, n_emb - 1), -1), nrow = n_emb)
      rot1 <- rotation_map(drop(rotation_lm(data = embedding[,exp_group_red == exp_group_red[1],drop=FALSE], design = matrix(1, nrow = ncol(mds_layout)),
                                            obs_embedding = mds_layout, base_point = diag(nrow = n_emb))), base_point = diag(nrow = n_emb))
      rot2 <- rotation_map(drop(rotation_lm(data = embedding[,exp_group_red == exp_group_red[1],drop=FALSE],
                                                design = matrix(1, nrow = ncol(mds_layout)), obs_embedding = reflec_mat %*% mds_layout,
                                                base_point = diag(nrow = n_emb))), base_point = diag(nrow = n_emb))
      if(sum((embedding[,exp_group_red == exp_group_red[1],drop=FALSE] - rot1 %*% mds_layout)^2) <
         sum((embedding[,exp_group_red == exp_group_red[1],drop=FALSE] - rot2  %*% reflec_mat %*% mds_layout)^2)){
        rot1 %*% mds_layout
      }else{
        rot2 %*% reflec_mat %*% mds_layout
      }
    }
  }else if(is.matrix(target_layout)){
    tl <- target_layout
    function(...){
      tl
    }
  }else{
    stopifnot(is.function(target_layout))
    target_layout
  }

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
    new_embedding <- normalize_vec_length(new_embedding)
    if(verbose) check_alignment_prerequisites(new_embedding, alignment_red)

    # Calculate the mean per alignment level and replicate for each experimental design group
    mean_embeddings <- target_layout(new_embedding, alignment_red, exp_group_red) |>
      duplicate_cols(each = n_exp_group_levels) |>
      normalize_vec_length()

    coef <- array(0, dim  = c(n_embedding, n_embedding, ncol(design_matrix)))
    new_embedding_cp <- new_embedding
    last_round_error <- Inf
    for(iter in seq_len(n_iter)){

      # Align data to centers
      coef <- coef + rotation_lm(data = mean_embeddings,  design = new_design, obs_embedding = new_embedding,
                                 base_point = base_point, tangent_regression = TRUE)

      # Rotate centers
      for(gr in exp_group_levels){
        new_embedding[,gr == exp_group_red] <- rotation_map(sum_tangent_vectors(coef, design_matrix[which(exp_group == gr)[1],]), base_point) %*% new_embedding_cp[,gr == exp_group_red]
      }

      error <- sum((new_embedding - mean_embeddings)^2, na.rm = TRUE)
      if(verbose) message("Iter: ", iter, "\tError: ", error)
      if(abs(error - last_round_error) / (error + 0.1) < tolerance){
        break
      }
      last_round_error <- error
    }

    for(gr in exp_group_levels){
      new_embedding[,gr == exp_group_red] <- rotation_map(sum_tangent_vectors(coef, design_matrix[which(exp_group == gr)[1],]), base_point) %*% new_embedding_cp[,gr == exp_group_red]
    }
    if(verbose) message("Final result\tError: ", sum((new_embedding - mean_embeddings)^2, na.rm = TRUE))

    # Rotate all datapoints
    for(gr in exp_group_levels){
      diffemb_embedding[,gr == exp_group] <- rotation_map(sum_tangent_vectors(coef, design_matrix[which(exp_group == gr)[1],]), base_point) %*% diffemb_embedding[,gr == exp_group]
    }
  }

  list(alignment_coefficients = -coef, diffemb_embedding = diffemb_embedding)
}


check_alignment_prerequisites <- function(embedding, alignment){
  n_obs <- ncol(embedding) / length(unique(alignment))
  dist_var <- matrixStats::colSds(mply_dbl(unique(alignment), \(al){
    stopifnot(sum(alignment == al) == n_obs)
    c(as.matrix(dist(t(embedding[,alignment == al,drop=FALSE]))))
  }, ncol = n_obs^2), na.rm = TRUE)
  df <- data.frame(obs1 = rep(seq(n_obs), times = n_obs), obs2 = rep(seq(n_obs), each = n_obs),
             dist_var = dist_var, row.names = NULL, stringsAsFactors = FALSE)
  res <- tapply(df$dist_var, df$obs1, mean, na.rm = TRUE)
  res <- data.frame(group = seq_len(n_obs), distance_from_diagonal = res)
  message("Distance from diagonal: \n",
          paste0(capture.output(res[order(-res$distance_from_diagonal)[seq_len(min(4, n_obs))],,drop=FALSE]), collapse = "\n"))
  invisible(res)
}






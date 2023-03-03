#' Enforce additional alignment of cell clusters beyond the direct differential embedding
#'
#' @param fit a `DiffEmbSeqFit` object
#' @param design_matrix the design matrix for the alignment
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#' @param ... additional parameters that are passed on to relevant functions
#'
#'
#' @export
align_neighbors <- function(fit, rotating = TRUE, stretching = TRUE,
                            data_matrix = assay(fit), cells_per_cluster = 20, mnn = 10,
                            design = fit$alignment_design_matrix, ridge_penalty = 0, verbose = TRUE){
  if(verbose) message("Find mutual nearest neighbors")
  design_matrix <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))$design_matrix
  mnn_groups <- get_mutual_neighbors(data_matrix, design_matrix, cells_per_cluster = cells_per_cluster, mnn = mnn)
  # if(verbose) message("Adjust latent positions using a '", method, "' transformation")
  correction <- correct_design_matrix_groups(fit, mnn_groups, fit$embedding, design, rotating = rotating, stretching = stretching, ridge_penalty = ridge_penalty)
  correct_fit(fit, correction)
}

#' @rdname align_neighbors
#' @export
align_harmony <- function(fit, rotating = TRUE, stretching = TRUE, ...,
                          design = fit$alignment_design_matrix,
                          ridge_penalty = 0, min_cluster_membership = 0.1, max_iter = 10, verbose = TRUE){
  if(verbose) message("Select cells that are considered close with 'harmony'")
  design_matrix <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))$design_matrix
  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)
  if(! requireNamespace("harmony", quietly = TRUE)){
    stop("'harmony' is not installed. Please install it from CRAN.")
  }
  # Ignore best practice and call private methods from harmony
  harm_obj <- harmony_init(fit$embedding, design_matrix, ..., verbose = verbose)
  for(idx in seq_len(max_iter)){
    harm_obj <- harmony_max_div_clustering(harm_obj)
    matches <- lapply(seq_len(harm_obj$K), \(row_idx) which(harm_obj$R[row_idx,] > min_cluster_membership))
    weights <- lapply(seq_len(harm_obj$K), \(row_idx) harm_obj$R[row_idx,matches[[row_idx]]])
    index_groups <- lapply(matches, \(idx) mm_groups[idx])
    # if(verbose) message("Adjust latent positions using a '", method, "' transformation")
    correction <- correct_design_matrix_groups(fit, list(matches = matches, index_groups = index_groups, weights = weights),
                                               harm_obj$Z_orig, design, rotating = rotating, stretching = stretching, ridge_penalty = ridge_penalty)
    harm_obj$Z_corr <- correction$embedding
    harm_obj$Z_cos <- t(t(harm_obj$Z_corr) / sqrt(colSums(harm_obj$Z_corr^2)))
    if(harm_obj$check_convergence(1)){
      if(verbose) message("Converged")
      break
    }
  }

  correct_fit(fit, correction)
}


#' @rdname align_neighbors
#' @export
align_by_template <- function(fit, rotating = TRUE, stretching = TRUE,
                              alignment_template, cells_per_cluster = 20, mnn = 10,
                              design = fit$alignment_design_matrix, ridge_penalty = 0, verbose = TRUE){
  stopifnot(is.matrix(alignment_template))
  stopifnot(ncol(alignment_template) == ncol(fit))
  if(verbose) message("Received template that puts similar cells close to each other")
  align_neighbors(fit, rotating = rotating, stretching = stretching, data_matrix = alignment_template,
                  cells_per_cluster = cells_per_cluster, mnn = mnn,
                  design = design, ridge_penalty = ridge_penalty, verbose = verbose)
}

#' @rdname align_neighbors
#' @export
align_by_grouping <- function(fit, rotating = TRUE, stretching = TRUE,
                              grouping, design = fit$alignment_design_matrix, ridge_penalty = 0, verbose = TRUE){
  if(verbose) message("Received sets of cells that are considered close")
  if(is.list(grouping)){
    # Check that it conforms to the expectation of the mnn_grouping
    if("matches" %in% names(grouping)){
      # Do nothing
    }else{
      grouping <- list(matches = grouping)
    }
    valid <- vapply(grouping$matches, \(entry){
      is.numeric(entry) && all(is.na(entry) | (entry > 0 & entry <= ncol(fit)))
    }, FUN.VALUE = logical(1L))
    stopifnot(all(valid))
  }else if(is.vector(grouping) || is.factor(grouping)){
    stopifnot(length(grouping) == ncol(fit))
    # Convert grouping into the right shape
    unique_values <- unique(grouping)
    unique_values <- unique_values[! is.na(unique_values)]
    matches <- lapply(unique_values, \(lvl) which(grouping == lvl))
    grouping <- list(matches = matches)
  }else{
    stop("Cannot handle 'grouping' of class", paste0(class(grouping), collapse = ", "))
  }

  if(!"index_groups" %in% names(grouping)){
    design_matrix <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))$design_matrix
    mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)
    grouping$index_groups <- lapply(grouping$matches, \(idx) mm_groups[idx])
  }

  correction <- correct_design_matrix_groups(fit, grouping, fit$embedding, design, rotating = rotating, stretching = stretching, ridge_penalty = ridge_penalty)
  correct_fit(fit, correction)
}


#' Find the rotation and stretching coefficients to make latent position of matching groups similar
#'
#' @param matching_groups a list with two (plus one optional) elements.
#'   \describe{
#'     \item{matches}{a list of vectors with the indices of the cells that are part of the match}
#'     \item{index_groups}{a list of vectors whose lengths correspond to `match`. The entries correspond to
#'     the separate conditions in the design matrix. *This is weird because it is redundant with the `design` argument.*}
#'     \item{weights}{an optional list of vectors whose lengths correspond to `match`. The entries correspond to
#'     the weight of that cell in the match.}
#'   }
#'
#' @keywords internal
correct_design_matrix_groups <- function(fit, matching_groups, embedding, design, rotating = TRUE, stretching = TRUE,
                                         ridge_penalty = 0, verbose = TRUE){

  n_embedding <- nrow(embedding)
  base_point <- diag(nrow = n_embedding)
  ridge_penalty <- handle_ridge_penalty_parameter(ridge_penalty)

  des <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))
  design_matrix <- des$design_matrix
  if(is.null(des$design_formula) && matrix_equals(fit$design_matrix, design_matrix)){
    design_formula <- fit$design
  }else{
    design_formula <- des$design_formula
  }
  D <- do.call(rbind, lapply(seq_along(matching_groups$matches), \(idx){
    do.call(rbind, lapply(unique(matching_groups$index_groups[[idx]]), \(igr){
      design_matrix[matching_groups$matches[[idx]][which(matching_groups$index_groups[[idx]] == igr)[1]],,drop=FALSE]
    }))
  }))
  Y <- do.call(cbind, lapply(seq_along(matching_groups$matches), \(idx){
    do.call(cbind, lapply(unique(matching_groups$index_groups[[idx]]), \(igr){
      if(is.null(matching_groups$weights)){
        matrixStats::rowMeans2(embedding, cols = matching_groups$matches[[idx]][matching_groups$index_groups[[idx]] == igr])
      }else{
        matrixStats::rowWeightedMeans(embedding[,matching_groups$matches[[idx]],drop=FALSE], w = matching_groups$weights[[idx]], cols = matching_groups$index_groups[[idx]] == igr)
      }
    }))
  }))

  # M is the mean of the aggregated groups Y
  M <- do.call(cbind, lapply(seq_along(matching_groups$matches), \(idx){
    Y_tmp <- do.call(cbind, lapply(unique(matching_groups$index_groups[[idx]]), \(igr){
      if(is.null(matching_groups$weights)){
        matrixStats::rowMeans2(embedding, cols = matching_groups$matches[[idx]][matching_groups$index_groups[[idx]] == igr])
      }else{
        matrixStats::rowWeightedMeans(embedding[,matching_groups$matches[[idx]],drop=FALSE], w = matching_groups$weights[[idx]], cols = matching_groups$index_groups[[idx]] == igr)
      }
    }))
    vec <- if(is.null(Y_tmp)){
      matrix(nrow = nrow(embedding), ncol = 0)
    }else{
      rowMeans(Y_tmp)
    }
    duplicate_cols(vec, length(unique(matching_groups$index_groups[[idx]])))
  }))

  weights <- if(is.null(matching_groups$weights)){
    unlist(lapply(seq_along(matching_groups$matches), \(idx) vapply(unique(matching_groups$index_groups[[idx]]), \(igr){
      sum(matching_groups$index_groups[[idx]] == igr)
    }, FUN.VALUE = 0.0)))
  }else{
    unlist(lapply(seq_along(matching_groups$matches), \(idx) vapply(unique(matching_groups$index_groups[[idx]]), \(igr){
      sum(matching_groups$weights[[idx]][matching_groups$index_groups[[idx]] == igr])
    }, FUN.VALUE = 0.0)))
  }


  if(rotating && ! stretching){
    rotation_coef <- rotation_lm(M, design = D, obs_embedding = Y, base_point = base_point, ridge_penalty = ridge_penalty$rotation, weights = weights)
    stretch_coef <- array(0, dim(rotation_coef))
  }else if(stretching && ! rotating){
    stretch_coef <- spd_lm(M, design = D, obs_embedding = Y, base_point = base_point, ridge_penalty = ridge_penalty$stretching, weights = weights)
    rotation_coef <- array(0, dim(stretch_coef))
  }else if(stretching && rotating){
    plm <- polar_lm(M, D, Y, base_point = base_point, ridge_penalty = ridge_penalty, weights = weights)
    stretch_coef <- plm$stretch_coef
    rotation_coef <- plm$rotation_coef
  }else{
    stretch_coef <- array(0, c(nrow(Y), nrow(Y), ncol(D)))
    rotation_coef <- array(0, dim(stretch_coef))
  }

  embedding <- apply_stretching(embedding, stretch_coef, design_matrix, base_point)
  embedding <- apply_rotation(embedding, rotation_coef, design_matrix, base_point)


  list(rotation_coefficients = -rotation_coef, stretch_coefficients = -stretch_coef,
       embedding = embedding, design_matrix = design_matrix, design_formula = design_formula)
}


polar_lm <- function(data, design, obs_embedding, base_point, ridge_penalty = 0, weights = NULL, max_iter = 10, tolerance = 1e-8){
  ridge_penalty <- handle_ridge_penalty_parameter(ridge_penalty)

  if(length(unique(get_groups(design, n_groups = ncol(design) * 10))) == ncol(design)){
    # Clever initialization
    tryCatch({
      pla <- polar_lm_analytic(data, design, obs_embedding, base_point)
      rotation_coef <- pla$rotation_coef
      stretch_coef <- pla$stretch_coef
    }, error = function(e){
      # Error might occur if some group has too few observations
      rotation_coef <<- rotation_lm(data, design = design, obs_embedding = obs_embedding, base_point = base_point, ridge_penalty = ridge_penalty$rotation, weights = weights)
      stretch_coef <<- array(0, dim(rotation_coef))
    })
  }else{
    # Naive initialization
    rotation_coef <- rotation_lm(data, design = design, obs_embedding = obs_embedding, base_point = base_point, ridge_penalty = ridge_penalty$rotation, weights = weights)
    stretch_coef <- array(0, dim(rotation_coef))
  }

  error <- error_last_round <- mean((apply_rotation(apply_stretching(obs_embedding, stretch_coef, design, base_point), rotation_coef, design, base_point) - data)^2)
  for(idx in seq_len(max_iter)){
    # Apply **inverse** of rotation to means before fitting stretching
    Mprime <- apply_rotation(data, -rotation_coef, design, base_point)
    stretch_coef <- spd_lm(Mprime, design = design, obs_embedding = obs_embedding, base_point = base_point, ridge_penalty = ridge_penalty$stretching, weights = weights)
    # Stretch the observations before fitting the rotation
    Yprime <- apply_stretching(obs_embedding, stretch_coef, design, base_point)
    rotation_coef <- rotation_lm(data, design = design, obs_embedding = Yprime, base_point = base_point, ridge_penalty = ridge_penalty$rotation, weights = weights)
    # Calculate error
    error <- mean((apply_rotation(apply_stretching(obs_embedding, stretch_coef, design, base_point), rotation_coef, design, base_point) - data)^2)
    if((is.na(error) || is.na(error_last_round)) || abs(error_last_round - error) / (error + 0.5) < tolerance){
      # Error can be NaN if nrow(embedding) == 0
      break
    }
    error_last_round <- error
  }
  list(rotation_coef = rotation_coef, stretch_coef = stretch_coef)
}

polar_lm_analytic <-  function(data, design, obs_embedding, base_point){
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  groups <- unique(mm_groups)

  if(length(groups) == ncol(design)){
    coefs <- lapply(groups, \(gr){
      # Simple linear regression
      beta <- t(coef(lm.fit(t(obs_embedding[,mm_groups == gr,drop=FALSE]), t(data[,mm_groups == gr,drop=FALSE]))))
      if(any(is.na(beta))){
        stop("Group ", gr, " contains too few observations")
      }else{
        # Polar decomposition where I force U to be a rotation
        svd <- svd(beta)
        diag_elem <- c(rep(1, times = ncol(beta) - 1), Matrix::det(svd$u %*% t(svd$v)))
        U <- svd$u %*% diag(diag_elem, nrow = length(diag_elem))  %*% t(svd$v)
        P <- project_psd(coef(lm.fit(U, beta)))

        list(rotation_coef = rotation_log(base_point, U), stretch_coef = spd_log(base_point, P))
      }
    })
    rot_coefs_matrix <- do.call(cbind, lapply(coefs, \(e) c(e$rotation_coef)))
    stretch_coefs_matrix <- do.call(cbind, lapply(coefs, \(e) c(e$stretch_coef)))
    artificial_design <- as.matrix(Matrix::sparseMatrix(i = seq_along(mm_groups), j = mm_groups, x = rep(1, length(mm_groups))))
    change_mat <- t(coef(lm.fit(design, artificial_design)))
    rot_coefs_matrix <- rot_coefs_matrix %*% change_mat
    stretch_coefs_matrix <- stretch_coefs_matrix %*% change_mat

    dims <- c(nrow(data), nrow(data), ncol(design))
    dimnames <- list(NULL, NULL, colnames(design))
    list(rotation_coef = array(rot_coefs_matrix, dims, dimnames), stretch_coef = array(stretch_coefs_matrix, dims, dimnames))
  }else{
    stop("More groups than design columns")
  }
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


get_mutual_neighbors <- function(data, design_matrix, cells_per_cluster = 20, mnn = 10){
  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)
  if(is.null(mm_groups)){
    stop("The model matrix contains too many groups. Is maybe one of the covariates continuous?\n",
         "This error could be removed, but this feature hasn't been implemented yet.")
  }
  groups <- unique(mm_groups)

  matches <- list()
  index_group <- list()
  if(cells_per_cluster <= 1){
    # Iterate over combinations of design matrix row groups
    for(idx1 in seq_along(groups)){
      for(idx2 in seq_along(groups)){
        if(idx1 < idx2){
          sel1 <- mm_groups == idx1
          sel2 <- mm_groups == idx2
          k1 <- min(mnn, ceiling(sum(sel1) / 4), ceiling(sum(sel2) / 4))
          mnn_list <- BiocNeighbors::findMutualNN(t(data[,sel1]), t(data[,sel2]), k1 = k1)
          matches <- c(matches, lapply(seq_along(mnn_list$first), \(idx3){
            c(which(sel1)[mnn_list$first[idx3]], which(sel2)[mnn_list$second[idx3]])
          }))
          index_group <- c(index_group, rep(list(c(1,2)), length(mnn_list$first)))
        }
      }
    }
  }else{
    # Cluster each design matrix row group
    clusters <- lapply(groups, \(gr){
      sel <- mm_groups == gr
      k <- max(1, round(sum(sel) / cells_per_cluster))
      cl <- kmeans(t(data[,sel,drop=FALSE]), centers = k, nstart = 10, iter.max = 100)
      list(Y = t(cl$centers), indices = which(sel), cluster = unname(cl$cluster))
    })
    # Iterate over combinations of design matrix row groups
    for(idx1 in seq_along(groups)){
      for(idx2 in seq_along(groups)){
        if(idx1 < idx2){
          # Find MNN for each cluster center
          cl1 <- clusters[[idx1]]
          cl2 <- clusters[[idx2]]
          k1 <- min(mnn, ceiling(ncol(cl1$Y) / 4), ceiling(ncol(cl2$Y) / 4))
          mnn_list <- BiocNeighbors::findMutualNN(t(cl1$Y), t(cl2$Y), k1 = k1)
          matches <- c(matches, lapply(seq_along(mnn_list$first), \(idx3){
            c(cl1$indices[cl1$cluster == mnn_list$first[idx3]], cl2$indices[cl2$cluster == mnn_list$second[idx3]])
          }))
          index_group <- c(index_group, lapply(seq_along(mnn_list$first), \(idx3){
            c(rep(1, length(cl1$indices[cl1$cluster == mnn_list$first[idx3]])),
              rep(2, length(cl2$indices[cl2$cluster == mnn_list$second[idx3]])))
          }))
        }
      }
    }
  }
  list(matches = matches, index_groups = index_group)
}



correct_fit <- function(fit, correction){
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  reducedDim(fit, "embedding") <- t(correction$embedding)
  if(! matrix_equals(correction$design_matrix, fit$alignment_design_matrix) ||
     is.null(correction$design_formula) != is.null(fit$alignment_design) ||
     correction$design_formula != fit$alignment_design){
    metadata(fit)[["alignment_rotation"]] <- correction$rotation_coefficients
    metadata(fit)[["alignment_stretching"]] <- correction$stretch_coefficients
    metadata(fit)[["alignment_design_matrix"]] <- correction$design_matrix
    metadata(fit)[["alignment_design"]] <- correction$design_formula
  }else{
    metadata(fit)[["alignment_rotation"]] <-  metadata(fit)[["alignment_rotation"]] + correction$rotation_coefficients
    metadata(fit)[["alignment_stretching"]] <-  metadata(fit)[["alignment_stretching"]] + correction$stretch_coefficients
  }
  fit
}

handle_ridge_penalty_parameter <- function(ridge_penalty){
  ridge_penalty <- lapply(ridge_penalty, \(x) x * 1.0)
  if(length(ridge_penalty) == 2){
    if(is.null(names(ridge_penalty))){
      names(ridge_penalty) <- c("rotation", "stretching")
    }else if(! all(names(ridge_penalty) %in% c("rotation", "stretching"))){
      stop("Names of ridge penalty must be 'rotation' and 'stretching'. They are: ", paste0(names(ridge_penalty), collapse = ", "))
    }
  }else if(length(ridge_penalty) == 1){
    if(is.null(names(ridge_penalty))){
      ridge_penalty <- list(rotation = ridge_penalty[[1]], stretching = ridge_penalty[[1]])
    }else{
      if(names(ridge_penalty) == "rotation"){
        ridge_penalty <- c(ridge_penalty, list(stretching = 0))
      }else if(names(ridge_penalty) == "stretching"){
        ridge_penalty <- c(ridge_penalty, list(rotation = 0))
      }else{
        stop("Cannot handle name of ridge penalty: ", names(ridge_penalty))
      }
    }
  }else{
    stop("Cannot handle ridge penalty of length: ", length(ridge_penalty))
  }
  ridge_penalty
}

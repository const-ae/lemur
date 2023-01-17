#' Enforce additional alignment of cell clusters beyond the direct differential embedding
#'
#' @param fit a `DiffEmbSeqFit` object
#' @param method how flexible should the alignment be
#' @param design_matrix the design matrix for the alignment
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#' @param ... additional parameters that are passed on to relevant functions
#'
#'
#' @export
align_neighbors <- function(fit, method = c("rotation", "stretching", "rotation+stretching"),
                            data_matrix = assay(fit), cells_per_cluster = 20, mnn = 10,
                            design_matrix = fit$design_matrix, verbose = TRUE){
  method <- match.arg(method)
  if(verbose) message("Find mutual nearest neighbors")
  mnn_groups <- get_mutual_neighbors(data_matrix, design_matrix, cells_per_cluster = cells_per_cluster, mnn = mnn)
  if(verbose) message("Adjust latent positions using a '", method, "' transformation")
  correction <- correct_design_matrix_groups(mnn_groups, fit$diffemb_embedding, design_matrix, method = method)
  correct_fit(fit, correction)
}

#' @rdname align_neighbors
#' @export
align_harmony <- function(fit, method = c("rotation", "stretching", "rotation+stretching"), ...,
                          design_matrix = fit$design_matrix, max_iter = 10, verbose = TRUE){
  method <- match.arg(method)
  if(verbose) message("Select cells that are considered close with 'harmony'")
  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)
  if(! requireNamespace("harmony", quietly = TRUE)){
    stop("'harmony' is not installed. Please install it from CRAN.")
  }
  # Ignore best practice and call private methods from harmony
  harm_obj <- harmony_init(fit$diffemb_embedding, design_matrix, verbose = verbose)
  for(idx in seq_len(max_iter)){
    harm_obj <- harmony_max_div_clustering(harm_obj)
    matches <- apply(harm_obj$R, 1, \(row) which(row > 0.1))
    index_groups <- lapply(matches, \(idx) mm_groups[idx])
    if(verbose) message("Adjust latent positions using a '", method, "' transformation")
    correction <- correct_design_matrix_groups(list(matches = matches, index_groups = index_groups),
                                               harm_obj$Z_orig, design_matrix, method = method)
    harm_obj$Z_corr <- correction$diffemb_embedding
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
align_scvi <- function(fit, method = c("rotation", "stretching", "rotation+stretching"), counts = assay(fit, "counts"),
                       n_hidden = 128, n_layers = 10,  dropout_rate = 1,
                       dispersion = c("gene", "gene-batch", "gene-label", "gene-cell"),
                       gene_likelihood = c("nb", "zinb", "poisson"), latent_distribution = c("normal", "ln"),
                       cells_per_cluster = 20, mnn = 10,
                       design_matrix = fit$design_matrix, verbose = TRUE){
  method <- match.arg(method)
  dispersion <- match.arg(dispersion)
  gene_likelihood <- match.arg(gene_likelihood)
  latent_distribution <- match.arg(latent_distribution)

  if(verbose) message("Select cells that are considered close with 'SCVI'")
  if(requireNamespace("reticulate", quietly = TRUE)){
    sc <- reticulate::import("scanpy", convert = FALSE)
    scvi <- reticulate::import("scvi", convert = FALSE)

    obs <- data.frame(condition = apply(fit$design_matrix, 1, \(row) paste0(row, collapse = "-")))

    adata <- sc$AnnData(X = t(as.matrix(assay(fit))),
                        obs = obs,
                        layers = list(counts = t(as.matrix(counts))))
    scvi$model$SCVI$setup_anndata(adata, layer = "counts", batch_key = "condition")

    vae <- scvi$model$SCVI(adata, n_hidden = as.integer(n_hidden), n_layers = as.integer(n_layers), dropout_rate = as.integer(dropout_rate),
                           dispersion = dispersion, gene_likelihood = gene_likelihood, latent_distribution = latent_distribution)
    vae$train()
    scvi_lat <- t(as.matrix(vae$get_latent_representation()))
  }else{
    stop("Cannot call 'SCVI' because the R package 'reticulate' is not available to call python functions.")
  }
  align_by_template(fit, method = method, alignment_template = scvi_lat,
                    cells_per_cluster = cells_per_cluster, mnn = mnn, design_matrix = design_matrix, verbose = verbose)
}

#' @rdname align_neighbors
#' @export
align_by_template <- function(fit, method = c("rotation", "stretching", "rotation+stretching"),
                              alignment_template, cells_per_cluster = 20, mnn = 10,
                              design_matrix = fit$design_matrix, verbose = TRUE){
  method <- match.arg(method)
  stopifnot(is.matrix(alignment_template))
  stopifnot(ncol(alignment_template) == ncol(fit))
  if(verbose) message("Received template that puts similar cells close to each other")
  align_neighbors(fit, method = method, data_mat = alignment_template,
                  cells_per_cluster = cells_per_cluster, mnn = mnn,
                  design_matrix = design_matrix, verbose = verbose)
}

#' @rdname align_neighbors
#' @export
align_by_grouping <- function(fit, method = c("rotation", "stretching", "rotation+stretching"),
                              grouping,
                              design_matrix = fit$design_matrix, verbose = TRUE){
  method <- match.arg(method)
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
    mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)
    grouping$index_groups <- lapply(matches, \(idx) mm_groups[idx])
  }

  correction <- correct_design_matrix_groups(grouping, fit$diffemb_embedding, design_matrix, method = method)
  correct_fit(fit, correction)
}



correct_design_matrix_groups <- function(matching_groups, diffemb_embedding, design_matrix, method = c("rotation", "stretching", "rotation+stretching"),
                                         max_iter = 10, n_iter = 1, tolerance = 1e-8,  verbose = TRUE, ...){
  method <- match.arg(method)

  n_embedding <- nrow(diffemb_embedding)
  base_point <- diag(nrow = n_embedding)

  D <- do.call(rbind, lapply(seq_along(matching_groups$matches), \(idx){
    do.call(rbind, lapply(unique(matching_groups$index_groups[[idx]]), \(igr){
      design_matrix[matching_groups$matches[[idx]][which(matching_groups$index_groups[[idx]] == igr)[1]],,drop=FALSE]
    }))
  }))
  Y <- do.call(cbind, lapply(seq_along(matching_groups$matches), \(idx){
    do.call(cbind, lapply(unique(matching_groups$index_groups[[idx]]), \(igr){
      rowMeans(diffemb_embedding[, matching_groups$matches[[idx]][matching_groups$index_groups[[idx]] == igr],drop=FALSE])
    }))
  }))
  M <- do.call(cbind, lapply(seq_along(matching_groups$matches), \(idx){
    duplicate_cols(rowMeans(diffemb_embedding[, matching_groups$matches[[idx]],drop=FALSE]), length(unique(matching_groups$index_groups[[idx]])))
  }))

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
    for(idx in seq_len(max_iter)){
      # Apply **inverse** of rotation to means before fitting stretching
      Mprime <- apply_rotation(M, -rotation_coef, D, base_point)
      stretch_coef <- spd_lm(Mprime, design = D, obs_embedding = Y, base_point = base_point)
      # Stretch the observations before fitting the rotation
      Yprime <- apply_stretching(Y, stretch_coef, D, base_point)
      rotation_coef <- rotation_lm(M, design = D, obs_embedding = Yprime, base_point = base_point)
      # Calculate error
      error <- mean((apply_rotation(apply_stretching(Y, stretch_coef, D, base_point), rotation_coef, D, base_point) - M)^2)
      if((is.na(error) || is.na(error_last_round)) || abs(error_last_round - error) / (error + 0.5) < tolerance){
        # Error can be NaN if nrow(diffemb_embedding) == 0
        break
      }
      error_last_round <- error
    }
  }

  diffemb_embedding <- apply_rotation(
    apply_stretching(diffemb_embedding, stretch_coef, design_matrix, base_point),
    rotation_coef, design_matrix, base_point)


  list(rotation_coefficients = -rotation_coef, stretch_coefficients = -stretch_coef,
       diffemb_embedding = diffemb_embedding, design_matrix = design_matrix)
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

  reducedDim(fit, "diffemb_embedding") <- t(correction$diffemb_embedding)
  if(! all(dim(correction$design_matrix) == dim(fit$alignment_design_matrix)) ||
     ! all(correction$design_matrix == fit$alignment_design_matrix)){
    metadata(fit)[["alignment_rotation"]] <- correction$rotation_coefficients
    metadata(fit)[["alignment_stretching"]] <- correction$stretch_coefficients
    metadata(fit)[["alignment_design_matrix"]] <- correction$design_matrix
  }else{
    metadata(fit)[["alignment_rotation"]] <-  metadata(fit)[["alignment_rotation"]] + correction$rotation_coefficients
    metadata(fit)[["alignment_stretching"]] <-  metadata(fit)[["alignment_stretching"]] + correction$stretch_coefficients
  }
  metadata(fit)[["bootstrap_samples"]] <- lapply(metadata(fit)[["bootstrap_samples"]], \(samp){
    reducedDim(samp, "diffemb_embedding") <- t(apply_rotation(
      apply_stretching(samp$diffemb_embedding, -correction$stretch_coefficients, correction$design_matrix, diag(nrow = samp$n_embedding)),
      -correction$rotation_coefficients, correction$design_matrix, diag(nrow = samp$n_embedding)))
    if(! all(dim(correction$design_matrix) == dim(samp$alignment_design_matrix)) ||
       ! all(correction$design_matrix == samp$alignment_design_matrix)){
      metadata(samp)[["alignment_rotation"]] <- correction$rotation_coefficients
      metadata(samp)[["alignment_stretching"]] <- correction$stretch_coefficients
      metadata(samp)[["alignment_design_matrix"]] <-   correction$design_matrix
    }else{
      metadata(samp)[["alignment_rotation"]] <- metadata(samp)[["alignment_rotation"]] + correction$rotation_coefficients
      metadata(samp)[["alignment_stretching"]] <- metadata(samp)[["alignment_stretching"]] + correction$stretch_coefficients
    }
    samp
  })
  fit
}



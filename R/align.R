#' Enforce additional alignment of cell clusters beyond the direct differential embedding
#'
#' @param fit a `lemur_fit` object
#' @param design a specification of the design (matrix or formula) that is used
#'   for the transformation. Default: `fit$design_matrix`
#' @param ridge_penalty specification how much the flexibility of the transformation
#'   should be regularized. This can be a single positive scalar or named list
#'   (`"stretching"` and `"rotation"`). Default: `0.01`
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#' @param ... additional parameters that are passed on to relevant functions
#' @param cells_per_cluster argument specific for `align_neighbors`. Before mutual nearest neighbor
#'   search the cells are clustered per condition. Default: `20`
#' @param mnn argument specific for `align_neighbors`. The number of mutual nearest neighbors.
#'   Default: `10`
#' @param data_matrix argument specific for `align_neighbors`. The data matrix that is used for
#'   the mutual nearest neighbor search. Default: `assay(fit)`
#' @param min_cluster_membership argument specific for `align_harmony`. The minimum probability
#'   from the soft clustering. Default: `0.1`
#' @param max_iter argument specific for `align_harmony`. The number of iterations. Default: `10`
#' @param template argument specific for `align_by_template`. A readily integrated dataset
#'   that is approximated using the `align_neighbors` function.
#' @param grouping argument specific for `align_by_grouping`. A manually provided list of groups
#'   that are considered similar.
#'
#'
#' @export
align_neighbors <- function(fit,
                            data_matrix = assay(fit), cells_per_cluster = 20, mnn = 10,
                            design = fit$alignment_design_matrix, ridge_penalty = 0.01, verbose = TRUE){
  if(verbose) message("Find mutual nearest neighbors")
  design_matrix <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))$design_matrix
  mnn_groups <- get_mutual_neighbors(data_matrix, design_matrix, cells_per_cluster = cells_per_cluster, mnn = mnn)
  # if(verbose) message("Adjust latent positions using a '", method, "' transformation")
  correction <- correct_design_matrix_groups(fit, mnn_groups, fit$embedding, design, ridge_penalty = ridge_penalty)
  correct_fit(fit, correction)
}

#' @rdname align_neighbors
#' @export
align_harmony <- function(fit, ...,
                          design = fit$alignment_design_matrix,
                          ridge_penalty = 0.01, min_cluster_membership = 0.001, max_iter = 10, verbose = TRUE){
  if(verbose) message("Select cells that are considered close with 'harmony'")
  if(is.null(attr(design, "ignore_degeneracy"))){
    # It doesn't matter for harmony if the design is degenerate
    attr(design, "ignore_degeneracy") <- TRUE
  }
  design_matrix <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))$design_matrix
  mm_groups <- get_groups(design_matrix, n_groups = ncol(design_matrix) * 10)
  if(! requireNamespace("harmony", quietly = TRUE)){
    stop("'harmony' is not installed. Please install it from CRAN.")
  }
  # Ignore best practice and call private methods from harmony
  harm_obj <- harmony_init(fit$embedding, design_matrix, ..., verbose = verbose)
  for(idx in seq_len(max_iter)){
    harm_obj <- harmony_max_div_clustering(harm_obj)
    threshold <- min(min_cluster_membership, max(harm_obj$R) * 0.5)
    matches <- lapply(seq_len(harm_obj$K), \(row_idx) which(harm_obj$R[row_idx,] > threshold))
    weights <- lapply(seq_len(harm_obj$K), \(row_idx) harm_obj$R[row_idx,matches[[row_idx]]])
    index_groups <- lapply(matches, \(idx) mm_groups[idx])
    # if(verbose) message("Adjust latent positions using a '", method, "' transformation")
    correction <- correct_design_matrix_groups(fit, list(matches = matches, index_groups = index_groups, weights = weights),
                                               harm_obj$Z_orig, design, ridge_penalty = ridge_penalty)
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
align_by_template <- function(fit, template,
                              cells_per_cluster = 20, mnn = 10,
                              design = fit$alignment_design_matrix, ridge_penalty = 0.01, verbose = TRUE){
  stopifnot(is.matrix(alignment_template))
  stopifnot(ncol(alignment_template) == ncol(fit))
  if(verbose) message("Received template that puts similar cells close to each other")
  align_neighbors(fit, data_matrix = alignment_template,
                  cells_per_cluster = cells_per_cluster, mnn = mnn,
                  design = design, ridge_penalty = ridge_penalty, verbose = verbose)
}

#' @rdname align_neighbors
#' @export
align_by_grouping <- function(fit,
                              grouping, design = fit$alignment_design_matrix, ridge_penalty = 0.01, verbose = TRUE){
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

  correction <- correct_design_matrix_groups(fit, grouping, fit$embedding, design, ridge_penalty = ridge_penalty)
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
correct_design_matrix_groups <- function(fit, matching_groups, embedding, design, ridge_penalty = 0.01, verbose = TRUE){

  n_embedding <- nrow(embedding)
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

  interact_design_matrix <- duplicate_cols(D, each = n_embedding)  * duplicate_cols(t(Y), times = ncol(D))
  alignment_coefs <- ridge_regression(M - Y, X = interact_design_matrix, ridge_penalty = ridge_penalty)
  alignment_coefs <- array(alignment_coefs, dim = c(n_embedding, n_embedding, ncol(D)))

  embedding <- apply_linear_transformation(embedding, alignment_coefs, design_matrix)


  list(alignment_coefs = alignment_coefs, embedding = embedding, design_matrix = design_matrix, design_formula = design_formula)
}


forward_linear_transformation <- function(alignment_coefficients, design_vector){
  diag(nrow = dim(alignment_coefficients)[1]) + sum_tangent_vectors(alignment_coefficients, design_vector)
}

reverse_linear_transformation <- function(alignment_coefficients, design_vector){
  n_embedding <- dim(alignment_coefficients)[1]
  if(n_embedding == 0){
    matrix(nrow = 0, ncol = 0)
  }else{
    solve(diag(nrow = n_embedding) + sum_tangent_vectors(alignment_coefficients, design_vector))
  }
}

apply_linear_transformation <- function(A, alignment_coefs, design){
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  groups <- unique(mm_groups)
  for(gr in groups){
    A[,mm_groups == gr] <- forward_linear_transformation(alignment_coefs,  design[which(mm_groups == gr)[1],])  %*% A[,mm_groups == gr]
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
    metadata(fit)[["alignment_coefficients"]] <- correction$alignment_coefs
    metadata(fit)[["alignment_design_matrix"]] <- correction$design_matrix
    metadata(fit)[["alignment_design"]] <- correction$design_formula
  }else{
    metadata(fit)[["alignment_coefficients"]] <- correction$alignment_coefs
  }
  fit
}

handle_ridge_penalty_parameter <- function(ridge_penalty){
  if(any(names(ridge_penalty) %in% c("rotation", "stretching"))){
    stop("The alignment function has changed and the rotation and stretching specification is now defunct")
  }
  ridge_penalty
}

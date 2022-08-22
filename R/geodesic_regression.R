
#' Solve d(P, exp_p(V * x))^2 for V
#'
#'
grassmann_geodesic_regression <- function(coordsystems, design, base_point, tangent_regression = FALSE){
  # Validate input
  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)
  n_emb <- ncol(base_point)

  coordsystems <- if(is.list(coordsystems)){
    coordsystems
  }else if(is.array(coordsystems)){
    stopifnot(length(dim(coordsystems)) == 3)
    destack_slice(coordsystems)
  }else{
    stop("Cannot handle coordsystems of type: ", paste0(class(coordsystems), collapse = ", "))
  }
  stopifnot(length(coordsystems) == nrow(design))
  stopifnot(all(vapply(coordsystems, \(emb) nrow(emb) == n_amb && ncol(emb) == n_emb, FUN.VALUE = logical(1L))))
  # stopifnot(all(vapply(coordsystems, \(emb) is_grassmann_element(emb), FUN.VALUE = logical(1L))))


  # Initialize with tangent regression (if possible)
  tangent_vecs <- lapply(coordsystems, \(emb) as.vector(grassmann_log(base_point, emb)))
  merged_vecs <- stack_cols(tangent_vecs)
  tangent_fit <- t(lm.fit(design, t(merged_vecs))$coefficients)
  coef <- stack_slice(lapply(seq_len(ncol(tangent_fit)), \(idx) matrix(tangent_fit[,idx], nrow = n_amb, ncol = n_emb)))
  dimnames(coef) <- list(NULL, NULL, colnames(tangent_fit))


  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}

#' Solve ||Y - exp_p(V * x) Y ||^2_2 for V
#'
#'
grassmann_lm <- function(data, design, base_point, tangent_regression = FALSE){
  nas <- apply(data, 2, anyNA) | apply(design, 1, anyNA)
  data <- data[,!nas,drop=FALSE]
  design <- design[!nas,,drop=FALSE]

  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)
  n_emb <- ncol(base_point)

  # Initialize with tangent regression
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  if(is.null(mm_groups)){
    stop("The model matrix contains too many groups. Is maybe one of the covariates continuous?\n",
         "This error could be removed, but this feature hasn't been implemented yet.")
  }
  if(any(table(mm_groups) < n_emb)){
    stop("Too few datapoints in some design matrix group.\n",
         "This error could be removed, but this feature hasn't been implemented yet.")
  }
  groups <- unique(mm_groups)
  reduced_design <- mply_dbl(groups, \(gr) design[which(mm_groups == gr)[1],], ncol = ncol(design))
  group_planes <- lapply(groups, \(gr) pca(data[,mm_groups == gr,drop=FALSE], n = n_emb)$coordsystem)
  coef <- grassmann_geodesic_regression(group_planes, design = reduced_design, base_point = base_point, tangent_regression = TRUE)
  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}

#' Solve d(R, exp_p(V * x))^2 for V
#'
#'
rotation_geodesic_regression <- function(rotations, design, base_point, tangent_regression = FALSE){
  # Validate input
  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)

  rotations <- if(is.list(rotations)){
    rotations
  }else if(is.array(rotations)){
    stopifnot(length(dim(rotations)) == 3)
    destack_slice(rotations)
  }else{
    stop("Cannot handle rotations of type: ", paste0(class(rotations), collapse = ", "))
  }
  stopifnot(length(rotations) == nrow(design))
  stopifnot(all(vapply(rotations, \(rot) nrow(rot) == n_amb && ncol(rot) == n_amb, FUN.VALUE = logical(1L))))
  # stopifnot(all(vapply(rotations, \(rot) is_rotation_element(rot), FUN.VALUE = logical(1L))))

  # TODO: If only two rotations, do procrustes to solve the problem.

  # Initialize with tangent regression (if possible)
  tangent_vecs <- lapply(rotations, \(rot){
    # TODO: only extract upper triangular matrix, because
    # the other elements are redundant.
    as.vector(rotation_log(base_point, rot))
  })
  merged_vecs <- stack_cols(tangent_vecs)
  tangent_fit <- t(lm.fit(design, t(merged_vecs))$coefficients)
  tangent_fit[is.na(tangent_fit)] <- 0
  coef <- stack_slice(lapply(seq_len(ncol(tangent_fit)), \(idx) matrix(tangent_fit[,idx], nrow = n_amb, ncol = n_amb)))
  dimnames(coef) <- list(NULL, NULL, colnames(tangent_fit))


  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}

#' Solve ||Y - exp_p(V * x) Z ||^2_2 for V (aka. Procrustes regression)
#'
#' Here data = t(grassmann_map(V * X)) Y
#'
rotation_lm <- function(data, design, obs_embedding, base_point, tangent_regression = FALSE){
  nas <- apply(data, 2, anyNA) | apply(design, 1, anyNA) | apply(obs_embedding, 2, anyNA)
  data <- data[,!nas,drop=FALSE]
  design <- design[!nas,,drop=FALSE]
  obs_embedding <- obs_embedding[,!nas,drop=FALSE]

  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)

  # Initialize with tangent regression
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  if(is.null(mm_groups)){
    stop("The model matrix contains too many groups. Is maybe one of the covariates continuous?\n",
         "This error could be removed, but this feature hasn't been implemented yet.")
  }
  # if(any(table(mm_groups) < n_emb)){
  #   stop("Too few datapoints in some design matrix group.\n",
  #        "This error could be removed, but this feature hasn't been implemented yet.")
  # }
  groups <- unique(mm_groups)
  reduced_design <- mply_dbl(groups, \(gr) design[which(mm_groups == gr)[1],], ncol = ncol(design))
  group_rot <- lapply(groups, \(gr){
    sel <- mm_groups == gr
    x <- obs_embedding[,sel,drop=FALSE]
    y <- data[,sel,drop=FALSE]
    data_rank <- min(qr(x)$rank, qr(y)$rank)
    if(data_rank == 0){
      stop("Not sure what to do now")
    }else if(data_rank < n_amb){
      # The system is underdetermined.

      # The rotation will be somewhere in the space spanned by x and y
      qr_space <- qr.Q(qr(cbind(x, y)), complete = TRUE)
      n_red <- min(data_rank * 2, n_amb)
      space <- qr_space[,seq_len(n_red),drop=FALSE]
      anti_space <- qr_space[,-seq_len(n_red),drop=FALSE]

      # Approximate the rotation using linear regression
      coef_lower_dim <- t(lm.fit(round(t(x) %*% space, 10), round(t(y) %*% space, 10))$coefficient)
      coef_is_na <- apply(coef_lower_dim, 2, anyNA)

      # Tricky bit! [x,y] has 2n dim, but we only have n observations.
      # To fill the remaining n, we want to find a minimal rotation (i.e. close to identify matrix).
      # We project the subset of the identity to the null-space of the coefficients
      coef_null_space_proj <- diag(nrow = n_red) - coef_lower_dim[,!coef_is_na,drop=FALSE] %*% t(coef_lower_dim[,!coef_is_na,drop=FALSE])
      coef_lower_dim[,coef_is_na] <- coef_null_space_proj  %*% diag(nrow = n_red)[,setdiff(seq_len(n_red), seq(sum(!coef_is_na))),drop=FALSE]

      # The following is a modified version of `project_rotation`, where I
      # resolve the indeterminancy of the SVD and force it to use the anti-space
      # This ensures that points in the null space of [x, y] are not moved by the rotation
      svd <- svd(coef_lower_dim)
      U <- cbind(space %*% svd$u, anti_space)
      V <- cbind(space %*% svd$v, anti_space)
      diag_elem <- c(rep(1, times = ncol(U) - 1), Matrix::det(U %*% t(V)))
      U %*% diag(diag_elem, nrow = length(diag_elem)) %*% t(V)
    }else{
      # coef <- t(lm.fit(t(x), t(y))$coefficient)
      # project_rotation(coef)
      procrustes_rotation(y, x)
    }
  })
  coef <- rotation_geodesic_regression(group_rot, design = reduced_design, base_point = base_point, tangent_regression = TRUE)
  # line search
  original_error <- sum(vapply(groups, \(gr){
    sel <- mm_groups == gr
    sum((data[,sel,drop=FALSE] - obs_embedding[,sel,drop=FALSE])^2)
  }, FUN.VALUE = 0.0))
  for(idx in 0:20){
    error <- sum(vapply(groups, \(gr){
      sel <- mm_groups == gr
      sum((data[,sel,drop=FALSE] - rotation_map(sum_tangent_vectors(0.5^idx * coef, reduced_design[groups == gr, ]), base_point) %*% obs_embedding[,sel,drop=FALSE])^2)
    }, FUN.VALUE = 0.0))
    if(error < original_error){
      coef <- 0.5^idx * coef
      break
    }
  }
  if(idx == 20) coef <- 0 * coef
  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}

#' Solve ||Y - R Z ||^2_2 for R in SO(n) (aka. Procrustes analysis)
#'
procrustes_rotation <- function(data, obs_embedding){
  n <- nrow(obs_embedding)
  stopifnot(nrow(data) == n)
  stopifnot(ncol(obs_embedding) == ncol(data))
  zy_dim <- qr(cbind(data, obs_embedding))$rank
  if(zy_dim < n){
    stop("This branch isn't working yet.")
    # There is extra space between outside [z,y] which
    # is not at all affected by the rotation
    qr_space <- qr.Q(qr(cbind(obs_embedding, data)), complete = TRUE)
    space <- qr_space[,seq_len(zy_dim),drop=FALSE]
    anti_space <- qr_space[,-seq_len(zy_dim),drop=FALSE]
    red_svd <- svd(t(space) %*% obs_embedding %*% t(data) %*% space)
  # browser()
    # sv0 <- red_svd$d < 1e-18
    # red_svd$u[sv0,sv0] <- diag(nrow=sum(sv0))

    U <- cbind(space %*% red_svd$u, anti_space)
    V <- cbind(space %*% red_svd$v, anti_space)
  }else{
    svd <- svd(obs_embedding %*% t(data))
    U <- svd$u
    V <- svd$v
  }
  diag_elem <- c(rep(1, times = ncol(U) - 1), Matrix::det(V %*% t(U)))
  V %*% diag(diag_elem, nrow = length(diag_elem)) %*% t(U)

}


get_groups <- function (model_matrix, n_groups) {
  if (!glmGamPoi:::lte_n_equal_rows(model_matrix, n_groups)) {
    NULL
  } else {
    glmGamPoi:::get_row_groups(model_matrix, n_groups = n_groups)
  }
}





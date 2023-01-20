
#######################
# Grassmann Manifold  #
#######################


#' Solve d(P, exp_p(V * x))^2 for V
#'
#'
grassmann_geodesic_regression <- function(coordsystems, design, base_point, weights = rep(1, length(coordsystems)), tangent_regression = FALSE){
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
  stopifnot(length(coordsystems) == n_obs)
  stopifnot(all(vapply(coordsystems, \(emb) nrow(emb) == n_amb && ncol(emb) == n_emb, FUN.VALUE = logical(1L))))
  # stopifnot(all(vapply(coordsystems, \(emb) is_grassmann_element(emb), FUN.VALUE = logical(1L))))
  stopifnot(length(weights) == n_obs)


  # Initialize with tangent regression (if possible)
  tangent_vecs <- lapply(coordsystems, \(emb) as.vector(grassmann_log(base_point, emb)))
  merged_vecs <- stack_cols(tangent_vecs)
  tangent_fit <- if(nrow(merged_vecs) == 0){
    matrix(nrow = 0, ncol = ncol(design))
  }else{
    t(lm.wfit(design, t(merged_vecs), w = weights)$coefficients)
  }
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
  group_planes <- lapply(groups, \(gr) pca(data[,mm_groups == gr,drop=FALSE], n = n_emb, center = FALSE)$coordsystem)
  elem_per_group <- vapply(groups, \(gr) sum(mm_groups == gr), FUN.VALUE = 0L)
  coef <- grassmann_geodesic_regression(group_planes, design = reduced_design, base_point = base_point, weights = elem_per_group, tangent_regression = TRUE)
  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}

######################
# Rotation Manifold  #
######################


#' Solve d(R, exp_p(V * x))^2 for V
#'
#'
rotation_geodesic_regression <- function(rotations, design, base_point, weights = 1, ridge_penalty = ridge_penalty, tangent_regression = FALSE){
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
  weights <- rep_len(weights, nrow(design))

  # Initialize with tangent regression
  tangent_vecs <-  lapply(rotations, \(rot){
    c(rotation_log(base_point, rot)[upper.tri(rot)])
  })
  merged_vecs <- stack_cols(tangent_vecs)
  if(nrow(merged_vecs) == 0){
    tangent_fit <- matrix(nrow = 0, ncol = ncol(merged_vecs))
  }else{
    if(all(ridge_penalty == 0)){
      tangent_fit <- t(lm.wfit(design, t(merged_vecs), w = weights)$coefficients)
    }else{
      tangent_fit <- ridge_regression(merged_vecs, X = design, ridge_penalty = ridge_penalty, weights = weights)
    }
  }
  tangent_fit[is.na(tangent_fit)] <- 0
  coef <- stack_slice(lapply(seq_len(ncol(tangent_fit)), \(idx){
    res <- matrix(0, nrow = n_amb, ncol = n_amb)
    res[upper.tri(res)] <- tangent_fit[,idx]
    res - t(res)
  }))
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
rotation_lm <- function(data, design, obs_embedding, base_point, ridge_penalty = 0, tangent_regression = FALSE){
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
  groups <- unique(mm_groups)
  reduced_design <- mply_dbl(groups, \(gr) design[which(mm_groups == gr)[1],], ncol = ncol(design))
  group_rot <- lapply(groups, \(gr){
    sel <- mm_groups == gr
    procrustes_rotation(data[,sel,drop=FALSE], obs_embedding[,sel,drop=FALSE])
  })
  group_sizes <- vapply(groups, \(gr) sum(mm_groups == gr), FUN.VALUE = integer(1L))
  coef <- rotation_geodesic_regression(group_rot, design = reduced_design, base_point = base_point, weights = group_sizes, ridge_penalty = ridge_penalty, tangent_regression = TRUE)
  # line search
  # original_error <- sum(vapply(groups, \(gr){
  #   sel <- mm_groups == gr
  #   sum((data[,sel,drop=FALSE] - obs_embedding[,sel,drop=FALSE])^2)
  # }, FUN.VALUE = 0.0))
  # for(idx in 0:20){
  #   if(idx != 0) stop("Halving step size")
  #   error <- sum(vapply(groups, \(gr){
  #     sel <- mm_groups == gr
  #     sum((data[,sel,drop=FALSE] - rotation_map(sum_tangent_vectors(0.5^idx * coef, reduced_design[groups == gr, ]), base_point) %*% obs_embedding[,sel,drop=FALSE])^2)
  #   }, FUN.VALUE = 0.0))
  #   if(error < original_error){
  #     coef <- 0.5^idx * coef
  #     break
  #   }
  # }
  # if(idx == 20) coef <- 0 * coef
  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}

#' Solve ||Y - R Z ||^2_2 for R in SO(n) (aka. Procrustes analysis)
#'
procrustes_rotation <- function(Y, Z){
  # This code is based on an answer by Mike Hawk on cross validated
  # https://stats.stackexchange.com/a/599015/130486
  stopifnot(nrow(Y) == nrow(Z))
  stopifnot(ncol(Y) == ncol(Z))

  if(nrow(Y) == 0){
    matrix(nrow = 0, ncol = 0)
  }else if(nrow(Y) == 1){
    matrix(1)
  }else if(nrow(Y) == 2){
    a <- sum(Z * Y)
    b <- sum(Y[1,] * Z[2,] - Y[2,] * Z[1,])
    theta <- atan2(a, b) - pi/2
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2)
  }else{
    M <- Z %*% t(Y)
    svd <- svd(M)
    sel <- abs(svd$d) > 1e-12
    if(any(! sel)){
      # Cross product is not full rank, recurse on subproblem
      U1 <- svd$u[,sel,drop=FALSE]
      V1 <- svd$v[,sel,drop=FALSE]
      U2 <- svd$u[,!sel,drop=FALSE]
      V2 <- svd$v[,!sel,drop=FALSE]
      fact <- Matrix::det(svd$v %*% t(svd$u))
      U2[,ncol(U2)] <- U2[,ncol(U2)] * fact
      sub_rot <- procrustes_rotation(t(V2), t(U2))
      Ucompl <- cbind(U1, U2 %*% t(sub_rot))
      Vcompl <- cbind(V1, V2)
      Vcompl %*% t(Ucompl)
    }else{
      U <- svd$u
      V <- svd$v
      fact <- Matrix::det(V %*% t(U))
      U[,ncol(U)] <- U[,ncol(U)] * fact
      V %*% t(U)
    }
  }
}


################################################
# Symmetric Positive Definite Matrix Manifold  #
################################################

#' Solve d(P, exp_p(V * x))^2 for V
#'
spd_geodesic_regression <- function(spd_matrices, design, base_point, weights = 1, ridge_penalty = 0, tangent_regression = FALSE){
  # Validate input
  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)

  spd_matrices <- if(is.list(spd_matrices)){
    spd_matrices
  }else if(is.array(spd_matrices)){
    stopifnot(length(dim(spd_matrices)) == 3)
    destack_slice(spd_matrices)
  }else{
    stop("Cannot handle SPD matrices of type: ", paste0(class(spd_matrices), collapse = ", "))
  }
  stopifnot(length(spd_matrices) == nrow(design))
  stopifnot(all(vapply(spd_matrices, \(spd) nrow(spd) == n_amb && ncol(spd) == n_amb, FUN.VALUE = logical(1L))))
  # stopifnot(all(vapply(spd_matrices, \(spd) is_spd_element(spd), FUN.VALUE = logical(1L))))
  weights <- rep_len(weights, nrow(design))

  # Initialize with tangent regression (if possible)
  tangent_vecs <- lapply(spd_matrices, \(spd){
    c(spd_log(base_point, spd)[upper.tri(spd, diag = TRUE)])
  })
  merged_vecs <- stack_cols(tangent_vecs)
  if(nrow(merged_vecs) == 0){
    tangent_fit <- matrix(nrow = 0, ncol = ncol(merged_vecs))
  }else{
    if(all(ridge_penalty == 0)){
      tangent_fit <- t(lm.wfit(design, t(merged_vecs), w = weights)$coefficients)
    }else{
      tangent_fit <- ridge_regression(merged_vecs, X = design, ridge_penalty = ridge_penalty, weights = weights)
    }
  }
  tangent_fit[is.na(tangent_fit)] <- 0
  coef <- stack_slice(lapply(seq_len(ncol(tangent_fit)), \(idx){
    res <- matrix(0, nrow = n_amb, ncol = n_amb)
    res[upper.tri(res, diag = TRUE)] <- tangent_fit[,idx]
    diag(res) <- diag(res) / 2
    res + t(res)
  }))
  dimnames(coef) <- list(NULL, NULL, colnames(tangent_fit))


  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}


#' Solve ||Y - exp_p(V * x) Z ||^2_2 for V (aka. SPD Procrustes regression)
#'
#' Here data = t(grassmann_map(V * X)) Y
#'
spd_lm <- function(data, design, obs_embedding, base_point, ridge_penalty = 0, tangent_regression = FALSE){
  nas <- apply(data, 2, anyNA) | apply(design, 1, anyNA) | apply(obs_embedding, 2, anyNA)
  data <- data[,!nas,drop=FALSE]
  design <- design[!nas,,drop=FALSE]
  obs_embedding <- obs_embedding[,!nas,drop=FALSE]

  n_obs <- nrow(design)
  n_coef <- ncol(design)
  n_amb <- nrow(base_point)
  if(n_amb == 0){
    return(array(dim = c(0, 0, n_coef)))
  }

  # Initialize with tangent regression
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  if(is.null(mm_groups)){
    stop("The model matrix contains too many groups. Is maybe one of the covariates continuous?\n",
         "This error could be removed, but this feature hasn't been implemented yet.")
  }

  groups <- unique(mm_groups)
  reduced_design <- mply_dbl(groups, \(gr) design[which(mm_groups == gr)[1],], ncol = ncol(design))
  group_spd <- lapply(groups, \(gr){
    sel <- mm_groups == gr
    x <- obs_embedding[,sel,drop=FALSE]
    y <- data[,sel,drop=FALSE]
    procrustes_spd(y, x)
  })
  group_sizes <- vapply(groups, \(gr) sum(mm_groups == gr), FUN.VALUE = integer(1L))

  coef <- spd_geodesic_regression(group_spd, design = reduced_design, base_point = base_point, weights = group_sizes, ridge_penalty = ridge_penalty, tangent_regression = TRUE)
  # line search
  original_error <- sum(vapply(groups, \(gr){
    sel <- mm_groups == gr
    sum((data[,sel,drop=FALSE] - obs_embedding[,sel,drop=FALSE])^2)
  }, FUN.VALUE = 0.0))
  for(idx in 0:20){
    error <- sum(vapply(groups, \(gr){
      sel <- mm_groups == gr
      sum((data[,sel,drop=FALSE] - spd_map(sum_tangent_vectors(0.5^idx * coef, reduced_design[groups == gr, ]), base_point) %*% obs_embedding[,sel,drop=FALSE])^2)
    }, FUN.VALUE = 0.0))
    if(error < original_error){
      coef <- 0.5^idx * coef
      break
    }
  }
  if(idx == 20){
    warning("SPD Procrustes regression did not converge. Set all coefficients to zero.")
    coef <- 0 * coef
  }
  if(tangent_regression){
    coef
  }else{
    # warning("Refine regression using Riemannian optimization. (Not yet implemented)")
    coef
  }
}





#' Solve ||B - A X ||^2_2 for A in SPD(n) (aka. Symmetric Positive Definite Procrustes analysis)
#'
procrustes_spd <- function(data, obs_embedding, maxiter = 1000, tolerance = 1e-8){
  # Implementation based on 'A semi-analytical approach for the positive semidefinite Procrustes problem'
  # by Gillis et al. (2018)
  # Code adapted from https://sites.google.com/site/nicolasgillis/code?authuser=0#h.p_ID_160
  n <- nrow(obs_embedding)
  stopifnot(nrow(data) == n)
  stopifnot(ncol(obs_embedding) == ncol(data))

  A <- init_procrustes_spd(data, obs_embedding)

  # Precomputing
  XXt <- obs_embedding %*% t(obs_embedding)
  eig <- eigen(XXt)
  x <- diag(eig$values, nrow = n)
  u <- eig$vectors
  Lx <- max(eig$values)
  mux <- min(eig$values)
  qx <- mux / Lx
  BXt <- data %*% t(obs_embedding)

  # Parameters and initialization
  alpha0 <- 0.1
  alpha <- alpha0
  beta <- numeric(0L)
  Y <- A
  i <- 1
  eps_last_round <- 0
  eps <- Inf

  for(i in  seq_len(maxiter)){
    Ap <- A
    alpha_i <- alpha[i]
    alpha[i+1] <- (sqrt((alpha_i^2 - qx)^2 + 4 * alpha_i^2) + (qx - alpha_i^2)) / 2
    beta[i] <- alpha_i * (1 - alpha_i) / (alpha_i^2 + alpha[i + 1])
    A <- project_psd(Y - (Y %*% XXt - BXt) / Lx)
    Y <- A + beta[i] * (A - Ap)

    eps_last_round <- eps
    eps <- sum((A - Ap)^2)
    if(abs(eps_last_round - eps) / (abs(eps) + 0.5) < tolerance){
      break
    }
  }
  A
}

init_procrustes_spd <- function(data, obs_embedding){
  n <- nrow(obs_embedding)
  # Init A (three different options)
  A1 <- project_psd(t(pracma::mldivide(t(obs_embedding), t(data))))
  A2 <- matrix(0, nrow = n, ncol = n)
  for(i in seq_len(n)){
    A2[i,i] <- max(0, sum(obs_embedding[i,] * data[i,]) / (sum(obs_embedding[i,]^2) + 1e-6))
  }
  A3 <- tryCatch({
    # This might fail on non-full rank input
    analytic_approx_procrustes_spd(data, obs_embedding)
  }, error = function(err){
    matrix(0, nrow = n, ncol = n)
  })

  e1 <- sum((A1 %*% obs_embedding - data)^2)
  e2 <- sum((A2 %*% obs_embedding - data)^2)
  e3 <- sum((A3 %*% obs_embedding - data)^2)
  if(e1 < e2 && e1 < e3){
    A1
  }else if(e2 < e3){
    A2
  }else{
    A3
  }
}

#' Solve ||Y - P Z ||^2_2 for P in SPD(n) (aka. Symmetric Positive Definite Procrustes analysis)
#'
analytic_approx_procrustes_spd <- function(data, obs_embedding){
  # Implementation based on 'SOLUTION OF SYMMETRIC POSITIVE SEMIDEFINITE PROCRUSTES PROBLEM' by Peng et al. (2019)
  # Seems to produce good but not necessarily optimal solutions
  n <- nrow(obs_embedding)
  stopifnot(nrow(data) == n)
  stopifnot(ncol(obs_embedding) == ncol(data))

  # Step 2: Eq. (1.2)
  z_svd <- svd(obs_embedding)
  sel <- abs(z_svd$d) >= 1e-10
  # The rank of Z
  r <- sum(sel)
  U1 <- z_svd$u[,sel,drop=FALSE]
  U2 <- z_svd$u[,!sel,drop=FALSE]
  V1 <- z_svd$v[,sel,drop=FALSE]
  sigma <- z_svd$d[sel]

  # Step 3: Eq. (2.21)
  phi_pre <- matrix(sigma^2, nrow = r, ncol = r)
  phi <- 1/(phi_pre + t(phi_pre))
  S_pre <- t(U1) %*% data %*% V1 %*% diag(sigma, nrow = r)
  S_hat <- phi * (S_pre + t(S_pre))
  eigen_decomp <- eigen(S_hat)


  # Step 4: Eq. (2.22)
  P11 <- with(eigen_decomp, vectors %*% diag(pmax(values, 0)) %*% t(vectors))

  if(r < n){
    # Step 5:
    P12 <- diag(1/sigma, nrow = r) %*% t(V1) %*% t(data) %*% U2
    # Step 6:
    if(qr(cbind(P11, P12))$rank != qr(P11)$rank) stop("The procrustres problem does not have a solution")
    # Step 7 (assuming P22 = 0)
    P22 <- t(P12) %*% solve(t(P11) %*% P11) %*% t(P11) %*% P12
    P <- rbind(cbind(P11, P12),
               cbind(t(P12), P22))
  }else{
    P <- P11
  }
  z_svd$u %*% P %*% t(z_svd$u)
}


project_psd <- function(Q){
  Q <- (Q + t(Q)) / 2
  eig <- eigen(Q)
  with(eig, vectors %*% diag(pmax(values, 1e-12), nrow = length(values)) %*% t(vectors))
}



get_groups <- function (design_matrix, n_groups) {
  if (!glmGamPoi:::lte_n_equal_rows(design_matrix, n_groups)) {
    NULL
  } else {
    glmGamPoi:::get_row_groups(design_matrix, n_groups = n_groups)
  }
}





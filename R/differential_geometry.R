
# Grassmann (Gr(n, k))

grassmann_map <- function(x, base_point){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L93
  svd <- svd(x)
  z <- base_point %*% svd$v %*% diag(cos(svd$d), nrow = length(svd$d)) %*% t(svd$v) +
    svd$u %*% diag(sin(svd$d), nrow = length(svd$d)) %*% t(svd$v)
  # qr.Q(qr(z))
  z
}

grassmann_log <- function(p, q){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L174
  n <- nrow(p)
  k <- ncol(p)
  stopifnot(nrow(q) == n, ncol(q) == k)

  z <- t(q) %*% p
  At <- t(q) - z %*% t(p)
  Bt <- lm.fit(z, At)$coefficients
  svd <- svd(t(Bt), k, k)
  svd$u %*% diag(atan(svd$d), nrow = k) %*% t(svd$v)
}

project_grassmann <- function(x){
  qr.Q(qr(x))
}

project_grassmann_tangent <- function(x, base_point){
  x - base_point %*% t(base_point) %*% x
}


random_grassmann_point <- function(n, k, ...){
  V <- randn(n, k, ...)
  project_grassmann(V)
}

random_grassmann_tangent <- function(p, ...){
  n <- nrow(p)
  k <- ncol(p)
  Z <- randn(n, k, ...)
  X <- project_grassmann_tangent(p, Z)
  norm_X <- sqrt(sum(X^2))
  X / norm_X
}

# Rotations (SO(n))

project_rotation <- function(x){
  svd <- svd(x)
  diag_elem <- c(rep(1, times = ncol(x) - 1), Matrix::det(svd$u %*% t(svd$v)))
  svd$u %*% diag(diag_elem, nrow = length(diag_elem)) %*% t(svd$v)
}

project_rotation_tangent <- function(x, base_point){
  skew(x)
}

rotation_map <- function(x, base_point){
  base_point %*% expm::expm(x)
}

rotation_log <- function(p, q){
  suppressWarnings({
    logm <- expm::logm(t(p) %*% q)
  })
  if(all(is.na(logm))){
    # The Higham08 algorithm failed. Try Eigen
    logm <- expm::logm(t(p) %*% q, method = "Eigen")
  }

  skew(logm)
}


# # Stiefel
#
# stiefel_map <- function(x, base_point){
#   # Implementation based on https://juliamanifolds.github.io/Manifolds.jl/latest/manifolds/stiefel.html#Base.exp-Tuple{Stiefel,%20Vararg{Any,%20N}%20where%20N}
#   n <- ncol(base_point)
#   ptx <- t(base_point) %*% x
#   I_n <- diag(nrow = n)
#   Zero_n <- matrix(0, nrow = n, ncol = n)
#   block <- cbind(rbind(ptx, I_n), rbind(-t(x) %*% x), ptx)
#   rbind(p, x) %*% expm::expm(block) %*% rbind(expm::expm(-ptx), Zero_n)
# }
#
# project_stiefel <- function(x){
#   qr.Q(qr(x))
# }


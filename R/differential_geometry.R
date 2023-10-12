
# Grassmann (Gr(n, k))

grassmann_map <- function(x, base_point){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L93
  if(ncol(base_point) == 0 || nrow(base_point) == 0){
    base_point
  }else if(any(is.na(x))){
    matrix(NA, nrow = nrow(x), ncol = ncol(x))
  }else{
    svd <- svd(x)
    z <- base_point %*% svd$v %*% diag(cos(svd$d), nrow = length(svd$d)) %*% t(svd$v) +
      svd$u %*% diag(sin(svd$d), nrow = length(svd$d)) %*% t(svd$v)
    # Calling `qr.Q(qr(z))` is problematic because it can flip the signs
    z
  }
}

grassmann_log <- function(p, q){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L174
  # The Grassmann manifold handbook proposes an alternative algorithm in section 5.2
  n <- nrow(p)
  k <- ncol(p)
  stopifnot(nrow(q) == n, ncol(q) == k)
  if(n == 0 || k == 0){
    p
  }else{
    z <- t(q) %*% p
    At <- t(q) - z %*% t(p)
    Bt <- lm.fit(z, At)$coefficients
    svd <- svd(t(Bt), k, k)
    svd$u %*% diag(atan(svd$d), nrow = k) %*% t(svd$v)
  }
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
  project_grassmann_tangent(Z, p)
}


grassmann_angle_from_tangent <- function(x, normalized = TRUE){
  # Conversion of tangent vector to angle taken from Proposition 5.1 of the Grassmann manifold handbook
  thetas <- (svd(x)$d / pi * 180)
  if(normalized){
    thetas <- thetas %% 180
    max(pmin(thetas, 180 - thetas))
  }else{
    thetas[1]
  }
}

grassmann_angle_from_points <- function(p, q){
  grassmann_angle_from_tangent(grassmann_log(p, q))
}



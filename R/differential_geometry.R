
# Grassmann (Gr(n, k))

grassmann_map <- function(x, base_point){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L93
  if(ncol(base_point) == 0 || nrow(base_point) == 0){
    base_point
  }else{
    svd <- svd(x)
    z <- base_point %*% svd$v %*% diag(cos(svd$d), nrow = length(svd$d)) %*% t(svd$v) +
      svd$u %*% diag(sin(svd$d), nrow = length(svd$d)) %*% t(svd$v)
    # qr.Q(qr(z))
    z
  }
}

grassmann_log <- function(p, q){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L174
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
  # norm_X <- sqrt(sum(X^2))
  # X / norm_X
}

# Rotations (SO(n))

project_rotation <- function(x){
  if(nrow(x) == 0){
    x
  }else{
    svd <- svd(x)
    diag_elem <- c(rep(1, times = ncol(x) - 1), Matrix::det(svd$u %*% t(svd$v)))
    svd$u %*% diag(diag_elem, nrow = length(diag_elem)) %*% t(svd$v)
  }
}

project_rotation_tangent <- function(x, base_point){
  skew(x)
}

rotation_map <- function(x, base_point){
  base_point %*% expm::expm(x)
}

rotation_log <- function(p, q){
  if(nrow(p) == 0){
    p
  }else{
    suppressWarnings({
      logm <- expm::logm(t(p) %*% q)
    })
    if(all(is.na(logm))){
      # The Higham08 algorithm failed. Try Eigen
      logm <- tryCatch({
         expm::logm(t(p) %*% q, method = "Eigen")
      }, error = function(error){
        for(idx in 1:10){
          tmp <- tryCatch({
            expm::logm(t(p) %*% q + randn(nrow(p), ncol(p), sd = 1e-12), method = "Eigen")
          }, error = function(error2) NULL)
          if(is.null(tmp)){
            # try once more
          }else{
            break
          }
        }
        tmp
      })
      if(is.null(logm)){
        stop("Cannot calculate the matrix logarithm. It may not exist")
      }
    }

    skew(logm)
  }
}


random_rotation_point <- function(n, ...){
  V <- randn(n, n, ...)
  project_rotation(V)
}

random_rotation_tangent <- function(p, ...){
  n <- nrow(p)
  stopifnot(n == ncol(p))
  Z <- randn(n, n, ...)
  project_rotation_tangent(Z, p)
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

# Sphere


project_sphere <- function(x){
  x
}

project_sphere_tangent <- function(x, base_point){
  base_point <- base_point / sqrt(sum(base_point^2))
  x - drop(t(base_point) %*% x) * base_point
}

sphere_map <- function(x, base_point){
  radius <- sqrt(sum(base_point^2))
  base_point <- base_point  / radius
  norm_x <- sqrt(sum(x^2))
  if(abs(norm_x < 1e-18)){
    radius * (cos(norm_x) * base_point + x)
  }else{
    radius * (cos(norm_x) * base_point + sin(norm_x) / norm_x * x )
  }
}

sphere_log <- function(p, q){
  tol <- 1e-12
  radius_p <- sqrt(sum(p^2))
  radius_q <- sqrt(sum(q^2))
  stopifnot(abs(radius_p - radius_q) < tol)
  p <- p / radius_p
  q <- q / radius_q

  # Following https://github.com/JuliaManifolds/Manifolds.jl/blob/9cdc063df740dd0579f12469e3e663935de2df0e/src/manifolds/Sphere.jl#L296-L309
  cosAngle <- drop(pmin(pmax(t(p) %*% q, -1), 1))
  res <- if(abs(cosAngle + 1) < tol){ # cosAngle \approx -1
    res <- matrix(0, nrow = length(p))
    if(abs(p[1] - 1) < tol && length(p) > 1){
      res[2] <- 1
    }else{
      res[1] <- 1
    }
    res <- res - drop(t(p) %*% res) * p
    res * pi / sqrt(sum(res^2))
  }else if(abs(cosAngle - 1) < tol){ # cosAngle \approx 1
    q
  }else{
    angle <- acos(cosAngle)
    (q - cosAngle * p) * angle / sin(angle)
  }
  project_sphere_tangent(res, p)
}


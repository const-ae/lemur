

make_synthetic_data <- function(n_genes = 30, n_cells = 500, n_centers = 4, n_lat = 2, treatment_effect = 0.1){
  n_lat <- min(n_lat, n_genes)
  centers <- duplicate_cols(randn(n_lat, n_centers, sd = 2), ceiling(n_cells / n_centers))[,seq_len(n_cells)]
  true_Z <- centers + rnorm(n_cells * n_lat, sd = 0.1)
  stopifnot(n_centers <= length(LETTERS))
  cell_type <- rep(LETTERS[seq_len(n_centers)], length.out = n_cells)
  condition <- sample(letters[1:3], n_cells, replace = TRUE)
  design_matrix <- model.matrix(~ condition)
  plane <- qr.Q(qr(randn(n_genes, n_genes)))[seq_len(n_genes), seq_len(n_lat)]

  true_P <- nullspace(plane) %*% cbind(intercept = randn(n_genes - n_lat, 1, sd = 0),
                                        beta1 = randn(n_genes - n_lat, 1, sd = treatment_effect),
                                        beta2 = randn(n_genes - n_lat, 1, sd = treatment_effect))
  true_P <- true_P + rnorm(prod(dim(true_P)), sd = 1e-8)
  true_Beta <- randn(n_genes, ncol(design_matrix), sd = 0.1)

  dir <- scale(randn(n_lat, 1), center = FALSE)
  Y <- true_Beta %*% t(design_matrix) + plane %*% true_Z + true_P %*% (t(design_matrix) * duplicate_rows(t(dir) %*% true_Z, 3))

  linear_effect <- LinearEmbeddingMatrix(design_matrix, true_Beta)
  linear_embedding <- LinearEmbeddingMatrix(t(true_Z), plane)
  interaction_embedding <- LinearEmbeddingMatrix(t(t(design_matrix) * duplicate_rows(t(dir) %*% true_Z, 3)), true_P)

  colnames(Y) <- paste0("cell_", seq_len(n_cells))
  rownames(Y) <- paste0("gene_", seq_len(n_genes))
  sce <- SingleCellExperiment(list(logcounts = Y), colData = data.frame(condition = condition, cell_type = cell_type),
                       reducedDims = list(linear_effect = linear_effect, linear_embedding = linear_embedding, interaction_embedding = interaction_embedding))

}


make_vectors <- function(n_genes, n_obs, sd = 0.1){
  x <- randn(n_genes, n_obs)
  x <- apply(x, 2, \(.x) .x / sqrt(sum(.x^2)))
  z <- randn(n_genes, 1)
  z <- apply(z, 2, \(.z) .z / sqrt(sum(.z^2)))
  bp <- diag(nrow=n_genes)
  pert <- project_rotation_tangent(randn(n_genes, n_genes, sd = sd), bp)
  y <- rotation_map(pert, bp) %*% x
  list(x=x, y=y, z=z, bp=bp, pert=pert)
}

make_vectors2 <- function(n_genes, n_obs, sd = 0.1){
  x <- randn(n_genes, n_obs)
  z <- randn(n_genes, 1)
  bp <- diag(nrow=n_genes)
  pert <- project_spd_tangent(randn(n_genes, n_genes, sd = sd), bp)
  y <- spd_map(pert, bp) %*% x
  list(x=x, y=y, z=z, bp=bp, pert=pert)
}

principal_angle <- function(A, B){
  acos(pmax(0, pmin(1, svd(t(qr.Q(qr(A))) %*% qr.Q(qr(B)))$d))) / pi * 180
}

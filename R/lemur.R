

#' Main function to fit the latent embedding multivariate regression (LEMUR) model
#'
#' @param data a matrix with observations in the columns and features in the rows.
#'   Or a `SummarizedExperiment` / `SingleCellExperiment` object
#' @param design a formula referring to global objects or column in the `colData` of `data`
#'   and `col_data` argument
#' @param col_data an optional data frame with `ncol(data)` rows.
#' @param n_embedding the dimension of the $k$-plane that is rotated through space.
#' @param linear_coefficient_estimator specify which estimator is used to center the conditions.
#'   `"linear"` runs simple regression it works well in many circumstances but can produce poor
#'   results if the composition of the cell types changes between conditions (e.g., one cell type
#'   disappears). `"cluster_median"` works similar as `"linear"` but is robust against compositional
#'   changes. `"zero"` skips the centering step which is also robust against compositional changes.
#'   However, expression changes affecting all cells equally are not regressed out.
#' @param use_assay if `data` is a `SummarizedExperiment` / `SingleCellExperiment` object,
#'   which assay should be used.
#' @param test_fraction the fraction of cells that are split of before the model fit to keep an
#'   independent set of test observations. Default: 20% (`0.2`).
#' @param ... additional parameters that are passed on to the internal function `lemur_impl`.
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#'
#' @return an object of class `lemur_fit` which extends [`SingleCellExperiment`]. Accordingly,
#'   all functions that work for `sce`'s also work for `lemur_fit`'s. In addition, we
#'   give easy access to the fitted values using the dollar notation (e.g., `fit$embedding`).
#'   For details see the [lemur_fit-class] help page.
#'
#' @references
#'   * Ahlmann-Eltze, C. & Huber, W. (2023). Analysis of multi-condition single-cell data with latent
#'   embedding multivariate regression. bioRxiv [https://doi.org/10.1101/2023.03.06.531268](https://doi.org/10.1101/2023.03.06.531268)
#'
#' @seealso [`align_by_grouping`], [`align_harmony`], [`test_de`], [`find_de_neighborhoods`]
#'
#' @export
lemur <- function(data, design = ~ 1, col_data = NULL,
                  n_embedding = 15,
                  linear_coefficient_estimator = c("linear", "cluster_median", "zero"),
                  use_assay = "logcounts",
                  test_fraction = 0.2,
                  ...,
                  verbose = TRUE){

  data_mat <- handle_data_parameter(data, on_disk = FALSE, assay = use_assay)
  col_data <- glmGamPoi:::get_col_data(data, col_data)
  des <- handle_design_parameter(design, data, col_data)
  al_des <- des

  if(! is(data, "SummarizedExperiment")){
    data <- SingleCellExperiment(assays = setNames(list(data_mat), use_assay), colData = col_data)
  }

  # Avoid problems when splitting into test and training data
  col_data <- convert_dataframe_cols_chr_to_fct(col_data)

  # Create indicator vector which cells are used for training and which for testing
  is_test_data <- rep(FALSE, ncol(data))
  if(test_fraction < 0 || test_fraction >= 1){
    stop("'test_fraction' must be at least 0 and smaller than 1.")
  }else{
    if(verbose) message("Storing ", round(test_fraction, 2) * 100, "% of the data (", round(ncol(data) * test_fraction), " cells)",
                        " as test data.")
    is_test_data[sample.int(ncol(data), size = round(ncol(data) * test_fraction), replace = FALSE)] <- TRUE
  }

  design_matrix <- des$design_matrix[!is_test_data,,drop=FALSE]
  res <- lemur_impl(data_mat[,!is_test_data,drop=FALSE], design_matrix, n_embedding = n_embedding,
                    linear_coefficient_estimator = linear_coefficient_estimator, verbose = verbose, ...)
  alignment_design <- if(matrix_equals(design_matrix, res$alignment_design_matrix)){
    des$design_formula
  }else{
    NULL
  }
  embedding <- matrix(NA, nrow = res$n_embedding, ncol = ncol(data))
  embedding[,!is_test_data] <- res$embedding
  embedding[,is_test_data] <- project_on_lemur_fit_impl(Y = data_mat[,is_test_data,drop=FALSE], design_matrix = des$design_matrix[is_test_data,,drop=FALSE],
                                                      alignment_design_matrix = al_des$design_matrix[is_test_data,,drop=FALSE],
                                                      coefficients = res$coefficients, linear_coefficients = res$linear_coefficients,
                                                      alignment_coefficients = res$alignment_coefficients, base_point = res$base_point)

  lemur_fit(data, col_data = col_data,
            row_data = if(is(data, "SummarizedExperiment")) rowData(data) else NULL,
            n_embedding = res$n_embedding,
            design = des$design_formula,
            design_matrix = des$design_matrix,
            linear_coefficients = res$linear_coefficients,
            base_point = res$base_point,
            coefficients = res$coefficients,
            embedding = embedding,
            alignment_coefficients = res$alignment_coefficients,
            alignment_design = alignment_design,
            alignment_design_matrix = al_des$design_matrix,
            use_assay = use_assay,
            is_test_data = is_test_data)
}


lemur_impl <- function(Y, design_matrix,
                       n_embedding = 15,
                       base_point = c("global_embedding", "mean"),
                       linear_coefficient_estimator = c("linear", "cluster_median", "zero"),
                       linear_coefficients = NULL,
                       coefficients = NULL,
                       embedding = NULL,
                       alignment_coefficients = NULL,
                       alignment_design_matrix = NULL,
                       n_iter = 10, tol = 1e-8,
                       reshuffling_fraction = 0,
                       verbose = TRUE){
  alignment_coef_fixed_but_embedding_fitted <- ! is.null(alignment_coefficients) && is.null(embedding)
  linear_coefficient_estimator <- match.arg(linear_coefficient_estimator)

  # Set reduced dimensions
  stopifnot(n_embedding >= 0)
  n_ambient_eff <- nrow(Y)
  n_embedding <- min(n_embedding, nrow(Y), ncol(Y))
  linear_coef_fixed <-  ! is.null(linear_coefficients)
  diffemb_coef_fixed <- ! is.null(coefficients)
  embedding_fixed <- ! is.null(embedding)
  alignment_coef_fixed <- ! is.null(alignment_coefficients)
  if(is.null(alignment_design_matrix)){
    alignment_design_matrix <- design_matrix
  }


  # if(reshuffling_fraction != 0){
  #   stopifnot(reshuffling_fraction > 0 && reshuffling_fraction <= 1)
  #   warning("Reshuffling elements from the design matrix is a beta feature, that was ",
  #           "designed to regularize the differential inference for unmatched populations.\n",
  #           "On the down-side a large 'reshuffle_fraction' can adversely affect the accuracy of the inference.",
  #           "Please use with care.", call. = FALSE)
  #   sel <- sample.int(nrow(design_matrix), size = round(nrow(design_matrix) * reshuffling_fraction), replace = FALSE)
  #   shuf_sel <- sample(sel)
  #   design_matrix[sel,] <- design_matrix[shuf_sel,]
  #   alignment_design_matrix[sel,] <- alignment_design_matrix[shuf_sel,]
  # }

  # Initialize values
  if(linear_coef_fixed){
    if(length(linear_coefficients) == 1){
      linear_coefficients <- matrix(linear_coefficients, nrow = n_ambient_eff, ncol = ncol(design_matrix))
    }
    stopifnot(nrow(linear_coefficients) == n_ambient_eff & ncol(linear_coefficients))

  }else{
    if(verbose) message("Regress out global effects using ", linear_coefficient_estimator, " method.")
    linear_coefficients <- estimate_linear_coefficient(Y = Y, design_matrix = design_matrix, method = linear_coefficient_estimator)
  }
  if(linear_coefficient_estimator == "zero"){
    Y_clean <- Y
  }else{
    Y_clean <- Y - linear_coefficients %*% t(design_matrix)
  }
  if(!is.matrix(base_point)){
    if(verbose) message("Find base point for differential embedding")
    base_point <- find_base_point(Y_clean, base_point, n_embedding = n_embedding)
  }

  initial_error <- sum(Y_clean^2)
  if(verbose) message("Fit differential embedding model")
  if(verbose) message("Initial error: ", sprintf("%.3g", initial_error))
  if(! diffemb_coef_fixed){
    coefficients <- array(0, dim = c(n_ambient_eff, n_embedding, ncol(design_matrix)))
  }
  if(! embedding_fixed){
    embedding <- project_data_on_diffemb(Y_clean, design = design_matrix,
                                         coefficients = coefficients, base_point = base_point)
  }

  if(! diffemb_coef_fixed){
    if(verbose) message("---Fit Grassmann linear model")
    coefficients <- grassmann_lm(Y_clean, design = design_matrix, base_point = base_point)
  }
  if(! embedding_fixed){
    embedding <- project_data_on_diffemb(Y_clean, design = design_matrix,
                                                 coefficients = coefficients, base_point = base_point)
  }
  if(verbose){
    residuals <- Y - project_diffemb_into_data_space(embedding, design = design_matrix, coefficients = coefficients, base_point = base_point) - linear_coefficients %*% t(design_matrix)
    error <- sum(residuals^2)
    message("Final error: ", sprintf("%.3g", error))
  }


  if(alignment_coef_fixed_but_embedding_fitted){
    # Rotate the embedding if it wasn't provided
    stop("Fixing 'alignment_coefficients' without fixing 'embedding' is not implemented")
  }else{
    alignment_coefficients <- array(0, c(n_embedding, n_embedding, ncol(alignment_design_matrix)))
  }

  # Make sure that axes are ordered by variance
  if(prod(dim(embedding)) > 0 && all(!is.na(embedding))){
    svd_emb <- svd(embedding)
    rot <- svd_emb$u
    base_point <- base_point %*% rot
    for(idx in seq_len(dim(coefficients)[3])) {
      coefficients[,,idx] <- coefficients[,,idx] %*% rot
    }
    embedding <- t(svd_emb$v) * svd_emb$d
  }

  list(n_embedding = n_embedding,
       design_matrix = design_matrix, data = Y,
       linear_coefficients = linear_coefficients,
       base_point = base_point,
       coefficients = coefficients,
       embedding = embedding,
       alignment_coefficients = alignment_coefficients)
}


find_base_point <- function(Y_clean, base_point, n_embedding){
  n_genes <- nrow(Y_clean)
  if(is.matrix(base_point)){
    stopifnot(nrow(base_point) == n_genes)
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
      pca(Y_clean, n = n_embedding, center = FALSE)$coordsystem
    }else if(base_point_meth == "mean"){
      stop("'base_point = \"mean\"' is not implemented. Please use 'global_embedding'.")
    }
  }
}


project_diffemb_into_data_space <- function(embedding, design, coefficients, base_point){
  n_genes <- nrow(base_point)
  res <- matrix(NA, nrow = n_genes, ncol = ncol(embedding))
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  for(gr in unique(mm_groups)){
    covars <- design[which(mm_groups == gr)[1], ]
    res[,mm_groups == gr] <- grassmann_map(sum_tangent_vectors(coefficients, covars), base_point) %*% embedding[,mm_groups == gr,drop=FALSE]
  }
  res
}

project_data_on_diffemb <- function(Y_clean, design, coefficients, base_point){
  n_emb <- ncol(base_point)
  res <- matrix(NA, nrow = n_emb, ncol = ncol(Y_clean))
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  for(gr in unique(mm_groups)){
    covars <- design[which(mm_groups == gr)[1], ]
    res[,mm_groups == gr] <- t(grassmann_map(sum_tangent_vectors(coefficients, covars), base_point)) %*% Y_clean[,mm_groups == gr,drop=FALSE]
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








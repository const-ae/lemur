#' Enforce additional alignment of cell clusters beyond the direct differential embedding
#'
#' @param fit a `lemur_fit` object
#' @param grouping argument specific for `align_by_grouping`. Either a vector which assigns
#'   each cell to one group or a matrix with `ncol(fit)` columns where the rows are a soft-assignment
#'   to a cluster (i.e., columns sum to `1`). `NA`'s are allowed.
#' @param design a specification of the design (matrix or formula) that is used
#'   for the transformation. Default: `fit$design_matrix`
#' @param ridge_penalty specification how much the flexibility of the transformation
#'   should be regularized. Default: `0.01`
#' @param max_iter argument specific for `align_harmony`. The number of iterations. Default: `10`
#' @param preserve_position_of_NAs argument specific for `align_by_grouping`.
#'   Boolean flag to decide if `NA`s in the `grouping` mean that these cells should stay where they are (if
#'   possible) or if they are free to move around. Default: `FALSE`
#' @param ... additional parameters that are passed on to relevant functions
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#'
#' @returns The `fit` object with the updated `fit$embedding` and `fit$alignment_coefficients`.
#'
#' @export
align_harmony <- function(fit, design = fit$alignment_design,
                          ridge_penalty = 0.01, max_iter = 10, ..., verbose = TRUE){
  if(verbose) message("Select cells that are considered close with 'harmony'")
  if(is.null(attr(design, "ignore_degeneracy"))){
    # It doesn't matter for harmony if the design is degenerate
    attr(design, "ignore_degeneracy") <- TRUE
  }
  design_matrix <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))$design_matrix
  act_design_matrix <- design_matrix[!fit$is_test_data,,drop=FALSE]


  mm_groups <- get_groups(act_design_matrix, n_groups = ncol(act_design_matrix) * 10)
  if(! requireNamespace("harmony", quietly = TRUE)){
    stop("'harmony' is not installed. Please install it from CRAN.")
  }
  training_fit <- fit$training_data
  # Ignore best practice and call private methods from harmony
  harm_obj <- harmony_init(training_fit$embedding, act_design_matrix, ..., verbose = verbose)
  for(idx in seq_len(max_iter)){
    harm_obj <- harmony_max_div_clustering(harm_obj)

    alignment <- align_impl(training_fit$embedding, harm_obj$R, act_design_matrix, ridge_penalty = ridge_penalty)

    harm_obj$Z_corr <- alignment$embedding
    harm_obj$Z_cos <- t(t(harm_obj$Z_corr) / sqrt(colSums(harm_obj$Z_corr^2)))
    if(harm_obj$check_convergence(1)){
      if(verbose) message("Converged")
      break
    }
  }

  correct_fit(fit, alignment$alignment_coefficients, design)
}

#' @rdname align_harmony
#' @export
align_by_grouping <- function(fit, grouping, design = fit$alignment_design,
                              ridge_penalty = 0.01, preserve_position_of_NAs = FALSE, verbose = TRUE){
  if(verbose) message("Received sets of cells that are considered close")

  if(is.list(grouping)){
    stop("'grouping' must be a vector/factor with distinct elements for each group or ",
         "a matrix with `ncol(fit)` columns and one row per group.")
  }else if(! is.matrix(grouping)){
    stopifnot(length(grouping) == ncol(fit))
    grouping <- grouping[! fit$is_test_data]
  }else{
    stopifnot(ncol(grouping) == ncol(fit))
    grouping <- grouping[,! fit$is_test_data,drop=FALSE]
  }

  design_matrix <- handle_design_parameter(design, fit, glmGamPoi:::get_col_data(fit, NULL))$design_matrix
  act_design_matrix <- design_matrix[!fit$is_test_data,,drop=FALSE]

  res <- align_impl(fit$training_data$embedding, grouping, act_design_matrix,
                    ridge_penalty = ridge_penalty, calculate_new_embedding = FALSE)

  correct_fit(fit, res$alignment_coefficients, design)
}


#' Align the points according to some grouping
#'
#' @returns A list with the new embedding and the coefficients
#'
#' @keywords internal
align_impl <- function(embedding, grouping, design_matrix, ridge_penalty = 0.01,
                       preserve_position_of_NAs = FALSE, calculate_new_embedding = TRUE){
  if(! is.matrix(grouping)){
    grouping_matrix <- one_hot_encoding(grouping)
  }else{
    stopifnot(ncol(grouping) == ncol(embedding))
    stopifnot(all(grouping >= 0, na.rm = TRUE))
    # Make sure the entries sum to 1 (and don't touch them if the column is all zero)
    col_sums <- MatrixGenerics::colSums2(grouping)
    col_sums[col_sums == 0] <- 1
    grouping_matrix <- t(t(grouping) / col_sums)
  }

  # NA's are converted to zero columns ensuring that `diff %*% grouping_matrix = 0`
  grouping_matrix[,MatrixGenerics::colAnyNAs(grouping_matrix)] <- 0
  if(! preserve_position_of_NAs){
    all_zero_col <- MatrixGenerics::colSums2(grouping_matrix) == 0
    grouping_matrix <- grouping_matrix[,! all_zero_col,drop=FALSE]
    embedding <- embedding[,! all_zero_col,drop=FALSE]
    design_matrix <- design_matrix[! all_zero_col,,drop=FALSE]
  }

  stopifnot(ncol(embedding) == ncol(grouping_matrix))
  stopifnot(ncol(embedding) == nrow(design_matrix))

  n_groups <- nrow(grouping_matrix)
  n_emb <- nrow(embedding)
  K <- ncol(design_matrix)

  conditions <- get_groups(design_matrix, ncol(design_matrix) * 10)
  conds <- unique(conditions)

  # Calculate mean per cell_type+condition
  cond_ct_means <- lapply(conds, \(co){
    t(mply_dbl(seq_len(n_groups), \(idx){
      MatrixGenerics::rowWeightedMeans(embedding, w = grouping_matrix[idx,], cols = conditions == co)
    }, ncol = n_emb))
  })

  # Calculate target as mean of `cond_ct_means` per cell_type
  target <- matrix(NA, nrow = n_emb, ncol = n_groups)
  for(idx in seq_len(n_groups)){
    target[,idx] <- colMeans(mply_dbl(conds, \(co) cond_ct_means[[co]][,idx], ncol = n_emb), na.rm = TRUE)
  }

  # Shift all cells by `ctc_mean - ct_target` (`new_pos`)
  new_pos <- embedding
  for(co in conds){
    diff <- target - cond_ct_means[[co]]
    new_pos[,conditions == co] <- new_pos[,conditions == co] + diff %zero_dom_mat_mult% grouping_matrix[,conditions == co]
  }

  # Approximate shift by regressing `new_pos ~ S(x) * orig_pos + offset(x)`
  interact_design_matrix <- duplicate_cols(design_matrix, each = n_emb + 1) * duplicate_cols(t(rbind(1, embedding)), times = ncol(design_matrix))
  # interact_design_matrix <- duplicate_cols(design_matrix, each = n_emb) * duplicate_cols(t(embedding), times = ncol(design_matrix))
  alignment_coefs <- ridge_regression(new_pos - embedding, interact_design_matrix, ridge_penalty = ridge_penalty)
  alignment_coefs <- array(alignment_coefs, dim = c(n_emb, n_emb + 1, ncol(design_matrix)))

  new_emb <-if(calculate_new_embedding){
    apply_linear_transformation(embedding, alignment_coefs, design_matrix)
  }else{
    NULL
  }
  list(alignment_coefficients = alignment_coefs, embedding = new_emb)
}

#' Take a vector and convert it to a one-hot encoded matrix
#'
#' @returns A matrix with `length(unique(groups))` rows and `length(groups)` columns.
#'
#' @keywords internal
one_hot_encoding <- function(groups){
  if(is.factor(groups)){
    levels <- levels(groups)
  }else{
    levels <- unique(groups)
  }

  res <- matrix(0, nrow = length(levels), ncol = length(groups), dimnames = list(levels, names(groups)))
  for(i in seq_along(levels)){
    if(is.na(levels[i])){
      # Do nothing
    }else{
      res[i, groups == levels[i]] <- 1
    }
  }
  res
}

forward_linear_transformation <- function(alignment_coefficients, design_vector){
  n_emb <- dim(alignment_coefficients)[1]
  if(n_emb == 0){
    sum_tangent_vectors(alignment_coefficients, design_vector)
  }else{
    cbind(0, diag(nrow = n_emb)) + sum_tangent_vectors(alignment_coefficients, design_vector)
  }
}

reverse_linear_transformation <- function(alignment_coefficients, design_vector){
  n_embedding <- dim(alignment_coefficients)[1]
  if(n_embedding == 0){
    matrix(nrow = 0, ncol = 0)
  }else{
    solve(diag(nrow = n_embedding) + sum_tangent_vectors(alignment_coefficients[,-1,,drop=FALSE], design_vector))
  }
}

apply_linear_transformation <- function(A, alignment_coefs, design){
  mm_groups <- get_groups(design, n_groups = ncol(design) * 10)
  groups <- unique(mm_groups)
  for(gr in groups){
    A[,mm_groups == gr] <- forward_linear_transformation(alignment_coefs,  design[which(mm_groups == gr)[1],])  %*% rbind(1, A[,mm_groups == gr,drop=FALSE])
  }
  A
}

correct_fit <- function(fit, alignment_coefs, design){
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  if(! all(fit$alignment_coefficients == 0)) stop("Can only apply alignment once")

  metadata(fit)[["alignment_coefficients"]] <- alignment_coefs
  des <- handle_design_parameter(design, fit, fit$colData)
  metadata(fit)[["alignment_design"]] <- des$design_formula
  metadata(fit)[["alignment_design_matrix"]] <- des$design_matrix
  reducedDim(fit, "embedding") <- t(apply_linear_transformation(fit$embedding, alignment_coefs, des$design_matrix))
  fit
}

handle_ridge_penalty_parameter <- function(ridge_penalty){
  if(any(names(ridge_penalty) %in% c("rotation", "stretching"))){
    stop("The alignment function has changed and the rotation and stretching specification is now defunct")
  }
  ridge_penalty
}

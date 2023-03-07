
#' Differential expression for each cell (or position in the latent emebedding)
#'
#' @param fit the result of calling [`lemur()`]
#' @param contrast Specification of the contrast: a call to `cond()` specifying a full observation
#'    (e.g. `cond(treatment = "A", sex = "male") - cond(treatment = "C", sex = "male")` to
#'    compare treatment A vs C for male observations). Unspecified factors default to the reference level.
#' @param embedding matrix of size `n_embedding` \eqn{\times} `n` that specifies where in the latent space
#'   the differential expression is tested. It defaults to the position of all cells from the original fit.
#' @param consider specify which part of the model are considered for the differential expression test.
#'
#' @return if `is.null(embedding)` the `fit` object with a new assay called `"DE"`. Otherwise
#'  return the matrix with the differential expression values
#'
#' @export
test_de <- function(fit,
                    contrast,
                    embedding = NULL,
                    consider = c("embedding+linear", "embedding", "linear")){
  if(is.null(embedding)){
    embedding <- fit$embedding
    use_provided_diff_emb <- FALSE
  }else{
    use_provided_diff_emb <- TRUE
  }

  consider <- match.arg(consider)
  with_lm <- consider == "embedding+linear" || consider == "linear"
  with_emb <- consider == "embedding+linear" || consider == "embedding"


  cntrst <- parse_contrast({{contrast}}, formula = fit$design)
  al_cntrst <- parse_contrast({{contrast}}, formula = fit$alignment_design)
  diff <- evaluate_contrast_tree(cntrst, al_cntrst, \(x, y){
    predict(fit, newdesign = x, alignment_design_matrix = y, embedding = embedding, with_linear_model = with_lm, with_differential_embedding = with_emb)
  })

  colnames(diff) <- colnames(embedding)
  rownames(diff) <- rownames(fit)
  if(use_provided_diff_emb){
    diff
  }else{
    assay(fit, "DE") <- diff
    fit
  }
}


#' Differential embedding for each condition
#'
#' @param fit the result of [`differential_embedding`]
#' @param contrast Specification of the contrast: a call to `cond()` specifying a full observation
#'    (e.g. `cond(treatment = "A", sex = "male") - cond(treatment = "C", sex = "male")` to
#'    compare treatment A vs C for male observations). Unspecified factors default to the reference level.
#' @param reduced_design an alternative specification of the null hypothesis.
#' @param consider specify which part of the model are considered for the differential expression test.
#' @param variance_est How or if the variance should be estimated. `'analytical'` is only compatible with `consider = "linear"`. `'resampling'` is the most flexible (to adapt the number
#'   of resampling iterations, set `n_resampling_iter`. Default: `100`)
#'
#' @return a data.frame
#'
#' @export
test_global <- function(fit,
                        contrast,
                        reduced_design = NULL,
                        consider = c("embedding+linear", "embedding", "linear"),
                        variance_est = c("analytical", "resampling", "none"), verbose = TRUE,
                        ...){


  variance_est <- match.arg(variance_est)
  full_design <- fit$design_matrix
  consider <- match.arg(consider)
  with_lm <- consider == "embedding+linear" || consider == "linear"
  with_emb <- consider == "embedding+linear" || consider == "embedding"

  # Implement with a likelihood ratio test (see glmGamPoi)
  if(is.null(reduced_design) == missing(contrast)){
    stop("Please provide either an alternative design (formula or matrix) or a contrast.")
  }else if(! missing(contrast)){
    cntrst <- parse_contrast({{contrast}}, formula = fit$design)
    if(inherits(cntrst, "contrast_relation") && cntrst$relation == "minus" &&
        inherits(cntrst$lhs, "model_vec") && inherits(cntrst$rhs, "model_vec")){
      lfc_diffemb <- grassmann_log(grassmann_map(sum_tangent_vectors(fit$coefficients, c(cntrst$lhs)), fit$base_point),
                                   grassmann_map(sum_tangent_vectors(fit$coefficients, c(cntrst$rhs)), fit$base_point))
      cntrst <- cntrst$lhs - cntrst$rhs
    }else{
      cntrst <- evaluate_contrast_tree(cntrst, cntrst, \(x, .) x) # Collapse tree
      lfc_diffemb <- sum_tangent_vectors(fit$coefficients, c(cntrst))
    }

    cntrst <- as.matrix(cntrst)
    if(nrow(cntrst) != ncol(full_design)){
      stop("The length of the contrast vector does not match the number of coefficients in the model (",
           ncol(full_design), ")\n", glmGamPoi:::format_matrix(cntrst))
    }
    # The modifying matrix of reduced_design has ncol(design_matrix) - 1 columns and rank.
    # The columns are all orthogonal to cntrst.
    # see: https://scicomp.stackexchange.com/a/27835/36204
    # Think about this as a rotation of of the design matrix. The QR decomposition approach
    # has the added benefit that the new columns are all orthogonal to each other, which
    # isn't necessary, but makes fitting more robust
    # The following is a simplified version of edgeR's glmLRT (line 159 in glmfit.R)
    qrc <- qr(cntrst)
    rot <- qr.Q(qrc, complete = TRUE)[,-1,drop=FALSE]
    reduced_design_mat <- full_design %*% rot

    lfc_linear_model <- fit$linear_coefficients %*% cntrst
  }else{
    reduced_design_mat <- handle_design_parameter(reduced_design, fit, fit$colData)$design_matrix
    if(ncol(reduced_design_mat) >= ncol(full_design)){
      stop("The reduced model is as complex (or even more complex) than ",
           "the 'fit' model. The 'reduced_design' should contain fewer terms ",
           "the original 'design'.")
    }
    rot <- matrix(lm.fit(full_design, reduced_design_mat)$coefficients, ncol = ncol(reduced_design_mat))
    if(any(abs(reduced_design_mat - full_design %*% rot) > 1e-10)){
      warning("Although the 'reduced_design' matrix has fewer columns than ",
              "'fit$design_matrix', it appears that the 'reduced_design' is not ",
              "nested in the 'fit$design_matrix'. Accordingly, the results of the ",
              "statistical test will be unreliable.")
    }
  }



  if(variance_est == "analytical"){
    if(with_emb){
      stop("Analytical differential embedding test is not implemented. You can set 'consider=\"linear\"")
    }else{  # only linear test
      resid_full <- as.matrix(t(fit$ambient_coordsystem) %*% residuals(fit, with_linear_model = TRUE, with_differential_embedding = FALSE, with_ambient_pca = TRUE))
      resid_red <- as.matrix(t(fit$ambient_coordsystem) %*% get_residuals_for_alt_fit(fit, reduced_design_mat = reduced_design_mat, with_linear_model = TRUE, with_differential_embedding = FALSE))
      pval <- multivar_wilks_ftest(RSS_full = resid_full %*% t(resid_full),
                                   RSS_red = resid_red %*% t(resid_red),
                                   n_features = fit$n_ambient, full_design, reduced_design_mat)
    }
  }else if(variance_est == "resampling"){
    if("n_resampling_iter" %in% ...names()){
      n_resampling_iter <- list(...)[["n_resampling_iter"]]
    }else{
      n_resampling_iter <- 99
    }
    if(verbose) message("Estimating null distribution of deviance using ", n_resampling_iter, " iterations.")
    # Applying the Freedman-Lane (1983) permutation method of the residuals
    # Fit the full
    deviance_ref <- sum(residuals(fit, with_linear_model = with_lm, with_differential_embedding = with_emb)^2)
    # Fit the reduced model
    resid_red <- get_residuals_for_alt_fit(fit, reduced_design_mat = reduced_design_mat, with_linear_model = with_lm, with_differential_embedding = with_emb)
    predict_red <- assay(fit, "expr") - resid_red
    deviance_red <- sum(resid_red^2)
    deviance_delta_null <- vapply(seq_len(n_resampling_iter), \(iter){
      new_Y <- predict_red + resid_red[,sample.int(ncol(resid_red), replace = FALSE),drop=FALSE]
      deviance_ref_new <- sum(get_residuals_for_alt_fit(fit, Y = new_Y, reduced_design_mat = full_design, with_linear_model = with_lm, with_differential_embedding = with_emb)^2)
      deviance_red_new <- sum(get_residuals_for_alt_fit(fit, Y = new_Y, reduced_design_mat = reduced_design_mat, with_linear_model = with_lm, with_differential_embedding = with_emb)^2)
      deviance_red_new - deviance_ref_new
    }, FUN.VALUE = numeric(1L))
    pval <- (sum((deviance_red - deviance_ref) < deviance_delta_null) + 1) / (n_resampling_iter + 1)
  }else{ # variance_est == "none"
    pval <- NA
  }


  if(! missing(contrast)){
    data.frame(contrast = rlang::as_label(rlang::enquo(contrast)),
               pval = pval,
               delta_diffemb = I(list(lfc_diffemb)),
               delta_linear = I(list(lfc_linear_model)),
               angle_degrees = grassmann_angle_from_tangent(lfc_diffemb, normalized = TRUE))
  }else{
    data.frame(full_design = if(! is.null(fit$design)) rlang::as_label(fit$design) else rlang::as_label(fit$design_matrix),
               reduced_design = rlang::as_label(reduced_design),
               pval = pval)
  }
}


multivar_wilks_ftest <- function(RSS_full, RSS_red, n_features, design_matrix_full, design_matrix_red){
  # Following https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Multivariate-Linear-Models.pdf
  stopifnot(nrow(design_matrix_full) == nrow(design_matrix_red))
  k1 <- ncol(design_matrix_full)
  k2 <- ncol(design_matrix_red)

  lambdas <- Re(eigen((RSS_red - RSS_full) %*% solve(RSS_full))$values)
  wilks_lambda <- prod(1/(1 + lambdas))
  df <- k1 - k2
  r <- nrow(design_matrix_full) - k1 - 1 - (n_features - df + 1) / 2
  u <- (n_features * df - 2) / 4
  t <- ifelse(n_features^2 + df^2 - 5 <= 0, 0, sqrt(n_features^2 * df^2 - 4) / (n_features^2 + df^2 - 5))
  fstat <- (1 - wilks_lambda^(1/t)) / wilks_lambda^(1/t) * (r * t - 2 * u) / (n_features * df)
  pf(fstat, df1 = n_features * df, df2 = r * t - 2 * u, lower.tail = FALSE)
}


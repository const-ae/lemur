
#' Differential expression for each cell (or position in the latent emebedding)
#'
#' @param fit the result of [`differential_embedding`]
#' @param contrast Specification of the contrast. This can be either
#'    * a numeric vector whose length matches the number of coefficients (e.g., `c(1, 0, -1, 0)`
#'       to compare the coefficients A vs C in a model with four treatment options)
#'    * the unquoted names of coefficients to compare (e.g., `treatmentA - treatmentC`)
#'    * a call to `fact()` specifying a full observation (e.g. `fact(treatment = "A", sex = "male") - fact(treatment = "C", sex = "male")` to
#'       compare treatment A vs C for male observations). Unspecified factors default to the reference level.
#'    * an extension of the previous two options where instead of subtracting the
#'       coefficients, they are compared directly (e.g. `fact(treatment = "A", sex = "male") < fact(treatment = "C", sex = "male")` or
#'       `fact(treatment = "A", sex = "male") == treatmentC`). This is the recommended approach, because `map(V1 - V2) != map(V1) - map(V2)`.
#' @param alignment_contrast same as `contrast` but applied to the `alignment_design_matrix`. This is for advanced use cases where
#'   separate experimental designs are used in the multi-condition PCA and the alignment step. Defaults to the `contrast` argument.
#' @param diffemb_embedding matrix of size `n_embedding` \eqn{\times} `n` that specifies where in the latent space
#'   the differential expression is tested. It defaults to the position of all cells from the original fit.
#' @param consider specify which part of the model are considered for the differential expression test.
#' @param variance_est How or if the variance should be estimated
#' @param return Should the function return a matrix or a long table of the results
#'
#' @export
test_differential_expression <- function(fit,
                                        contrast,
                                        alignment_contrast = {{contrast}},
                                        diffemb_embedding = NULL,
                                        consider = c("embedding+linear", "embedding", "linear"),
                                        variance_est = c("bootstrap", "none"),
                                        return = c("table", "matrix")){
  variance_est <- match.arg(variance_est)
  return <- match.arg(return)
  if(variance_est == "bootstrap" && is.null(fit$bootstrap_samples)){
    stop("No bootstrap samples available. Please call 'estimate_variance()' before calling 'test_differential_expression()'.")
  }
  if(is.null(diffemb_embedding)){
    diffemb_embedding <- fit$diffemb_embedding
    use_provided_diff_emb <- FALSE
  }else{
    use_provided_diff_emb <- TRUE
  }

  if(return == "table"){
    feature_names <- if(is.null(rownames(fit))){
      rep(paste0("feature_", seq_len(nrow(fit))), times = ncol(diffemb_embedding))
    }else{
      rep(rownames(fit), times = ncol(diffemb_embedding))
    }
    obs_names <- if(is.null(colnames(diffemb_embedding))){
      rep(paste0("obs_", seq_len(ncol(diffemb_embedding))), each = nrow(fit))
    }else{
      rep(colnames(diffemb_embedding), each = nrow(fit))
    }
  }
  consider <- match.arg(consider)
  with_lm <- consider == "embedding+linear" || consider == "linear"
  with_emb <- consider == "embedding+linear" || consider == "embedding"


  cntrst <- parse_contrast({{contrast}}, coefficient_names = colnames(fit$design_matrix), formula = fit$design)
  al_cntrst <- parse_contrast({{alignment_contrast}}, coefficient_names = colnames(fit$alignment_design_matrix), formula = fit$alignment_design)
  diff <- if(inherits(cntrst, "contrast_relation") && inherits(al_cntrst, "contrast_relation")){
    predict(fit, newdesign = cntrst$lhs, alignment_design_matrix = al_cntrst$lhs, diffemb_embedding = diffemb_embedding, with_linear_model = with_lm, with_differential_embedding = with_emb) -
      predict(fit, newdesign = cntrst$rhs, alignment_design_matrix = al_cntrst$rhs, diffemb_embedding = diffemb_embedding, with_linear_model = with_lm, with_differential_embedding = with_emb)
  }else if(inherits(cntrst, "contrast_relation") != inherits(al_cntrst, "contrast_relation")){
    stop("Both 'contrast' and 'alignment_contrast' must contain an equality relation ('==') or neither.")
  }else{
    predict(fit, newdesign = cntrst, alignment_design_matrix = al_cntrst, diffemb_embedding = diffemb_embedding, with_linear_model = with_lm, with_differential_embedding = with_emb)
  }
  if(variance_est == "bootstrap"){
    # Welfords's online algorithm to calculate mean and sd of bootstrap estimates
    msq_init <- mean_init <- matrix(0, nrow = nrow(fit), ncol = ncol(diffemb_embedding))
    preds <- fold_left(list(mean = mean_init, msq = msq_init, iter = 1))(fit$bootstrap_samples, \(elem, accum){
      bs_diffemb_embedding <- if(use_provided_diff_emb){
        diffemb_embedding
      }else{
        elem$diffemb_embedding
      }
      diff <- if(inherits(cntrst, "contrast_relation")){
        linearCoefficients <- elem$linear_coefficients
        predict(elem, newdesign = cntrst$lhs, alignment_design_matrix = al_cntrst$lhs, diffemb_embedding = bs_diffemb_embedding, with_linear_model = with_lm,
                with_differential_embedding = with_emb) -
          predict(elem, newdesign = cntrst$rhs, alignment_design_matrix = al_cntrst$rhs, diffemb_embedding = bs_diffemb_embedding, with_linear_model = with_lm,
                  with_differential_embedding = with_emb)
      }else{
        predict(elem, newdesign = cntrst, alignment_design_matrix = al_cntrst, diffemb_embedding = bs_diffemb_embedding, with_linear_model = with_lm,
                with_differential_embedding = with_emb)
      }
      delta <- diff - accum$mean
      accum$mean <- accum$mean + delta / accum$iter
      delta2 <- diff - accum$mean
      accum$msq <- accum$msq + delta * delta2
      accum$iter <- accum$iter + 1
      accum
    })
    # The point estimate using the original data is better
    # than the mean of the bootstrap samples (although they are close anyways)
    # diff <- preds$mean
    sd <- sqrt(preds$msq / (length(fit$bootstrap_samples) - 1))

    if(return == "matrix"){
      colnames(sd) <- colnames(diff) <- colnames(diffemb_embedding)
      rownames(sd) <- rownames(diff) <- rownames(fit)

      list(diff = diff, sd = sd)
    }else{
      pval <- if(! inherits(cntrst, "contrast_relation") || cntrst$relation == "equal"){
        pmin(pnorm(diff / sd, lower.tail = TRUE), pnorm(diff / sd, lower.tail = FALSE)) * 2
      }else if(cntrst$relation == "less_than"){
        pnorm(diff / sd, lower.tail = FALSE)
      }else if(cntrst$relation == "greather_than"){
        pnorm(diff / sd, lower.tail = TRUE)
      }
      adj_pval <- p.adjust(pval, "BH")
      data.frame(feature = feature_names,
                 obs = obs_names,
                 pval = c(pval),
                 adj_pval = c(adj_pval),
                 diff = c(diff),
                 adj_diff = c(diff / sd),
                 sd = c(sd), stringsAsFactors = FALSE, row.names = NULL)
    }
  }else{ # variance_est == "none"
    if(return == "matrix"){
      colnames(diff) <- colnames(diffemb_embedding)
      rownames(diff) <- rownames(fit)
      diff
    }else{
      data.frame(feature = feature_names,
                 obs = obs_names,
                 diff = c(diff), stringsAsFactors = FALSE, row.names = NULL)
    }
  }
}


#' Differential embedding for each condition
#'
#' @param fit the result of [`differential_embedding`]
#' @param contrast Specification of the contrast. This can be either
#'    * a numeric vector whose length matches the number of coefficients (e.g., `c(1, 0, -1, 0)`
#'       to compare the coefficients A vs C in a model with four treatment options)
#'    * the unquoted names of coefficients to compare (e.g., `treatmentA - treatmentC`)
#'    * a call to `fact()` specifying a full observation (e.g. `fact(treatment = "A", sex = "male") - fact(treatment = "C", sex = "male")` to
#'       compare treatment A vs C for male observations). Unspecified factors default to the reference level.
#'    * an extension of the previous two options where instead of subtracting the
#'       coefficients, they are compared directly (e.g. `fact(treatment = "A", sex = "male") < fact(treatment = "C", sex = "male")` or
#'       `fact(treatment = "A", sex = "male") == treatmentC`). This is the recommended approach, because `map(V1 - V2) != map(V1) - map(V2)`.
#' @param reduced_design an alternative specification of the null hypothesis.
#' @param consider specify which part of the model are considered for the differential expression test.
#' @param variance_est How or if the variance should be estimated. `'bootstrap'` is only compatible with a constrast,
#'   `'analytical'` is only compatible with `consider = "linear"`. `'resampling'` is the most flexible (to adapt the number
#'   of resampling iterations, set `n_resampling_iter`. Default: `100`)
#'
#' @return a data.frame
#'
#' @export
test_differential_embedding <- function(fit,
                                        contrast,
                                        reduced_design = NULL,
                                        consider = c("embedding+linear", "embedding", "linear"),
                                        variance_est = c("bootstrap", "analytical", "resampling", "none"), verbose = TRUE,
                                        ...){


  variance_est <- match.arg(variance_est)
  if(variance_est == "bootstrap" && is.null(fit$bootstrap_samples)){
    stop("No bootstrap samples available. Please call 'estimate_variance()' before calling 'test_differential_expression()'.")
  }else if(variance_est == "bootstrap" && missing(contrast)){
    stop("Boostrap test is only compatible with 'contrast' argument. Not with a 'reduced_design'.")
  }
  full_design <- fit$design_matrix
  consider <- match.arg(consider)
  with_lm <- consider == "embedding+linear" || consider == "linear"
  with_emb <- consider == "embedding+linear" || consider == "embedding"

  # Implement with a likelihood ratio test (see glmGamPoi)
  if(is.null(reduced_design) == missing(contrast)){
    stop("Please provide either an alternative design (formula or matrix) or a contrast.")
  }else if(! missing(contrast)){
    cntrst <- parse_contrast({{contrast}}, coefficient_names = colnames(full_design), formula = fit$design)
    if(inherits(cntrst, "contrast_relation")){
      if(cntrst$relation != "equal"){
        stop("differential embedding test can only be two-sided.")
      }
      lfc_diffemb <- grassmann_log(grassmann_map(sum_tangent_vectors(fit$diffemb_coefficients, c(cntrst$lhs)), fit$diffemb_basepoint),
                                   grassmann_map(sum_tangent_vectors(fit$diffemb_coefficients, c(cntrst$rhs)), fit$diffemb_basepoint))
      cntrst <- cntrst$lhs - cntrst$rhs
    }else{
      lfc_diffemb <- sum_tangent_vectors(fit$diffemb_coefficients, c(cntrst))
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
  }else if(variance_est == "bootstrap"){
    # Get all fit values and check if they are different from zero
    vals <- matrix(nrow = 0L, ncol = length(fit$bootstrap_samples))
    n_ambient_eff <- min(fit$n_ambient, nrow(fit))
    if(with_lm){
      vals <- rbind(vals, t(mply_dbl(fit$bootstrap_samples, \(bs) c(bs$linear_coefficients %*% cntrst), ncol = n_ambient_eff)))
    }
    if(with_emb){
      vals <- rbind(vals, t(mply_dbl(fit$bootstrap_samples, \(bs) c(sum_tangent_vectors(bs$diffemb_coefficients, cntrst)), ncol = n_ambient_eff * fit$n_embedding)))
    }
    # Filter out zero variance obs
    # vals <- vals[matrixStats::rowSds(vals) > 1e-6,,drop=FALSE]

    # Do an adapted version of the Hotelling z-test (based on the Mahalanobis distance)
    # that can handle n_bootstrap < n_ambient * n_coef * (n_embedding + 1)
    # see https://stats.stackexchange.com/a/514628/130486
    # zstat <- drop(t(rowMeans(vals)) %*% corpcor::invcov.shrink(t(vals), verbose = FALSE) %*% rowMeans(vals))
    zstat <- drop(t(rowMeans(vals)) %*% diag(1/(matrixStats::rowVars(vals)+1e-5), nrow = nrow(vals)) %*% rowMeans(vals))
    pval <- pchisq(zstat, df = nrow(vals), lower.tail = FALSE)
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


test_differential_abundance <- function(fit,
                                        contrast,
                                        reduced_design = NULL,
                                        variance_est = c("none", "bootstrap")){

  stop("Not implemented")

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


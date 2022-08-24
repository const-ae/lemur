

test_differential_expression <- function(fit,
                                        contrast,
                                        reduced_design = NULL,
                                        diffemb_embedding = fit$diffemb_embedding,
                                        consider = c("embedding+linear", "embedding", "linear"),
                                        variance_est = c("bootstrap", "none"),
                                        return = c("table", "matrix")){
  variance_est <- match.arg(variance_est)
  return <- match.arg(return)
  if(variance_est == "bootstrap" && is.null(fit$bootstrap_samples)){
    stop("No bootstrap samples available. Please call 'estimate_variance()' before calling 'test_differential_expression()'.")
  }
  if(! is.null(reduced_design)){
    stop("'reduced_design' is not implemented for 'test_differential_expression'.")
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
  diff <- if(inherits(cntrst, "contrast_relation")){
    predict(fit, newdesign = cntrst$lhs, diffemb_embedding = diffemb_embedding, with_linear_model = with_lm, with_differential_embedding = with_emb) -
      predict(fit, newdesign = cntrst$rhs, diffemb_embedding = diffemb_embedding, with_linear_model = with_lm, with_differential_embedding = with_emb)
  }else{
    predict(fit, newdesign = cntrst, diffemb_embedding = diffemb_embedding, with_linear_model = with_lm, with_differential_embedding = with_emb)
  }
  if(variance_est == "bootstrap"){
    preds <- vapply(fit$bootstrap_samples, \(bs){
      test_differential_expression(bs, cntrst, reduced_design = reduced_design, diffemb_embedding = diffemb_embedding,
                                   consider = consider, variance_est = "none", return = "matrix")
    }, FUN.VALUE = array(0.0, c(nrow(fit), ncol(diffemb_embedding))))
    diff <- apply(preds, c(1,2), mean)
    sd <- apply(preds, c(1,2), sd)

    if(return == "matrix"){
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
      diff
    }else{
      data.frame(feature = feature_names,
                 obs = obs_names,
                 diff = c(diff), stringsAsFactors = FALSE, row.names = NULL)
    }
  }
}


test_differential_embedding <- function(fit,
                                        contrast,
                                        reduced_design = NULL,
                                        consider = c("embedding+linear", "embedding", "linear"),
                                        variance_est = c("bootstrap", "analytical", "none"), verbose = TRUE){


  variance_est <- match.arg(variance_est)
  full_design <- fit$design_matrix
  two_sided <- TRUE
  dir_less_than <- NA
  consider <- match.arg(consider)
  with_lm <- consider == "embedding+linear" || consider == "linear"
  with_emb <- consider == "embedding+linear" || consider == "embedding"

  # Implement with a likelihood ratio test (see glmGamPoi)
  if(is.null(reduced_design) == missing(contrast)){
    stop("Please provide either an alternative design (formula or matrix) or a contrast.")
  }else if(! missing(contrast)){
    cntrst <- parse_contrast({{contrast}}, coefficient_names = colnames(fit$design_matrix), formula = fit$design)
    if(inherits(cntrst, "contrast_relation")){
      two_sided <- cntrst$relation == "equal"
      direction_less_than <- if(two_sided) NA else cntrst$relation == "less_than"
      cntrst <- cntrst$rhs - cntrst$lhs
    }
    cntrst <- as.matrix(cntrst)
    if(nrow(cntrst) != ncol(fit$design_matrix)){
      stop("The length of the contrast vector does not match the number of coefficients in the model (",
           ncol(fit$design_matrix), ")\n", glmGamPoi:::format_matrix(cntrst))
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
    reduced_design_mat <- fit$design_matrix %*% rot
    lfc_diffemb <- sum_tangent_vectors(fit$diffemb_coefficients, c(cntrst))
    lfc_linear_model <- fit$linear_coefficients %*% cntrst
  }else{
    reduced_design_mat <- handle_design_parameter(reduced_design, fit, fit$colData)$model_matrix
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
    warning("Analytical calculation of p-values doesn't work yet.")
    fit_alt <- differential_embedding_impl(matrix(nrow = nrow(fit), ncol = 0), design = reduced_design_mat,
                                           n_ambient = fit$n_ambient, n_embedding = fit$n_embedding,
                                           alignment = fit$alignment_method, base_point = fit$diffemb_basepoint,
                                           amb_pca = list(coordsystem = fit$ambient_coordsystem,
                                                          embedding = t(fit$ambient_coordsystem) %*% (assay(fit, "expr") - fit$ambient_offset),
                                                          offset = fit$ambient_offset),
                                           verbose = verbose)
    dev_full <- sum((predict(fit, with_linear_model = with_lm, with_differential_embedding = with_emb) - assay(fit, "expr"))^2)
    dev_red <- sum((predict_impl(object = NULL, diffemb_embedding = fit_alt$diffemb_embedding,
                                 with_linear_model = with_lm, with_differential_embedding = with_emb,
                                 n_ambient = fit_alt$n_ambient, n_embedding = fit_alt$n_embedding,
                                 design_matrix = fit_alt$design_matrix, design = fit_alt$design,
                                 ambient_coordsystem = fit_alt$ambient_coordsystem, ambient_offset = fit_alt$ambient_offset,
                                 linear_coefficients = fit_alt$linear_coefficients, diffemb_coefficients = fit_alt$diffemb_coefficients,
                                 diffemb_basepoint = fit_alt$diffemb_basepoint, alignment_coefficients = fit_alt$alignment_coefficients) - assay(fit, "expr"))^2)
    log_likelihood_ratio <- dev_red - dev_full
    df_test <- (ncol(fit$design_matrix) - ncol(fit_alt$design_matrix)) * (ifelse(with_lm, nrow(fit$design_matrix), 0) +
                                                                            ifelse(with_emb, fit$n_ambient * fit$n_embedding, 0))
    if(! two_sided) warning("Only two-sided tests are implemented.")
    pval <- pchisq(log_likelihood_ratio, df = df_test, lower.tail = FALSE)
  }else if(variance_est == "bootstrap"){
    stop("I am not sure what I need to do here")
  }else{ # variance_est == "none"
    pval <- NA
  }


  if(! missing(contrast)){
    data.frame(contrast = rlang::as_label(rlang::enquo(contrast)),
               pval = pval,
               change_diffemb = I(list(lfc_diffemb)),
               change_linear = I(list(lfc_linear_model)),
               dev_full = dev_full, dev_red = dev_red, df = df_test)
  }else{
    data.frame(full_design = if(! is.null(fit$design)) rlang::as_label(fit$design) else paste0("design matrix with ", ncol(fit$design_matrix), " columns"),
               reduced_design = rlang::as_label(reduced_design),
               pval = pval, dev_full = dev_full, dev_red = dev_red, df = df_test)
  }
}


test_differential_abundance <- function(fit,
                                        contrast,
                                        reduced_design = NULL,
                                        variance_est = c("none", "bootstrap")){

  stop("Not implemented")

}




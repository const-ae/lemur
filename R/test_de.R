

test_differences <- function(fit,
                             contrast,
                             reduced_design = NULL,
                             pseudobulk_by = NULL,
                             mle = FALSE,
                             ...){

  stop("Not implemented")

}


test_differential_expression <- function(fit,
                                        contrast,
                                        reduced_design = NULL,
                                        consider = c("embedding+linear", "embedding", "linear"),
                                        pseudobulk_by = NULL,
                                        variance_est = c("bootstrap", "none"),
                                        return = c("table", "matrix")){
  variance_est <- match.arg(variance_est)
  return <- match.arg(return)
  if(variance_est == "bootstrap" && is.null(fit$bootstrap_samples)){
    stop("No bootstrap samples available. Please call 'estimate_variance()' before calling 'test_differential_expression()'.")
  }
  if(return == "table"){
    feature_names <- if(is.null(rownames(fit))){
      paste0("feature_", seq_len(nrow(fit)))
    }else{
      rep(rownames(fit), times = ncol(fit))
    }
    obs_names <- if(is.null(colnames(fit))){
      paste0("obs_", seq_len(nrow(fit)))
    }else{
      rep(colnames(fit), each = nrow(fit))
    }
  }
  consider <- match.arg(consider)
  with_lm <- consider == "embedding+linear" || consider == "linear"
  with_emb <- consider == "embedding+linear" || consider == "embedding"


  cntrst <- parse_contrast({{contrast}}, coefficient_names = colnames(fit$design_matrix), formula = fit$design)
  diff <- if(inherits(cntrst, "contrast_relation")){
    predict(fit, newdesign = cntrst$lhs, with_linear_model = with_lm, with_differential_embedding = with_emb) -
      predict(fit, newdesign = cntrst$rhs, with_linear_model = with_lm, with_differential_embedding = with_emb)
  }else{
    predict(fit, newdesign = cntrst, with_linear_model = with_lm, with_differential_embedding = with_emb)
  }
  if(variance_est == "bootstrap"){
    preds <- vapply(fit$bootstrap_samples, \(bs){
      test_differential_expression(bs, cntrst, consider = consider, variance_est = "none", return = "matrix")
    }, FUN.VALUE = array(0.0, dim(fit)))
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
                                        pseudobulk_by = NULL,
                                        variance_est = c("bootstrap", "analytical", "none")){

  stop("Not implemented")
  # Implement with a likelihood ratio test (see glmGamPoi)

}


test_differential_abundance <- function(fit,
                                        contrast,
                                        reduced_design = NULL,
                                        pseudobulk_by = NULL,
                                        variance_est = c("none", "bootstrap")){

  stop("Not implemented")

}






edger_fit <- function(counts, design, offset, col_data = NULL,
                      abundance.trend = TRUE, robust = TRUE){

  if(is.matrix(design)){
    design_matrix <- design
  }else{
    design_matrix <- model.matrix(design, data = col_data)
  }

  edger_y <- edgeR::DGEList(counts = counts)
  edger_y <- edgeR::scaleOffset(edger_y, offset = offset)
  edger_y <- edgeR::estimateDisp(edger_y, design_matrix)
  edgeR::glmQLFit(edger_y, design_matrix, abundance.trend = abundance.trend, robust = robust)
}


edger_test_de <- function(edger_fit, contrast, design = NULL){
  if(rlang::is_quosure(contrast)){
    cntrst <- parse_contrast({{contrast}}, design)
    cntrst <- evaluate_contrast_tree(cntrst, cntrst, \(x, .) x)
  }else{
    cntrst <- contrast
  }

  edger_fit <- edgeR::glmQLFTest(edger_fit, contrast = cntrst)
  edger_res <- edgeR::topTags(edger_fit, n = nrow(edger_fit$counts), sort.by = "none")$table
  data.frame(name = rownames(edger_res), pval = edger_res$PValue, adj_pval = edger_res$FDR,
             f_statistic = edger_res$F, df1 = edger_fit$df.test, df2 = edger_fit$df.total, lfc = edger_res$logFC)
}


limma_fit <- function(values, design, col_data){
  model_matrix <- model.matrix(design, data = col_data)
  if(! is_contrast_estimable(cntrst, model_matrix)){
    stop("The contrast is not estimable from the model_matrix")
  }

  suppressWarnings({
    # limma warns about missing values. Here we expect missing values though.
    lm_fit <- limma::lmFit(values, model_matrix)
  })
  lm_fit
}

limma_test_de <- function(lm_fit, contrast, design, values = NULL, shrink = TRUE, trend = TRUE, robust = TRUE){
  cntrst <- parse_contrast({{contrast}}, formula = design)
  cntrst <- matrix(evaluate_contrast_tree(cntrst, cntrst, \(x, .) x), ncol = 1)

  lm_fit <- limma::contrasts.fit(lm_fit, contrasts = cntrst)
  if(shrink){
    lm_fit <- tryCatch({
      limma::eBayes(lm_fit, trend = trend, robust = robust)
    }, error = function(err){
      limma::eBayes(lm_fit, trend = FALSE, robust = TRUE)
    })
  }else{
    lm_fit <- limma_eBayes_without_shrinkage(lm_fit)
  }
  tt <- limma::topTable(lm_fit, number = nrow(lm_fit$coefficients), adjust.method = "BH", sort.by = "none")
  if(! is.null(values)){
    for(row in which(MatrixGenerics::rowAnyNAs(values))){
      # limma can return misleading results if there missing values in the wrong
      # places (see https://support.bioconductor.org/p/9150300/)
      if(! is_contrast_estimable(cntrst, lm_fit$design[!is.na(values[row,]),,drop=FALSE])){
        tt[row, c("P.Value", "adj.P.Val", "t", "logFC")] <- NA_real_
      }
    }
  }
  tt
}



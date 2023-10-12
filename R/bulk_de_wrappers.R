

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
  cntrst <- parse_contrast({{contrast}}, design, simplify = TRUE)

  edger_fit <- edgeR::glmQLFTest(edger_fit, contrast = cntrst)
  edger_res <- edgeR::topTags(edger_fit, n = nrow(edger_fit$counts), sort.by = "none")$table
  data.frame(name = rownames(edger_res), pval = edger_res$PValue, adj_pval = edger_res$FDR,
             f_statistic = edger_res$F, df1 = edger_fit$df.test, df2 = edger_fit$df.total, lfc = edger_res$logFC)
}


limma_fit <- function(values, design, col_data = NULL){
  if(is.matrix(design)){
    design_matrix <- design
  }else{
    design_matrix <- model.matrix(design, data = col_data)
  }
  if(! is_contrast_estimable(cntrst, design_matrix)){
    stop("The contrast is not estimable from the design_matrix")
  }

  suppressWarnings({
    # limma warns about missing values. Here we expect missing values though.
    lm_fit <- limma::lmFit(values, design_matrix)
  })
  lm_fit
}

limma_test_de <- function(lm_fit, contrast, design, values = NULL, shrink = TRUE, trend = TRUE, robust = TRUE){
  cntrst <- matrix(parse_contrast({{contrast}}, formula = design, simplify = TRUE), ncol = 1)
  cntrst <- evaluate_contrast_tree(cntrst, cntrst, \(x, .) x)

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
  data.frame(name = rownames(tt), pval = tt$P.Value, adj_pval = tt$adj.P.Val, t_statistic = tt$t, lfc = tt$logFC)
}



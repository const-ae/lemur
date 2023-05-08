
#' Project new data onto the latent spaces of an existing lemur fit
#'
#' @param fit an `lemur_fit` object
#' @param data a matrix with observations in the columns and features in the rows.
#'   Or a `SummarizedExperiment` / `SingleCellExperiment` object. The features must
#'   match the features in `fit`.
#' @param col_data col_data an optional data frame with `ncol(data)` rows.
#' @param use_assay if `data` is a `SummarizedExperiment` / `SingleCellExperiment` object,
#'   which assay should be used.
#' @param design,alignment_design the design formulas or design matrices that are used
#'   to project the data on the correct latent subspace. Both default to the designs
#'   from the `fit` object.
#' @param return which data structure is returned.
#'
#'
#' @returns either a matrix with the low-dimensional embeddings of the `data` or
#'   an object of class `lemur_fit` wrapping that embedding.
#'
#'
#' @export
project_on_lemur_fit <- function(fit, data, col_data = NULL, use_assay = "logcounts",
                                 design = fit$design, alignment_design = fit$alignment_design,
                                 return = c("matrix", "lemur_fit")){
  return <- match.arg(return)
  Y <- handle_data_parameter(data, on_disk = FALSE, assay = use_assay)
  col_data <- glmGamPoi:::get_col_data(data, col_data)
  des <- handle_design_parameter(design, data, col_data)
  design_matrix <- des$design_matrix
  al_des <- handle_design_parameter(design, data, col_data)
  al_design_matrix <- al_des$design_matrix
  Y_clean <- Y - fit$linear_coefficients %*% t(design_matrix)
  embedding <- project_data_on_diffemb(Y_clean, design = design_matrix, coefficients = fit$coefficients, base_point = fit$base_point)
  embedding <- apply_linear_transformation(embedding, fit$alignment_coefficients, al_design_matrix)
  colnames(embedding) <- colnames(data)

  if(return == "matrix"){
    embedding
  }else if(return == "lemur_fit"){
    lemur_fit(data, col_data = col_data,
              row_data = if(is(data, "SummarizedExperiment")) rowData(data) else NULL,
              n_embedding = fit$n_embedding,
              design = des$design_formula, design_matrix = design_matrix,
              linear_coefficients = fit$linear_coefficients,
              base_point = fit$base_point,
              coefficients = fit$coefficients,
              embedding = embedding,
              alignment_coefficients = fit$alignment_coefficients,
              alignment_design = al_des$design_formula,
              alignment_design_matrix = al_design_matrix,
              use_assay = use_assay, test_data = NULL)
  }
}


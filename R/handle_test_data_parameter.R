

handle_test_data_parameter <- function(fit, test_data, test_data_col_data, continuous_assay_name){
  if(is(test_data, "SummarizedExperiment")){
    if(! continuous_assay_name %in% assayNames(test_data)){
      stop("Cannot find assay '", continuous_assay_name, "' in the assays of 'independent_data'")
    }
    # Nothing else to be done here
  }else if(is.list(test_data)){
    if(! continuous_assay_name %in% names(test_data)){
      stop("Cannot find assay '", continuous_assay_name, "' in the names of 'independent_data'")
    }
    if(is.null(test_data_col_data)){
      stop("'independent_data_col_data' must not be NULL")
    }
    test_data <- SingleCellExperiment(assays = test_data)
  }else if(is.matrix(test_data)){
    message("'independent_data' is a matrix treating it as continuous values")
    if(is.null(test_data_col_data)){
      stop("'independent_data_col_data' must not be NULL")
    }
    test_data <- SingleCellExperiment(assays = setNames(list(test_data), continuous_assay_name))
  }else if(is.null(test_data)){
    # This is necessary to satisfy model.matrix in 'project_on_lemur_fit'
    col_data_copy <- fit$colData
    character_cols <- vapply(col_data_copy, is.character, logical(1L))
    col_data_copy[character_cols] <- lapply(col_data_copy[character_cols], as.factor)
    test_data <- SingleCellExperiment(assays = setNames(list(matrix(nrow = nrow(fit), ncol = 0) * 1.0), continuous_assay_name),
                                      colData = col_data_copy[integer(0L),,drop=FALSE])
  }else{
    stop("Cannot handle 'indepdendet_data' of type: ", paste0(class(test_data), collapse = ", "))
  }

  colData(test_data) <- S4Vectors::DataFrame(glmGamPoi:::get_col_data(test_data, test_data_col_data))

  test_data
}

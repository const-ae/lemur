
handle_design_parameter <- function(design, data, col_data){
  n_samples <- ncol(data)

  ignore_degeneracy <- isTRUE(attr(design, "ignore_degeneracy"))

  # Handle the design parameter
  if(is.matrix(design)){
    design_matrix <- design
    design_formula <- NULL
  }else if((is.vector(design) || is.factor(design))){
    if(length(design) != n_samples){
      if(length(design) == 1 && design == 1){
        stop("The specified design vector length (", length(design), ") does not match ",
             "the number of samples: ", n_samples, "\n",
             "Did you maybe mean: `design = ~ 1`?")
      }else{
        stop("The specified design vector length (", length(design), ") does not match ",
             "the number of samples: ", n_samples)
      }
    }
    tmp <- glmGamPoi:::convert_chr_vec_to_model_matrix(design, NULL)
    design_matrix <- tmp$model_matrix
    design_formula <- NULL
  }else if(inherits(design,"formula")){
    tmp <- convert_formula_to_design_matrix(design, col_data)
    design_matrix <- tmp$design_matrix
    design_formula <- tmp$formula
    attr(design_formula, "constructed_from") <- "formula"
  }else{
    stop("design argment of class ", class(design), " is not supported. Please ",
         "specify a `design_matrix`, a `character vector`, or a `formula`.")
  }

  if(nrow(design_matrix) != ncol(data)) stop("Number of rows in col_data does not match number of columns of data.")
  if(! is.null(rownames(design_matrix)) &&
     ! all(rownames(design_matrix) == as.character(seq_len(nrow(design_matrix)))) && # That's the default rownames
     ! is.null(colnames(data))){
    if(! all(rownames(design_matrix) == colnames(data))){
      if(setequal(rownames(design_matrix), colnames(data))){
        # Rearrange the rows to match the columns of data
        design_matrix <- design_matrix[colnames(data), ,drop=FALSE]
      }else{
        stop("The rownames of the design_matrix / col_data do not match the column names of data.")
      }
    }
  }

  if(any(matrixStats::rowAnyNAs(design_matrix))){
    stop("The design matrix contains 'NA's for sample ",
         paste0(head(which(DelayedMatrixStats::rowAnyNAs(design_matrix))), collapse = ", "),
         ". Please remove them before you call 'lemur()'.")
  }

  if(ncol(design_matrix) >= n_samples && ! ignore_degeneracy){
    stop("The design_matrix has more columns (", ncol(design_matrix),
         ") than the there are samples in the data matrix (", n_samples, " columns).\n",
         "Too few replicates / too many coefficients to fit model.\n",
         "The head of the design matrix: \n", glmGamPoi:::format_matrix(head(design_matrix, n = 3)))
  }

  # Check rank of design_matrix
  qr_mm <- qr(design_matrix)
  if(qr_mm$rank < ncol(design_matrix) && n_samples > 0  && ! ignore_degeneracy){
    is_zero_column <- DelayedMatrixStats::colCounts(design_matrix, value = 0) == nrow(design_matrix)
    if(any(is_zero_column)){
      stop("The model matrix seems degenerate ('matrix_rank(design_matrix) < ncol(design_matrix)'). ",
           "Column ", paste0(head(which(is_zero_column), n=10), collapse = ", "), " contains only zeros. \n",
           "The head of the design matrix: \n", glmGamPoi:::format_matrix(head(design_matrix, n = 3)))
    }else{
      stop("The model matrix seems degenerate ('matrix_rank(design_matrix) < ncol(design_matrix)'). ",
           "Some columns are perfectly collinear. Did you maybe include the same coefficient twice?\n",
           "The head of the design matrix: \n", glmGamPoi:::format_matrix(head(design_matrix, n = 3)))
    }
  }

  rownames(design_matrix) <- colnames(data)
  validate_design_matrix(design_matrix, data)
  # design_matrix <- add_attr_if_intercept(design_matrix)
  list(design_matrix = design_matrix, design_formula = design_formula)
}



convert_formula_to_design_matrix <- function(formula, col_data){
  attr(col_data, "na.action") <- "na.pass"
  tryCatch({
    mf <- model.frame(formula, data = col_data, drop.unused.levels = FALSE)
    terms <- attr(mf, "terms")
    attr(terms, "xlevels") <- stats::.getXlevels(terms, mf)
    mm <- stats::model.matrix.default(terms, mf)
    attr(terms, "contrasts") <- attr(mm, "contrasts")
  }, error = function(e){
    # Try to extract text from error message
    match <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match) == 2){
      stop("Error while parsing the formula (", formula, ").\n",
           "Variable '", match[2], "' not found in col_data or global environment. Possible variables are:\n",
           paste0(colnames(col_data), collapse = ", "), call. = FALSE)
    }else{
      stop(e$message)
    }
  })

  # Otherwise every copy of the model stores the whole global environment!
  attr(terms, ".Environment") <- c()
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  list(formula = terms, design_matrix = mm)
}

validate_design_matrix <- function(matrix, data){
  stopifnot(is.matrix(matrix))
  stopifnot(nrow(matrix) == ncol(data))
}


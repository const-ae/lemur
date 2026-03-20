
handle_design_parameter <- function(design, data, col_data, verbose = FALSE){
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
    col_data <- add_global_variables_to_col_data(design, col_data)
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

  if(verbose && is.null(design_formula)){
    message("The 'design' was not specified with a formula. This means that you cannot use 'cond(...)' in 'compute_contrasts(...)'.")
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
  list(design_matrix = design_matrix, design_formula = design_formula, col_data = col_data)
}



convert_formula_to_design_matrix <- function(formula, col_data){
  attr(col_data, "na.action") <- "na.pass"
  tryCatch({
    mf <- model.frame(formula, data = col_data, drop.unused.levels = FALSE)
    terms <- attr(mf, "terms")
    # xlevels is used for reconstructing the model matrix
    attr(terms, "xlevels") <- stats::.getXlevels(terms, mf)
    # vars_xlevels is used to check input to cond(...)
    attr(terms, "vars_xlevels") <- xlevels_for_formula_vars(terms, col_data)
    mm <- stats::model.matrix.default(terms, mf)
    attr(terms, "contrasts") <- attr(mm, "contrasts")
  }, error = function(e){
    # Try to extract text from error message
    match <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match) == 2){
      stop("Problem parsing the formula (", formula, ").\n",
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

add_global_variables_to_col_data <- function(formula, col_data){
  # Check if var is global and put it into col_data
  formula_env <- attr(formula, ".Environment")
  if(is.null(formula_env)) formula_env <- rlang::empty_env()
  for(gv in setdiff(all.vars(formula), colnames(col_data))){
    value <- rlang::eval_tidy(rlang::sym(gv), data = NULL, env = formula_env)
    is_vector_type <- (is(col_data, "DFrame") && is(value, "Vector")) || vctrs::obj_is_vector(value)
    if(is_vector_type){
      has_correct_length <- NROW(value) == 1 || NROW(value) == nrow(col_data)
      if(! has_correct_length){
        stop("Trying store global variables from formula in colData, however '", gv, "' ",
             "has length ", NROW(value), ", but it needs to be 1 or nrow(col_data) (", nrow(col_data), ").")
      }
    }else{
      stop("Trying store global variables from formula in colData, however '", gv, "' ",
           "is of type ",  class(value)[1]," and not a vector-type.")
    }
    col_data[[gv]] <- value
  }
  col_data
}

xlevels_for_formula_vars <- function(formula, data){
  # if(! is.null( attr(formula, "xlevels"))){
  #   attr(formula, "xlevels")
  if(! is.null( attr(formula, "vars_xlevels"))){
    attr(formula, "vars_xlevels")
  }else{
    # For all character / factor vars get xlevel
    all_vars <- all.vars(formula)
    formula_env <- attr(formula, ".Environment")
    if(is.null(formula_env)) formula_env <- rlang::empty_env()
    xlev <- lapply(all_vars, \(v){
      value <- rlang::eval_tidy(rlang::sym(v), data = as.data.frame(data), env = formula_env)
      if(is.character(value)){
        levels(as.factor(value))
      }else if(is.factor(value)){
        levels(value)
      }else{
        NULL
      }
    })
    names(xlev) <- all_vars
    xlev[!vapply(xlev, is.null, TRUE)]
  }
}

validate_design_matrix <- function(matrix, data){
  stopifnot(is.matrix(matrix))
  stopifnot(nrow(matrix) == ncol(data))
}


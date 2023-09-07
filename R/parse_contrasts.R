
# This function is adapted from proDA
parse_contrast <- function(contrast, formula) {

  if(missing(contrast)){
    stop("No contrast argument was provided!")
  }
  covar <- all.vars(formula)

  cnt_capture <- rlang::enquo(contrast)
  if(is.null(formula)){
    data_mask <- NULL
  }else{
    data_mask <- create_contrast_data_mask(formula)
  }

  tryCatch({
    res <- rlang::eval_tidy(cnt_capture, data = data_mask)
    if(! is.numeric(res)){
      if(is.character(res)){
        # If contrast was a string, eval will just spit it out the same way
        res <- rlang::eval_tidy(rlang::parse_expr(res), data = data_mask)
      }
    }
  }, error = function(e){
    # Try to extract text from error message
    match <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match) == 2){
      stop("Object '", match[2], "' not found. Please specify the contrast using:\n",
           "'cond(", paste0(paste0(covar, " = ?"), collapse = ", "), ") - ",
           "cond(", paste0(paste0(covar, " = ?"), collapse = ", "), ")'", call. = FALSE)
    }else{
      stop(e$message)
    }
  })
  res
}


create_contrast_data_mask <- function(formula){
  top <- rlang::new_environment(list(
    cond = function(...){
      .cond(formula, list(...))
    },
    "+" = .plus, "-" = .minus, "/" = .divide, "*" = .multiply,
    "==" = .equal,
    "<" = .lt, "<=" = .lt,
    ">" = .gt, ">=" = .gt
  ))
  bottom <- rlang::new_environment(parent = top)
  data_mask <- rlang::new_data_mask(bottom = bottom, top = top)
  data_mask$.cntrst <- rlang::as_data_pronoun(bottom)
  data_mask
}


.cond <- function(formula, level_sets = list()){
  if(is.null(formula)){
    stop("You called 'cond()' inside the contrast, however the original model ",
         "was not specified with a formula. Thus 'cond()' doesn't work and you ",
         "need to specify the contrast using the column names of the design matrix.")
  }
  if(is.null(attr(formula, "xlevels"))){
    warning("The formula has no 'xlevels' attribute. This is supicious and might indicate a bug.")
  }
  if(any(names(level_sets) == "")){
    stop("All arguments to 'cond()' must be named.")
  }
  if(any(duplicated(names(level_sets)))){
    stop("All arguments to 'cond()' must be unique.")
  }

  covar <- all.vars(formula)
  new_dat <- as.list(rep(0, length(covar)))
  names(new_dat) <- covar
  xlevels <- attr(formula, "xlevels")
  for(n in names(xlevels)){
    new_dat[[n]] <- factor(xlevels[[n]][1], levels = xlevels[[n]])
  }
  for(n in names(level_sets)){
    if(! n %in% names(new_dat)){
      stop("Setting the level of '", n, "' failed. You can only set the level of the following variables: ", paste0(covar, collapse = ", "))
    }
    if(length(level_sets[[n]]) != 1){
      stop("Each argument to 'cond()' must be length one. '", n, "' has length ", length(level_sets[[n]]))
    }
    if(n %in% names(xlevels)){
      if(! level_sets[[n]] %in% xlevels[[n]]){
        stop("You are trying to set '", n, "=", level_sets[[n]], "'. However only the following values for ", n,
             " are valid: ", paste0(xlevels[[n]], collapse = ", "))
      }
      new_dat[[n]] <- factor(level_sets[[n]], levels = xlevels[[n]])
    }else{
      new_dat[[n]] <- level_sets[[n]]
    }
  }
  res <- drop(model.matrix(formula, new_dat, contrasts.arg = attr(formula, "contrasts")))
  attr(res, "assign") <- NULL
  attr(res, "contrasts") <- NULL
  class(res) <- "model_vec"
  res
}


evaluate_contrast_tree <-function(c1, c2, FUN){
  stopifnot(all(class(c2) == class(c1)))
  if(inherits(c1, "contrast_relation")){
    stopifnot(c1$relation == c2$relation)
    if(c1$relation == "minus" && is.null(c1$rhs)){ # Unary minus
      - evaluate_contrast_tree(c1$lhs, c2$lhs, FUN = FUN)
    }else if(c1$relation == "minus" && is.null(c1$rhs)){ # Unary plus
      + evaluate_contrast_tree(c1$lhs, c2$lhs, FUN = FUN)
    }else if(c1$relation == "minus"){
      evaluate_contrast_tree(c1$lhs, c2$lhs, FUN = FUN) - evaluate_contrast_tree(c1$rhs, c2$rhs, FUN = FUN)
    }else if(c1$relation == "plus"){
      evaluate_contrast_tree(c1$lhs, c2$lhs, FUN = FUN) + evaluate_contrast_tree(c1$rhs, c2$rhs, FUN = FUN)
    }else if(c1$relation == "multiply"){
      evaluate_contrast_tree(c1$lhs, c2$lhs, FUN = FUN) * evaluate_contrast_tree(c1$rhs, c2$rhs, FUN = FUN)
    }else if(c1$relation == "divide"){
      evaluate_contrast_tree(c1$lhs, c2$lhs, FUN = FUN) / evaluate_contrast_tree(c1$rhs, c2$rhs, FUN = FUN)
    # }else if(c1$relation == "equal"){
    #   evaluate_contrast_tree(c1$lhs, c2$lhs, FUN = FUN) - evaluate_contrast_tree(c1$rhs, c2$rhs, FUN = FUN)
    }else if(c1$relation %in% c("equal", "less_than", "greater_than")){
      stop("(In)equalities are not allowed in contrasts")
    }else{
      stop("Canot handle contrast relationship of type: ", c1$relation)
    }
  }else if(inherits(c1, "model_vec")){
    FUN(c1, c2)
  }else{
    stopifnot(all(c1 == c2))
    c1
  }
}


.divide <- function(x, y){
  res <- list(lhs = x, rhs = y, relation =  "divide")
  class(res) <- "contrast_relation"
  res
}

.multiply <- function(x, y){
  res <- list(lhs = x, rhs = y, relation =  "multiply")
  class(res) <- "contrast_relation"
  res
}

.plus <- function(x, y){
  if(missing(y)){
    # Unary plus
    y <- NULL
  }
  res <- list(lhs = x, rhs = y, relation =  "plus")
  class(res) <- "contrast_relation"
  res
}

.minus <- function(x, y){
  if(missing(y)){
    # Unary minus
    y <- NULL
  }
  res <- list(lhs = x, rhs = y, relation =  "minus")
  class(res) <- "contrast_relation"
  res
}

.equal <- function(x, y){
  res <- list(lhs = x, rhs = y, relation =  "equal")
  class(res) <- "contrast_relation"
  res
}

.lt <- function(x, y){
  res <- list(lhs = x, rhs = y, relation =  "less_than")
  class(res) <- "contrast_relation"
  res
}

.gt <- function(x, y){
  res <- list(lhs = x, rhs = y, relation = "greater_than")
  class(res) <- "contrast_relation"
  res
}



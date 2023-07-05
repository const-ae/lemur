
handle_data_parameter <- function(data, on_disk, assay){
  if(is.matrix(data)){
    if(! is.numeric(data)){
      stop("The data argument must consist of numeric values and not of ", mode(data), " values")
    }
    if(is.null(on_disk) || isFALSE(on_disk)){
      data_mat <- data
    }else if(isTRUE(on_disk)){
      data_mat <- HDF5Array::writeHDF5Array(data)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else if(is(data, "DelayedArray")){
    if(is.null(on_disk) || isTRUE(on_disk)){
      data_mat <- data
    }else if(isFALSE(on_disk)){
      data_mat <- as.matrix(data)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else if(is(data, "SummarizedExperiment")){
    data_mat <- handle_data_parameter(SummarizedExperiment::assay(data, assay), on_disk)
  }else if(canCoerce(data, "SummarizedExperiment")){
    se <- as(data, "SummarizedExperiment")
    data_mat <- handle_data_parameter(SummarizedExperiment::assay(se, assay), on_disk)
  }else if(is(data, "dgCMatrix") || is(data, "dgTMatrix")) {
    if(isTRUE(on_disk)){
      data_mat <- HDF5Array::writeHDF5Array(data)
    }else if(isFALSE(on_disk)){
      data_mat <- as.matrix(data)
    }else{
      stop("glmGamPoi does not yet support sparse input data of type '", class(data),"'. ",
           "Please explicitly set the 'on_disk' parameter to force a conversion to a dense format either in-memory ('on_disk = FALSE') ",
           "or on-disk ('on_disk = TRUE')")
    }
  }else{
    stop("Cannot handle data of class '", class(data), "'.",
         "It must be of type dense matrix object (i.e., a base matrix or DelayedArray),",
         " or a container for such a matrix (for example: SummarizedExperiment).")
  }
  data_mat
}


convert_dataframe_cols_chr_to_fct <- function(col_data){
  character_cols <- vapply(col_data, is.character, logical(1L))
  col_data[character_cols] <- lapply(col_data[character_cols], as.factor)
  col_data
}


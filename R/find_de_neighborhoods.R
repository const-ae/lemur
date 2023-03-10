find_random_directions <- function(fit, n){
  point_pairs <- sapply(seq_len(n), \(idx) sample.int(ncol(fit), 2))

  dirs <- lapply(seq_len(ncol(point_pairs)), \(idx){
    sel <- point_pairs[,idx]
    vec <- fit$embedding[,sel[1]] - fit$embedding[,sel[2]]
    vec / sqrt(sum(vec^2))
  })
  attr(dirs, "point_pairs") <- point_pairs
  dirs
}

find_de_neighborhoods_with_random_projections_zscore <- function(fit, de_mat, random_projections = 50, include_complement = FALSE){
  stopifnot(all(dim(fit) == dim(de_mat)))
  n_genes <- nrow(fit)
  n_cells <- ncol(fit)

  if(is.numeric(random_projections) && length(random_projections) == 1){
    dirs <- find_random_directions(fit, random_projections)
  }else if(is.list(random_projections)){
    dirs <- random_projections
  }else{
    stop("Cannot handle 'random_projections' of type: ", class(random_projections))
  }
  emb <- fit$embedding
  proj <- do.call(rbind, lapply(dirs, \(dir) c(coef(lm.fit(x = matrix(dir, ncol = 1), emb)))))
  rownames(proj) <- paste0("proj_", seq_len(nrow(proj)))


  # for(region_idx in seq_len(n_regions_per_gene))
  region_idx <- 1

  # Find the most helpful projection for each gene
  correlation <- cor(t(proj), t(de_mat))
  best_proj <- t(apply(correlation, 2, \(col) order(-abs(col))))

  de_region_indices <- lapply(seq_len(n_genes), \(gene_idx){
    pr <- proj[best_proj[gene_idx, 1], ,drop=TRUE]
    order_pr <- order(pr)
    max_idx <- cumz_which_abs_max(de_mat[gene_idx,order_pr])
    rev_max_idx <- cumz_which_abs_max(rev(de_mat[gene_idx,order_pr]))
    if(max_idx$max > rev_max_idx$max){
      which(pr < pr[order_pr][max_idx$idx])
    }else{
      which(pr > rev(pr[order_pr])[rev_max_idx$idx])
    }
  })
  region_means <- vapply(seq_len(n_genes), \(gene_idx) mean(de_mat[gene_idx, de_region_indices[[gene_idx]]]), FUN.VALUE = numeric(1L))

  feature_names <- if(is.null(rownames(de_mat))){
    paste0("feature_", seq_len(nrow(de_mat)))
  }else{
    rownames(de_mat)
  }
  result <- data.frame(name = feature_names,
                       region = as.character(rep(region_idx, times = n_genes)),
                       indices = I(de_region_indices),
                       n_cells = lengths(de_region_indices),
                       mean = region_means)
  if(include_complement){
    all_indices <- seq_len(n_cells)
    complement_regions <- lapply(de_region_indices, \(indices) setdiff(all_indices, indices))
    complement_region_means <- vapply(seq_len(nrow(fit)), \(gene_idx) mean(de_mat[gene_idx, complement_regions[[gene_idx]]]), FUN.VALUE = numeric(1L))
    complement_result <- data.frame(name = feature_names,
                                    region = paste0("complement_to_", rep(region_idx, times = n_genes)),
                                    indices = I(complement_regions),
                                    n_cells = lengths(complement_regions),
                                    mean = complement_region_means)
    rbind(result, complement_result)
  }else{
    result
  }
}



#' Find differential expression neighborhoods
#'
#' @param fit the `lemur_fit` generated by `lemur()`
#' @param de_mat the matrix with the differential expression values. Defaults
#'   to an assay called `"DE"` that is produced by `lemur::test_de()`.
#' @param counts a count matrix with the same dimension as `fit`.
#'   Optional, if provided the function calculates the pseudobulked differential
#'   expression.
#' @param group_by if count is provided you need to specify how the pseudobulk
#'   samples are formed.
#' @param contrast a specification which contrast to fit. This defaults to the
#'   `contrast` argument that was used for `test_de` and is stored in `fit$contrast`.
#' @param design the design to use for the fit. Default: `fit$design`
#' @param random_projections the number of random projections to use
#' @param include_complement a boolean to specify if the complement to the identified
#'   neighborhood is also returned
#' @param verbose Should the method print information during the fitting. Default: `TRUE`.
#'
#' @return a data frame with one entry per gene / neighborhood containing the name
#'   of the neighborhood, the cell indices included in the neighborhood, the number of
#'   cells, and the mean. If counts is not `NULL`, the data frame in addition contains
#'   the output columns from `glmGamPoi::test_de()`
#'
#' @export
find_de_neighborhoods <- function(fit, de_mat = assay(fit, "DE"), counts = NULL, group_by,
                                  contrast = fit$contrast, design = fit$design,
                                  random_projections = 50, include_complement = TRUE,
                                  verbose = TRUE){

  tryCatch(force(de_mat), error = function(e){
    stop("Cannot find assay `\"DE\"`. Please first call 'lemur::test_de()' to calculate the differential expression matrix.")
  })

  de_regions <- find_de_neighborhoods_with_random_projections_zscore(fit, de_mat, random_projections = random_projections, include_complement = include_complement)

  if(is.null(counts)){
    de_regions
  }else{
    stopifnot(all(dim(fit) == dim(counts)))
    if(! is.null(rownames(fit)) && is.null(rownames(counts))){
      rownames(counts) <- rownames(fit)
    }
    if(any(rownames(fit) != rownames(counts))){
      stop("The rownames of fit and counts don't match.")
    }
    if(rlang::quo_is_null(rlang::enquo(contrast))){
      stop("The contrast argument is 'NULL'. Please specify.")
    }
    tryCatch({
      if(inherits(contrast, "contrast_relation")){
        contrast <- evaluate_contrast_tree(contrast, contrast, \(x, y) x)
      }
    }, error = function(e){
      # Do nothing. The 'contrast' is probably an unquoted expression
    })


    separator <- "-__-"
    # mask <- Matrix::sparseMatrix(i = integer(0L), j = integer(0L), x = numeric(0L), dims = c(nrow(de_regions),  ncol(fit)),
    #                      dimnames = list(paste0(de_regions$name, separator, de_regions$region), colnames(fit)),
    #                      repr = "T")
    mask <- matrix(0, nrow = nrow(de_regions),  ncol = ncol(fit),
                   dimnames = list(paste0(de_regions$name, separator, de_regions$region), colnames(fit)))
    for(idx in seq_len(nrow(de_regions))){
      mask[idx,de_regions$indices[[idx]]] <- 1
    }

    mask <- if(utils::packageVersion("Matrix") >= "1.4.2"){
      # See email from Martin Maechler from 2022-08-12
      as(as(as(mask, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    }else{
      # This approach is deprecated since 1.4.2 and triggers warnings
      as(mask, "dgTMatrix")
    }

    if(is.null(rownames(counts))){
      rownames(counts) <- paste0("feature_", seq_len(nrow(de_mat)))
    }
    masked_counts <- as(unname(counts[de_regions$name,,drop=FALSE]) * mask, "CsparseMatrix")
    masked_size_factors <- Matrix::sparseMatrix(i = mask@i, j = mask@j, x = mask@x * MatrixGenerics::colSums2(counts)[mask@j + 1L], index1 = FALSE, repr = "C")

    masked_sce <- SingleCellExperiment::SingleCellExperiment(list(masked_counts = masked_counts, masked_size_factors = masked_size_factors),
                                                             colData = SingleCellExperiment::colData(fit))
    region_psce <- glmGamPoi::pseudobulk(masked_sce, group_by = {{group_by}},
                                         aggregation_functions = list("masked_counts" = "rowSums2",
                                                                      "masked_size_factors" = "rowSums2"),
                                         verbose = verbose)

    glm_regions <- glmGamPoi::glm_gp(region_psce, design = design, use_assay = "masked_counts", verbose = verbose,
                                     offset = log(assay(region_psce, "masked_size_factors") + 1e-10),
                                     size_factors = FALSE, overdispersion = TRUE)
    de_res <- glmGamPoi::test_de(glm_regions, contrast = {{contrast}})
    cbind(de_regions, de_res[,-1])
  }
}

#' @importFrom glmGamPoi vars
#' @export
glmGamPoi::vars

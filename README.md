
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Latent Embedding Multivariate Regression (LEMUR)

<!-- badges: start -->
<!-- badges: end -->

The goal of lemur is to enable easy analysis of multi-condition
single-cell data. Lemur fits latent embedding regression, which means
that it tries to find the PCA embedding for each condition and
parameterizes the transition from one embedding to the other. For this
task, lemur uses geodesic regression on the Grassmann manifold, which is
solved efficiently using tangent-space linear modelling. The result is
an interpretable model of the gene expression for arbitrary experimental
designs that can be expressed using a design matrix.

## Installation

You can install the released version of lemur from Github

``` r
devtools::install_github("const-ae/lemur")
```

## Example

Make some data

``` r
mat <- matrix(rnorm(30 * 500), nrow = 30, ncol = 500)
col_data <- data.frame(condition = sample(letters[1:2], 500, replace = TRUE))
```

Fit the model

``` r
fit <- lemur::lemur(mat, design = ~ condition, col_data = col_data, n_ambient = 10, n_embedding = 2)
#> Fit ambient PCA
#> Regress out global effects
#> Find base point for differential embedding
#> Fit differential embedding model
#> -Iteration: 0    error: 6.5e+03
#> ---Fit Grassmann linear model
#> ---Update linear regression
#> -Iteration: 1    error: 4.85e+03
#> ---Fit Grassmann linear model
#> ---Update linear regression
#> -Iteration: 2    error: 4.85e+03
#> Converged
fit
#> class: lemur_fit_obj 
#> dim: 30 500 
#> metadata(12): n_ambient n_embedding ... alignment_design
#>   alignment_design_matrix
#> assays(1): expr
#> rownames: NULL
#> rowData names(0):
#> colnames: NULL
#> colData names(1): condition
#> reducedDimNames(2): linearFit embedding
#> mainExpName: NULL
#> altExpNames(0):
```

Let’s look at the coefficients

``` r
round(fit$coefficients, 5)
#> , , Intercept
#> 
#>           [,1]     [,2]
#>  [1,] -0.00077  0.00096
#>  [2,] -0.00308  0.00715
#>  [3,]  0.33556  0.09131
#>  [4,]  0.54514  0.85025
#>  [5,] -0.44563 -0.21422
#>  [6,]  0.34480 -0.01032
#>  [7,] -0.12499 -0.12281
#>  [8,]  0.00143  0.38604
#>  [9,]  0.19342  0.17167
#> [10,]  0.05171  0.05353
#> 
#> , , conditionb
#> 
#>           [,1]     [,2]
#>  [1,]  0.00114  0.00056
#>  [2,]  0.00163  0.01919
#>  [3,] -0.58263  0.43780
#>  [4,] -0.46199 -1.20915
#>  [5,]  0.65308 -0.05886
#>  [6,] -0.57180 -0.07554
#>  [7,]  0.12282  0.23320
#>  [8,]  0.03196 -0.58948
#>  [9,] -0.24322  0.05567
#> [10,]  0.05912 -0.07358
plot(t(fit$embedding), col = as.factor(col_data$condition))
```

<img src="man/figures/README-show-fit-results-1.png" width="100%" />

``` r
de_mat <- lemur::test_de(fit, contrast = cond(condition = "a") == cond(condition = "b"))
de_mat[1:5, 1:10]
#>             [,1]         [,2]        [,3]        [,4]        [,5]       [,6]
#> [1,]  0.54273594  0.227746766 -0.07380608  0.13953205 -0.40450075  0.6948219
#> [2,] -0.28660964 -0.087217712  0.25464570  0.05567957  0.25241819 -0.2824847
#> [3,]  0.32969179  0.118329964 -0.25143408 -0.03739997 -0.23873435  0.3204130
#> [4,]  0.10420045  0.006754924 -0.09103860 -0.02313281 -0.18703047  0.1482544
#> [5,]  0.04336833 -0.026404799 -0.40268196 -0.22439783 -0.04226467 -0.1287427
#>             [,7]       [,8]        [,9]       [,10]
#> [1,]  0.56702390 -0.2620746  0.05170869  0.01677785
#> [2,] -0.29499019  0.2234453 -0.11071287  0.06736877
#> [3,]  0.33823367 -0.2110121  0.14982926 -0.04656486
#> [4,]  0.11150553 -0.1447949 -0.04367860 -0.05913853
#> [5,]  0.03441165 -0.1361353  0.21118642 -0.11626663
```

Show the gene expression changes on the latent embedding

``` r
pal <- scales::col_numeric(scales::viridis_pal()(50), NULL)
plot(t(fit$embedding), col = pal(de_mat[10,]), pch = 16, cex = 1.3)
```

<img src="man/figures/README-visualize-differential-expression-1.png" width="100%" />

``` r
counts <- matrix(rpois(30 * 500, lambda = 2.4), nrow = 30, ncol = 500)
SummarizedExperiment::colData(fit)$patient_id <- sample(c("A", "B", "C"), size = 500, replace = TRUE)
de_regions <- lemur::find_de_neighborhoods(fit, de_mat, counts = counts, group_by = glmGamPoi::vars(patient_id, condition), 
                             contrast = cond(condition = "a") == cond(condition = "b"))
#> Aggregating assay 'masked_counts' using 'rowSums2'.
#> Aggregating assay 'masked_size_factors' using 'rowSums2'.
#> Make initial dispersion estimate
#> Make initial beta estimate
#> Estimate beta
#> Estimate dispersion
#> Fit dispersion trend
#> Shrink dispersion estimates
#> Estimate beta again
head(de_regions)
#>        name region      indices n_cells       mean         pval     adj_pval
#> 1 feature_1      1 1, 6, 7,....      98  0.6973596 0.000000e+00 0.000000e+00
#> 2 feature_2      1      52, 222       2 -0.8845059 4.174344e-31 6.591070e-31
#> 3 feature_3      1      52, 222       2  0.9624984 4.377569e-26 5.836759e-26
#> 4 feature_4      1 5, 11, 1....     123 -0.2583617 0.000000e+00 0.000000e+00
#> 5 feature_5      1           21       1  0.6887578 3.496154e-23 4.034023e-23
#> 6 feature_6      1          191       1  1.9229842 9.999814e-01 9.999814e-01
#>    f_statistic df1       df2        lfc
#> 1 4.799248e+03   1 485790184  -4.715161
#> 2 1.345344e+02   1 485790184 -10.000000
#> 3 1.115974e+02   1 485790184  -4.624491
#> 4 6.959848e+03   1 485790184  -5.000324
#> 5 9.835558e+01   1 485790184  -4.977280
#> 6 5.454098e-10   1 485790184 -10.000000
```

# Session Info

``` r
sessionInfo()
#> R version 4.2.1 RC (2022-06-17 r82503)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur ... 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.10                 compiler_4.2.1             
#>  [3] GenomeInfoDb_1.34.9         highr_0.10                 
#>  [5] XVector_0.38.0              DelayedMatrixStats_1.20.0  
#>  [7] MatrixGenerics_1.10.0       bitops_1.0-7               
#>  [9] tools_4.2.1                 zlibbioc_1.44.0            
#> [11] SingleCellExperiment_1.20.0 digest_0.6.31              
#> [13] viridisLite_0.4.1           evaluate_0.20              
#> [15] lifecycle_1.0.3             lattice_0.20-45            
#> [17] rlang_1.0.6                 Matrix_1.5-3               
#> [19] DelayedArray_0.24.0         cli_3.6.0                  
#> [21] rstudioapi_0.14             yaml_2.3.7                 
#> [23] expm_0.999-7                xfun_0.37                  
#> [25] fastmap_1.1.1               GenomeInfoDbData_1.2.9     
#> [27] lemur_0.0.6                 knitr_1.42                 
#> [29] vctrs_0.5.2                 S4Vectors_0.36.2           
#> [31] IRanges_2.32.0              stats4_4.2.1               
#> [33] grid_4.2.1                  Biobase_2.58.0             
#> [35] R6_2.5.1                    rmarkdown_2.20             
#> [37] irlba_2.3.5.1               farver_2.1.1               
#> [39] sparseMatrixStats_1.10.0    scales_1.2.1               
#> [41] matrixStats_0.63.0          htmltools_0.5.4            
#> [43] BiocGenerics_0.44.0         GenomicRanges_1.50.2       
#> [45] SummarizedExperiment_1.28.0 colorspace_2.1-0           
#> [47] glmGamPoi_1.11.6            RCurl_1.98-1.10            
#> [49] munsell_0.5.0
```

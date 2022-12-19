
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DiffEmbSeq

<!-- badges: start -->
<!-- badges: end -->

The goal of DiffEmbSeq is to enable easy analysis of multi-condition
single-cell data. DiffEmbSeq fits a differential embedding model, which
means that it tries to find the PCA embedding for each condition and
parameterizes the transition from one embedding to the other. For this
task, DiffEmbSeq uses geodesic regression on the Grassmann manifold,
which is solved efficiently using tangent-space linear modelling. The
result is an interpretable model of the gene expression for arbitrary
experimental designs that can be expressed using a design matrix.

## Installation

You can install the released version of DiffEmbSeq from
<https://git.embl.de/ahlmanne/DiffEmbSeq>

To get read access you have set the following *access token* which is
valid until the 25.02.2023.

``` r
Sys.setenv(GITLAB_PAT = "glpat-vE2cB9xNSjsU53d4z16f")
remotes::install_gitlab(repo = "ahlmanne/DiffEmbSeq", host = "https://git.embl.de/")
```

## Example

Make some data

``` r
mat <- matrix(rnorm(30 * 500), nrow = 30, ncol = 500)
col_data <- data.frame(condition = sample(letters[1:2], 500, replace = TRUE))
```

Fit the model

``` r
fit <- DiffEmbSeq::differential_embedding(mat, design = ~ condition, col_data = col_data,
                                          n_ambient = 10, n_embedding = 2)
#> Fit ambient PCA
#> Regress out global effects
#> Find base point for differential embedding
#> Fit differential embedding model
#> -Iteration: 0    error: 6.54e+03
#> ---Fit Grassmann linear model
#> ---Update linear regression
#> -Iteration: 1    error: 4.88e+03
#> ---Fit Grassmann linear model
#> ---Update linear regression
#> -Iteration: 2    error: 4.88e+03
#> Converged
fit
#> class: DiffEmbFit 
#> dim: 30 500 
#> metadata(11): n_ambient n_embedding ... alignment_stretching ''
#> assays(1): expr
#> rownames: NULL
#> rowData names(0):
#> colnames: NULL
#> colData names(1): condition
#> reducedDimNames(2): linearFit diffemb_embedding
#> mainExpName: NULL
#> altExpNames(0):
```

Let’s look at the coefficients

``` r
round(fit$diffemb_coefficients, 5)
#> , , Intercept
#> 
#>           [,1]     [,2]
#>  [1,]  0.00562 -0.00553
#>  [2,]  0.01415 -0.01122
#>  [3,] -0.06522 -0.16661
#>  [4,]  1.10901  0.02092
#>  [5,]  0.33172  0.41483
#>  [6,] -0.57218 -0.83102
#>  [7,] -0.10213 -0.18001
#>  [8,]  0.14411 -0.13790
#>  [9,] -0.27461 -0.08411
#> [10,]  0.12748 -0.41115
#> 
#> , , conditionb
#> 
#>           [,1]     [,2]
#>  [1,] -0.00387  0.00878
#>  [2,] -0.00964  0.01997
#>  [3,]  0.31304  0.27604
#>  [4,] -1.37932  0.38055
#>  [5,] -0.43357 -0.25998
#>  [6,]  0.37321  1.24946
#>  [7,]  0.08710  0.08453
#>  [8,] -0.10357  0.37645
#>  [9,]  0.25211  0.18000
#> [10,] -0.14092  0.29485
plot(t(fit$diffemb_embedding), col = as.factor(col_data$condition))
```

<img src="man/figures/README-show-fit-results-1.png" width="100%" />

Bootstrap to get an estimate of the parameter variance

``` r
fit <- DiffEmbSeq::estimate_variance(fit, n_bootstrap_samples = 30)
#> Start bootstrap iteration 1
#> Start bootstrap iteration 2
#> Start bootstrap iteration 3
#> Start bootstrap iteration 4
#> Start bootstrap iteration 5
#> Start bootstrap iteration 6
#> Start bootstrap iteration 7
#> Start bootstrap iteration 8
#> Start bootstrap iteration 9
#> Start bootstrap iteration 10
#> Start bootstrap iteration 11
#> Start bootstrap iteration 12
#> Start bootstrap iteration 13
#> Start bootstrap iteration 14
#> Start bootstrap iteration 15
#> Start bootstrap iteration 16
#> Start bootstrap iteration 17
#> Start bootstrap iteration 18
#> Start bootstrap iteration 19
#> Start bootstrap iteration 20
#> Start bootstrap iteration 21
#> Start bootstrap iteration 22
#> Start bootstrap iteration 23
#> Start bootstrap iteration 24
#> Start bootstrap iteration 25
#> Start bootstrap iteration 26
#> Start bootstrap iteration 27
#> Start bootstrap iteration 28
#> Start bootstrap iteration 29
#> Start bootstrap iteration 30
```

``` r
res <- DiffEmbSeq::test_differential_expression(fit, contrast = fact(condition = "a") == fact(condition = "b"),
                                                consider = "embedding+linear", variance_est = "bootstrap", 
                                                return = "table")
head(res)
#>     feature   obs       pval  adj_pval       diff   adj_diff        sd
#> 1 feature_1 obs_1 0.42323380 0.9601503  0.4774037  0.8008236 0.5961410
#> 2 feature_2 obs_1 0.62604229 0.9994854 -0.2489013 -0.4873049 0.5107711
#> 3 feature_3 obs_1 0.34931102 0.9280894 -0.4225459 -0.9359265 0.4514734
#> 4 feature_4 obs_1 0.02360279 0.4168384 -0.9903750 -2.2635341 0.4375348
#> 5 feature_5 obs_1 0.50884418 0.9938363  0.1461087  0.6606385 0.2211629
#> 6 feature_6 obs_1 0.34926914 0.9280894 -0.2416553 -0.9360079 0.2581766
```

Show the gene expression changes on the latent embedding

``` r
pal <- scales::col_numeric(scales::viridis_pal()(50), NULL)
plot(t(fit$diffemb_embedding), col = pal(res[res$feature == "feature_10",]$diff),
     pch = 16, cex = 1.3)
```

<img src="man/figures/README-visualize-differential-expression-1.png" width="100%" />

# Session Info

``` r
sessionInfo()
#> R version 4.1.1 (2021-08-10)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.7                  DiffEmbSeq_0.0.5           
#>  [3] compiler_4.1.1              GenomeInfoDb_1.28.4        
#>  [5] highr_0.9                   XVector_0.32.0             
#>  [7] MatrixGenerics_1.9.0        bitops_1.0-7               
#>  [9] tools_4.1.1                 zlibbioc_1.38.0            
#> [11] SingleCellExperiment_1.14.1 digest_0.6.27              
#> [13] viridisLite_0.4.0           lifecycle_1.0.0            
#> [15] evaluate_0.14               lattice_0.20-44            
#> [17] rlang_1.0.3                 Matrix_1.3-4               
#> [19] DelayedArray_0.18.0         cli_3.3.0                  
#> [21] rstudioapi_0.13             yaml_2.2.1                 
#> [23] parallel_4.1.1              expm_0.999-6               
#> [25] xfun_0.26                   fastmap_1.1.0              
#> [27] GenomeInfoDbData_1.2.6      stringr_1.4.0              
#> [29] knitr_1.34                  S4Vectors_0.30.0           
#> [31] IRanges_2.26.0              stats4_4.1.1               
#> [33] grid_4.1.1                  Biobase_2.52.0             
#> [35] R6_2.5.1                    rmarkdown_2.11             
#> [37] irlba_2.3.3                 farver_2.1.0               
#> [39] magrittr_2.0.1              scales_1.1.1               
#> [41] matrixStats_0.60.1          htmltools_0.5.2            
#> [43] BiocGenerics_0.38.0         GenomicRanges_1.44.0       
#> [45] SummarizedExperiment_1.22.0 colorspace_2.0-2           
#> [47] stringi_1.7.4               glmGamPoi_1.11.0           
#> [49] munsell_0.5.0               RCurl_1.98-1.4
```

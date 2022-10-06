
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
valid until the 05.11.2022.

``` r
Sys.setenv(GITLAB_PAT = "glpat-KvTnh36v4jitb_yQaEch")
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
#> Find base point for differential embedding
#> Fit differential embedding model
#> -Iteration: 0    error: 6.47e+03
#> ---Fit Grassmann linear model
#> ---Update linear regression
#> -Iteration: 1    error: 4.85e+03
#> ---Fit Grassmann linear model
#> ---Update linear regression
#> -Iteration: 2    error: 4.85e+03
#> Converged
fit
#> class: DiffEmbFit 
#> dim: 30 500 
#> metadata(10): n_ambient n_embedding ... alignment_coefficients ''
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
#>  [1,]  0.00000  0.00000
#>  [2,]  0.00000  0.00000
#>  [3,] -0.51930  0.61979
#>  [4,] -0.32508  0.06652
#>  [5,]  0.20514 -0.29096
#>  [6,]  0.22391  0.01630
#>  [7,] -0.24185 -0.07071
#>  [8,]  0.17054 -0.16438
#>  [9,] -0.07646 -0.06271
#> [10,] -0.03062  0.00817
#> 
#> , , conditionb
#> 
#>           [,1]     [,2]
#>  [1,]  0.00000  0.00000
#>  [2,]  0.00000  0.00000
#>  [3,]  0.58501 -1.38694
#>  [4,]  0.60686  0.51942
#>  [5,] -0.20612  0.10382
#>  [6,] -0.54786 -0.03385
#>  [7,]  0.34128  0.27822
#>  [8,] -0.16600  0.07207
#>  [9,]  0.11197  0.20350
#> [10,] -0.03057  0.01884
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
#>     feature   obs      pval  adj_pval        diff   adj_diff        sd
#> 1 feature_1 obs_1 0.8689930 0.9956089 -0.04006239 -0.1649377 0.2428940
#> 2 feature_2 obs_1 0.5717281 0.9750561  0.17651967  0.5655083 0.3121434
#> 3 feature_3 obs_1 0.6291247 0.9837155 -0.16292456 -0.4829591 0.3373465
#> 4 feature_4 obs_1 0.7890687 0.9899962  0.06962060  0.2675201 0.2602444
#> 5 feature_5 obs_1 0.6699281 0.9878032  0.08134729  0.4262467 0.1908456
#> 6 feature_6 obs_1 0.8437783 0.9931445 -0.04825848 -0.1970630 0.2448886
```

Show the gene expression changes on the latent embedding

``` r
pal <- scales::col_numeric(scales::viridis_pal()(50), NULL)
plot(t(fit$diffemb_embedding), col = pal(res[res$feature == "feature_10",]$diff),
     pch = 16, cex = 1.3)
```

<img src="man/figures/README-visualize-differential-expression-1.png" width="100%" />

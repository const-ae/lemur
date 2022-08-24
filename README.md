
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
valid until the 31.08.2022.

``` r
Sys.setenv(GITLAB_PAT = "glpat-_2uHibDPXsgcs9kYxMsM")
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
#> Fit Grassmann linear model
#> Align points
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
#>  [1,]  0.00201 -0.00393
#>  [2,] -0.00119  0.00207
#>  [3,]  0.06324  0.33867
#>  [4,] -0.23304 -0.46719
#>  [5,]  0.08602  0.46378
#>  [6,] -0.13781  0.09502
#>  [7,] -0.09884 -0.14099
#>  [8,] -0.13983  0.24415
#>  [9,]  0.15734 -0.05495
#> [10,] -0.03631 -0.03797
#> 
#> , , conditionb
#> 
#>           [,1]     [,2]
#>  [1,] -0.01751  0.01812
#>  [2,]  0.01143 -0.01220
#>  [3,]  0.19583 -1.18850
#>  [4,]  1.05134  0.32174
#>  [5,]  0.13315 -0.93826
#>  [6,]  0.78204 -0.32090
#>  [7,] -0.01098  0.25181
#>  [8,] -0.15715  0.20281
#>  [9,] -0.16325 -0.04712
#> [10,] -0.22024  0.10941
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
#>     feature   obs       pval  adj_pval          diff     adj_diff        sd
#> 1 feature_1 obs_1 0.99885395 0.9998219  0.0003341028  0.001436355 0.2326046
#> 2 feature_2 obs_1 0.33439082 0.9998219 -0.2372139626 -0.965307483 0.2457393
#> 3 feature_3 obs_1 0.07314758 0.9998219  0.5001220937  1.791908762 0.2791002
#> 4 feature_4 obs_1 0.44902426 0.9998219  0.2108653023  0.757042729 0.2785382
#> 5 feature_5 obs_1 0.63679548 0.9998219 -0.0851523931 -0.472183976 0.1803373
#> 6 feature_6 obs_1 0.94997553 0.9998219 -0.0183343559 -0.062737513 0.2922391
```

Show the gene expression changes on the latent embedding

``` r
pal <- scales::col_numeric(scales::viridis_pal()(50), NULL)
plot(t(fit$diffemb_embedding), col = pal(res[res$feature == "feature_10",]$diff),
     pch = 16, cex = 1.3)
```

<img src="man/figures/README-visualize-differential-expression-1.png" width="100%" />


<!-- README.md is generated from README.Rmd. Please edit that file -->

# Latent Embedding Multivariate Regression (LEMUR)

<!-- badges: start -->
<!-- badges: end -->

The goal of lemur is to enable easy analysis of multi-condition
single-cell data. Lemur fits latent embedding regression, which means it
tries to find the PCA embedding for each condition and parameterizes the
transition from one embedding to another. For this task, lemur uses
geodesic regression on the Grassmann manifold, which is solved
efficiently using tangent-space linear modeling. The result is an
interpretable model of the gene expression for arbitrary experimental
designs that can be expressed using a design matrix.

![Schematic of the matrix decomposition at the core of
LEMUR](man/figures/equation_schematic.png)

## Installation

Lemur depends on features from
[`glmGamPoi`](https://github.com/const-ae/glmGamPoi) which are only
available in the development version, so please install `glmGamPoi` from
GitHub before proceeding with lemur’s installation:

``` r
devtools::install_github("const-ae/glmGamPoi")
```

You can install the released version of lemur from Github:

``` r
devtools::install_github("const-ae/lemur")
```

## A word of caution

This package is being actively developed, and I am still making breaking
changes. I would be delighted if you decide to try out the package, and
please open an issue if you think you found a bug, have an idea for a
cool feature, or have any questions about how LEMUR works. Consider this
a *beta* release with the goal to gather feedback but be aware that code
written against the current version of lemur might not work in the
future.

## Quickstart

``` r
library(lemur)
library(SingleCellExperiment)

fit <- lemur(sce, design = ~ patient_id + condition, n_embedding = 15)
fit <- align_harmony(fit)   # This step is optional
fit <- test_de(fit, contrast = cond(condition = "ctrl") - cond(condition = "panobinostat"))
nei <- find_de_neighborhoods(fit, independent_matrix = counts(sce), group_by = vars(patient_id, condition))
```

## Example

We will demonstrate `lemur` using a dataset by Zhao et al. (2021). The
data consists of tumor biopsies from five glioblastomas which were
treated using panobinostat or used as a control. Accordingly, we will
analyze ten samples (patient-treatment combinations) using a paired
experimental design.

We start by loading some packages which are necessary to analyze the
data:

``` r
library(tidyverse)
library(SingleCellExperiment)
library(lemur)
set.seed(1)
```

We use a reduced-size version of the glioblastoma data that ships with
the `lemur` package.

``` r
lemur::glioblastoma_example_data
#> class: SingleCellExperiment 
#> dim: 300 5000 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(300): ENSG00000210082 ENSG00000118785 ... ENSG00000167468
#>   ENSG00000139289
#> rowData names(6): gene_id symbol ... strand. source
#> colnames(5000): CGCCAGAGCGCA AGCTTTACTGCG ... TGAACAGTGCGT TGACCGGAATGC
#> colData names(10): patient_id treatment_id ... sample_id id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Initially, the data separates by the known covariates `patient_id` and
`condition`.

``` r
orig_umap <- uwot::umap(as.matrix(t(logcounts(glioblastoma_example_data))))

as_tibble(orig_umap) %>%
  bind_cols(as_tibble(colData(glioblastoma_example_data))) %>%
  ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(color = patient_id, shape = condition), size = 0.5) +
    labs(title = "UMAP of logcounts")
#> Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
#> `.name_repair` is omitted as of tibble 2.0.0.
#> ℹ Using compatibility `.name_repair`.
```

<img src="man/figures/README-raw_umap-1.png" width="100%" />

We fit the LEMUR model by calling `lemur()`. We provide the experimental
design using a formula. The elements of the formula can refer to columns
of the `colData` of the `SingleCellExperiment` object. We also set the
number of latent dimensions, which has a similar interpretation as the
number of dimensions in PCA. Optionally, we can further align
corresponding cells using manually annotated cell types
(`align_by_grouping`) or an automated alignment procedure (e.g.,
`align_harmony`, `align_neighbors`).

``` r
fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition, n_embedding = 15, verbose = FALSE)

# We can regularize the alignment either using ridge regression
# or by allowing only rotations or stretching
fit <- align_harmony(fit, stretching = FALSE)
#> Select cells that are considered close with 'harmony'

fit
#> class: lemur_fit 
#> dim: 300 5000 
#> metadata(10): n_embedding design ... alignment_design_matrix row_mask
#> assays(1): expr
#> rownames(300): ENSG00000210082 ENSG00000118785 ... ENSG00000167468
#>   ENSG00000139289
#> rowData names(6): gene_id symbol ... strand. source
#> colnames(5000): CGCCAGAGCGCA AGCTTTACTGCG ... TGAACAGTGCGT TGACCGGAATGC
#> colData names(10): patient_id treatment_id ... sample_id id
#> reducedDimNames(2): linearFit embedding
#> mainExpName: NULL
#> altExpNames(0):
```

The `lemur()` function returns an object that extends
`SingleCellExperiment` and thus supports subsetting and all the usual
data accessor methods (e.g., `nrow`, `assay`, `colData`, `rowData`). In
addition, `lemur` overloads the `$` operator to allow easy access to
additional fields produced by the LEMUR model. For example, the
low-dimensional embedding can be accessed using `fit$embedding`:

``` r
umap <- uwot::umap(t(fit$embedding))

as_tibble(umap) %>%
  bind_cols(as_tibble(fit$colData)) %>%
  ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(color = patient_id, shape = condition), size = 0.5) +
    labs(title = "UMAP of latent space from LEMUR")
```

<img src="man/figures/README-lemur_umap-1.png" width="100%" />

The `test_de` function takes a `lemur_fit_obj` and returns with a new
assay `"DE"` with the predicted difference between two conditions
specified in the `contrast`. Note that `lemur` implements a special
notation for contrasts. Instead of providing a contrast vector or design
matrix column names, you provide for each *condition* the levels, and
`lemur` automatically forms the contrast vector. This makes the contrast
more readable.

``` r
fit <- test_de(fit, contrast = cond(condition = "panobinostat") - cond(condition = "ctrl"))
```

We can pick any gene and show the differential expression pattern on the
UMAP plot:

``` r
# EEF1A1
gene_sel <- "ENSG00000156508"

as_tibble(umap) %>%
  mutate(expr = assay(fit, "DE")[gene_sel,]) %>%
  ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(color = expr)) +
    scale_color_gradient2() +
    labs(title = "Differential expression on UMAP plot")
```

<img src="man/figures/README-umap_de-1.png" width="100%" />

Alternatively, we can use the matrix of differential expression values
(`assay(fit, "DE")`) to guide the selection of cell neighborhoods that
show consistent differential expression. If we provide a count matrix,
the function uses a pseudobulked differential expression test to confirm
the gene expression differences on the count level.

``` r
neighborhoods <- find_de_neighborhoods(fit, independent_matrix = counts(glioblastoma_example_data),
                                      group_by = vars(patient_id, condition),
                                      include_complement = FALSE, verbose = FALSE)

as_tibble(neighborhoods) %>%
  arrange(pval) %>%
  left_join(as_tibble(rowData(fit)), by = c("name" = "gene_id"))
#> # A tibble: 300 × 16
#>    name      selec…¹ indices n_cells sel_s…²    pval adj_p…³ f_sta…⁴   df1   df2
#>    <chr>     <lgl>   <I<lis>   <int>   <dbl>   <dbl>   <dbl>   <dbl> <int> <dbl>
#>  1 ENSG0000… TRUE    <int>      4769    90.6 1.16e-5 0.00348   124.      1  6.91
#>  2 ENSG0000… TRUE    <int>      3887   202.  7.06e-5 0.0106     70.8     1  6.91
#>  3 ENSG0000… TRUE    <int>      3918   -95.5 2.19e-4 0.0189     49.3     1  6.91
#>  4 ENSG0000… TRUE    <int>      3560   136.  3.16e-4 0.0189     43.8     1  6.91
#>  5 ENSG0000… TRUE    <int>      1633    69.2 3.90e-4 0.0189     40.8     1  6.91
#>  6 ENSG0000… TRUE    <int>      3888  -150.  3.94e-4 0.0189     40.7     1  6.91
#>  7 ENSG0000… TRUE    <int>      2900  -216.  4.42e-4 0.0189     39.2     1  6.91
#>  8 ENSG0000… TRUE    <int>      4775   -66.3 5.37e-4 0.0201     36.7     1  6.91
#>  9 ENSG0000… TRUE    <int>      4240    97.0 7.17e-4 0.0220     33.3     1  6.91
#> 10 ENSG0000… TRUE    <int>      2578    53.7 7.34e-4 0.0220     33.0     1  6.91
#> # … with 290 more rows, 6 more variables: lfc <dbl>, symbol <chr>,
#> #   chromosome <fct>, gene_length <int>, strand. <fct>, source <fct>, and
#> #   abbreviated variable names ¹​selection, ²​sel_statistic, ³​adj_pval,
#> #   ⁴​f_statistic
```

We can now specifically select regions with significant differential
expression:

``` r
# HLA-DRB1
sel_gene <- "ENSG00000196126"

as_tibble(umap) %>%
  mutate(expr = assay(fit, "DE")[sel_gene,]) %>%
  ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(color = expr)) +
    scale_color_gradient2() +
    labs(title = "Differential expression on UMAP plot")
```

<img src="man/figures/README-umap_de2-1.png" width="100%" />

To plot the boundaries of the differential expression neighborhood, we
create a helper dataframe and use the `geom_density2d` function from
`ggplot2`:

``` r
neighborhood_coordinates <- neighborhoods %>%
  filter(name == sel_gene) %>%
  mutate(cell_id = map(indices, \(idx) colnames(fit)[idx])) %>%
  unnest(c(indices, cell_id)) %>%
  left_join(as_tibble(umap, rownames = "cell_id"), by = "cell_id") %>%
  dplyr::select(name, cell_id, V1, V2)

as_tibble(umap) %>%
  mutate(expr = assay(fit, "DE")[sel_gene,]) %>%
  ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(color = expr)) +
    scale_color_gradient2() +
    geom_density2d(data = neighborhood_coordinates, breaks = 0.1, 
                   contour_var = "ndensity", color = "black") +
    labs(title = "Differential expression with neighborhood boundary")
```

<img src="man/figures/README-umap_de3-1.png" width="100%" />

We make a volcano plot of the differential expression results to better
understand the expression differences across all genes.

``` r
neighborhoods %>%
  ggplot(aes(x = lfc, y = -log10(pval))) +
    geom_point(aes(color  = adj_pval < 0.1)) +
    labs(title = "Vulcano plot of the neighborhoods")
```

<img src="man/figures/README-volcano_plot-1.png" width="100%" />

``` r
neighborhoods %>%
  ggplot(aes(x = n_cells, y = -log10(pval))) +
    geom_point(aes(color  = adj_pval < 0.1)) +
    labs(title = "Neighborhood size vs neighborhood significance")
```

<img src="man/figures/README-volcano_plot-2.png" width="100%" />

# FAQ

> I have already integrated my data using Harmony / MNN / Seurat. Can I
> call `lemur` directly with the aligned data?

No. You need to call `lemur` with the unaligned data so that it can
learn how much the expression of each gene changes between conditions.
What you can do, though, is to use the aligned data as a template for
lemur’s alignment step (`align_template`).

> Can I call lemur with
> [sctransformed](https://github.com/satijalab/sctransform) instead of
> log-transformed data?

Yes. You can call lemur with any variance stabilized count matrix. Based
on a [previous
project](https://www.biorxiv.org/content/10.1101/2021.06.24.449781v4), I
recommend to use log-transformation, but other methods will work just
fine.

> My data appears less integrated after calling `lemur()` than before.
> What is happening?!

This is a known issue and can be caused if the data has large
compositional shifts (for example, if one cell type disappears). The
problem is that the initial linear regression step, which centers the
conditions relative to each other, overcorrects and introduces a
consistent shift in the latent space. I am currently trying to develop a
generic solution for this problem, but until this is implemented, you
can manually fix the regression coefficients:

``` r
fit <- lemur(sce, design = ~ patient_id + condition, n_embedding = 15, linear_coefficient_estimator = "zero")
```

> The conditions still separate if I plot the data using UMAP / tSNE.
> Even after calling `align_harmony` / `align_neighbors`. What should I
> do?

You can try to increase `n_embedding`. If this still does not help,
there is little use in inferring differential expression neighborhoods.
But as I haven’t encountered such a dataset yet, I would like to try it
out myself. If you can share the data publicly, please open an issue.

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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] lubridate_1.9.2             forcats_1.0.0              
#>  [3] stringr_1.5.0               dplyr_1.1.0                
#>  [5] purrr_1.0.1                 readr_2.1.4                
#>  [7] tidyr_1.3.0                 tibble_3.1.8               
#>  [9] ggplot2_3.4.1               tidyverse_2.0.0            
#> [11] SingleCellExperiment_1.20.0 SummarizedExperiment_1.28.0
#> [13] Biobase_2.58.0              GenomicRanges_1.50.2       
#> [15] GenomeInfoDb_1.34.9         IRanges_2.32.0             
#> [17] S4Vectors_0.36.2            BiocGenerics_0.44.0        
#> [19] MatrixGenerics_1.10.0       matrixStats_0.63.0         
#> [21] lemur_0.0.12               
#> 
#> loaded via a namespace (and not attached):
#>  [1] splines_4.2.1             DelayedMatrixStats_1.20.0
#>  [3] expm_0.999-7              highr_0.10               
#>  [5] GenomeInfoDbData_1.2.9    yaml_2.3.7               
#>  [7] pillar_1.8.1              lattice_0.20-45          
#>  [9] glue_1.6.2                digest_0.6.31            
#> [11] XVector_0.38.0            colorspace_2.1-0         
#> [13] cowplot_1.1.1             htmltools_0.5.4          
#> [15] Matrix_1.5-3              pkgconfig_2.0.3          
#> [17] zlibbioc_1.44.0           scales_1.2.1             
#> [19] tzdb_0.3.0                timechange_0.2.0         
#> [21] generics_0.1.3            farver_2.1.1             
#> [23] ellipsis_0.3.2            withr_2.5.0              
#> [25] harmony_0.1.1             cli_3.6.0                
#> [27] magrittr_2.0.3            evaluate_0.20            
#> [29] fansi_1.0.4               MASS_7.3-58.2            
#> [31] tools_4.2.1               hms_1.1.2                
#> [33] lifecycle_1.0.3           munsell_0.5.0            
#> [35] DelayedArray_0.24.0       irlba_2.3.5.1            
#> [37] isoband_0.2.7             compiler_4.2.1           
#> [39] rlang_1.0.6               grid_4.2.1               
#> [41] RCurl_1.98-1.10           rstudioapi_0.14          
#> [43] RcppAnnoy_0.0.20          glmGamPoi_1.11.7         
#> [45] bitops_1.0-7              labeling_0.4.2           
#> [47] rmarkdown_2.20            gtable_0.3.1             
#> [49] codetools_0.2-19          R6_2.5.1                 
#> [51] knitr_1.42                fastmap_1.1.1            
#> [53] uwot_0.1.14               utf8_1.2.3               
#> [55] stringi_1.7.12            Rcpp_1.0.10              
#> [57] vctrs_0.5.2               tidyselect_1.2.0         
#> [59] xfun_0.37                 sparseMatrixStats_1.10.0
```

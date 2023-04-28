
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
set.seed(42)
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
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
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

# We can regularize the alignment using ridge penalty
fit <- align_harmony(fit, ridge_penalty = 0.5)
#> Select cells that are considered close with 'harmony'

fit
#> class: lemur_fit 
#> dim: 300 5000 
#> metadata(9): n_embedding design ... alignment_design_matrix row_mask
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
sel_gene <- "ENSG00000156508"

as_tibble(umap) %>%
  mutate(expr = assay(fit, "DE")[sel_gene,]) %>%
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
  left_join(as_tibble(rowData(fit)), by = c("name" = "gene_id")) %>%
  dplyr::select(name, symbol, n_cells, pval, adj_pval)
#> # A tibble: 300 × 5
#>    name            symbol n_cells      pval adj_pval
#>    <chr>           <chr>    <int>     <dbl>    <dbl>
#>  1 ENSG00000245532 NEAT1     3831 0.0000152  0.00244
#>  2 ENSG00000187193 MT1X      4032 0.0000163  0.00244
#>  3 ENSG00000125148 MT2A      3293 0.0000675  0.00586
#>  4 ENSG00000147588 PMP2      3683 0.0000782  0.00586
#>  5 ENSG00000169715 MT1E      3479 0.000287   0.0142 
#>  6 ENSG00000156508 EEF1A1    2604 0.000327   0.0142 
#>  7 ENSG00000177700 POLR2L    3466 0.000334   0.0142 
#>  8 ENSG00000198668 CALM1     4087 0.000379   0.0142 
#>  9 ENSG00000175899 A2M       3613 0.000563   0.0186 
#> 10 ENSG00000069275 NUCKS1    3851 0.000672   0.0186 
#> # ℹ 290 more rows
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
  dplyr::filter(name == sel_gene) %>%
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
#> R version 4.3.0 (2023-04-21)
#> Platform: x86_64-apple-darwin20 (64-bit)
#> Running under: macOS Big Sur 11.7
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Berlin
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] lubridate_1.9.2             forcats_1.0.0              
#>  [3] stringr_1.5.0               dplyr_1.1.2                
#>  [5] purrr_1.0.1                 readr_2.1.4                
#>  [7] tidyr_1.3.0                 tibble_3.2.1               
#>  [9] ggplot2_3.4.2               tidyverse_2.0.0            
#> [11] SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.0
#> [13] Biobase_2.60.0              GenomicRanges_1.52.0       
#> [15] GenomeInfoDb_1.36.0         IRanges_2.34.0             
#> [17] S4Vectors_0.38.0            BiocGenerics_0.46.0        
#> [19] MatrixGenerics_1.12.0       matrixStats_0.63.0         
#> [21] lemur_0.0.15               
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.3              xfun_0.39                
#>  [3] lattice_0.21-8            tzdb_0.3.0               
#>  [5] vctrs_0.6.2               tools_4.3.0              
#>  [7] bitops_1.0-7              generics_0.1.3           
#>  [9] fansi_1.0.4               highr_0.10               
#> [11] pkgconfig_2.0.3           Matrix_1.5-4             
#> [13] sparseMatrixStats_1.12.0  lifecycle_1.0.3          
#> [15] GenomeInfoDbData_1.2.10   farver_2.1.1             
#> [17] compiler_4.3.0            munsell_0.5.0            
#> [19] codetools_0.2-19          glmGamPoi_1.12.0         
#> [21] htmltools_0.5.5           RCurl_1.98-1.12          
#> [23] yaml_2.3.7                crayon_1.5.2             
#> [25] pillar_1.9.0              MASS_7.3-59              
#> [27] uwot_0.1.14               DelayedArray_0.25.0      
#> [29] tidyselect_1.2.0          digest_0.6.31            
#> [31] stringi_1.7.12            labeling_0.4.2           
#> [33] splines_4.3.0             cowplot_1.1.1            
#> [35] fastmap_1.1.1             grid_4.3.0               
#> [37] colorspace_2.1-0          cli_3.6.1                
#> [39] harmony_0.1.1             magrittr_2.0.3           
#> [41] utf8_1.2.3                withr_2.5.0              
#> [43] DelayedMatrixStats_1.22.0 scales_1.2.1             
#> [45] timechange_0.2.0          rmarkdown_2.21           
#> [47] XVector_0.40.0            hms_1.1.3                
#> [49] evaluate_0.20             knitr_1.42               
#> [51] RcppAnnoy_0.0.20          irlba_2.3.5.1            
#> [53] rlang_1.1.0               isoband_0.2.7            
#> [55] Rcpp_1.0.10               glue_1.6.2               
#> [57] rstudioapi_0.14           R6_2.5.1                 
#> [59] zlibbioc_1.46.0
```

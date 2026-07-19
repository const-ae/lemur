## Regenerates the golden test fixtures for R/max_diversity_clustering.R,
## used by tests/testthat/test-max_diversity_clustering.R to validate the
## pure-R re-implementation against the real `harmony` package, version 1.2.4
## (the last version before upstream's breaking changes).
##
## The fixture .rds files are NOT committed to the repo (see .gitignore) -
## they are sizable binary files, and CI/autobuild systems don't have (and
## shouldn't need) the patched harmony build this script depends on. The
## corresponding tests skip gracefully when the fixtures are missing. Run
## this script manually (`Rscript data-raw/generate_max_diversity_clustering_fixtures.R`
## from the package root) whenever the fixtures need to be (re-)created, and
## commit nothing but this script.
##
## This script is self-contained: it downloads harmony 1.2.4's exact source
## from the CRAN archive (the last version before upstream's breaking
## changes), patches in a test-only hook, and builds it into an isolated
## library - no local harmony checkout needed. It DOES need internet access
## and a working C++ toolchain (the same one used to build lemur itself).
##
## The patch adds an `arma::uvec fixed_update_order` field to the C++ harmony
## object: when set to a length-N permutation, update_R() uses it verbatim as
## the cell processing order instead of calling Armadillo's shuffle(). Real
## harmony's block-shuffle uses Armadillo's own RNG, which R's set.seed()
## can't reach - injecting a fixed, shared permutation is what makes the
## iterative clustering loop (not just the deterministic k-means init step)
## bit-exactly comparable between harmony and lemur's pure-R port.
##
## For the real-data fixture, this script also needs the muscData,
## ExperimentHub, irlba, and sparseMatrixStats packages (all Suggests-only,
## dev-time tooling - not added to lemur's own DESCRIPTION) and internet
## access to download/cache the Kang et al. 2018 dataset.

harmony_version <- "1.2.4"
harmony_url <- sprintf("https://cran.r-project.org/src/contrib/Archive/harmony/harmony_%s.tar.gz", harmony_version)
harmony_lib <- file.path(tempdir(), "lemur-harmony-fixture-lib")
harmony_src_dir <- file.path(tempdir(), "lemur-harmony-fixture-src")
out_dir <- "tests/testthat/fixtures/max_diversity_clustering"

stopifnot(
  "Run this script from the lemur package root (Rscript data-raw/generate_max_diversity_clustering_fixtures.R)" =
    file.exists("DESCRIPTION") && file.exists("R/max_diversity_clustering.R")
)

## ---- 1. Download harmony 1.2.4's exact source from the CRAN archive ------
unlink(harmony_src_dir, recursive = TRUE)
dir.create(harmony_src_dir, recursive = TRUE)
tarball <- file.path(harmony_src_dir, "harmony.tar.gz")
utils::download.file(harmony_url, tarball, quiet = TRUE, mode = "wb")
utils::untar(tarball, exdir = harmony_src_dir)
harmony_pkg_dir <- file.path(harmony_src_dir, "harmony")
stopifnot(
  "Download/extraction of the harmony source failed" = dir.exists(harmony_pkg_dir),
  "Downloaded harmony source is not version 1.2.4" =
    any(grepl(paste0("^Version:\\s*", harmony_version, "$"),
              readLines(file.path(harmony_pkg_dir, "DESCRIPTION"))))
)

## ---- 2. Patch in the fixed_update_order test-only hook -------------------
patch_file <- function(path, replacements){
  txt <- paste(readLines(path), collapse = "\n")
  for(r in replacements){
    stopifnot("patch target text not found - has the harmony source changed?" = grepl(r$old, txt, fixed = TRUE))
    txt <- sub(r$old, r$new, txt, fixed = TRUE)
  }
  writeLines(txt, path)
}

patch_file(file.path(harmony_pkg_dir, "src", "harmony.h"), list(
  list(
    old = "  VECTYPE Pr_b, theta, N_b, sigma, lambda;",
    new = paste(
      "  VECTYPE Pr_b, theta, N_b, sigma, lambda;",
      "",
      "  // TEST-ONLY HOOK (not part of upstream harmony): when this has length N,",
      "  // update_R() uses it verbatim as the cell processing order instead of",
      "  // calling shuffle(). Used to generate reproducible golden fixtures for",
      "  // lemur's pure-R port, where Armadillo's RNG can't be matched from R.",
      "  arma::uvec fixed_update_order;",
      sep = "\n"
    )
  )
))

patch_file(file.path(harmony_pkg_dir, "src", "harmony.cpp"), list(
  list(
    old = paste(
      "  // Generate the 0,N-1 indices",
      "  uvec indices = linspace<uvec>(0, N - 1, N);",
      "  update_order = shuffle(indices);",
      sep = "\n"
    ),
    new = paste(
      "  // Generate the 0,N-1 indices",
      "  uvec indices = linspace<uvec>(0, N - 1, N);",
      "  if (fixed_update_order.n_elem == N) {",
      "    // TEST-ONLY HOOK: see declaration in harmony.h",
      "    update_order = fixed_update_order;",
      "  } else {",
      "    update_order = shuffle(indices);",
      "  }",
      sep = "\n"
    )
  ),
  list(
    old = paste(
      "      .field(\"B_vec\", &harmony::B_vec)",
      "      .field(\"alpha\", &harmony::alpha)",
      "      ;",
      sep = "\n"
    ),
    new = paste(
      "      .field(\"B_vec\", &harmony::B_vec)",
      "      .field(\"alpha\", &harmony::alpha)",
      "      .field(\"fixed_update_order\", &harmony::fixed_update_order)",
      "      ;",
      sep = "\n"
    )
  )
))

## ---- 3. Build the patched harmony 1.2.4 into an isolated library ---------
dir.create(harmony_lib, showWarnings = FALSE, recursive = TRUE)
install_out <- system2("R", c("CMD", "INSTALL", "--no-multiarch",
                               paste0("--library=", shQuote(harmony_lib)), shQuote(harmony_pkg_dir)),
                       stdout = TRUE, stderr = TRUE)
cat(install_out, sep = "\n")
if(! dir.exists(file.path(harmony_lib, "harmony"))) stop("harmony installation failed, see output above.")

library(harmony, lib.loc = harmony_lib)
stopifnot(as.character(packageVersion("harmony", lib.loc = harmony_lib)) == "1.2.4")
library(Matrix)
library(vctrs)

get_groups <- function(design_matrix) vctrs::vec_group_id(unclass(design_matrix))

## Faithful copy of lemur's (removed) R/harmony_wrapper.R::harmony_init,
## pointed at the real harmony package, used only to generate these fixtures.
## Constructs the Rcpp module object directly (methods::new("Rcpp_harmony"))
## rather than via harmony::RunHarmony(..., max.iter = 0, return_object = TRUE)
## the way lemur's old wrapper did: RunHarmony's real parameter is `max_iter`
## (underscore), so `max.iter = 0` was silently swallowed as an unrecognized
## `...` arg and a full, wasted 10-round Harmony run happened before being
## immediately overwritten. Bypassing it gives a clean, side-effect-free
## reference for the algorithm itself.
harmony_init_ref <- function(embedding, design_matrix,
                              theta = 2, lambda = 1, sigma = 0.1,
                              nclust = min(round(ncol(embedding) / 30), 100),
                              tau = 0, block.size = 0.1, max.iter.cluster = 200,
                              epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, verbose = FALSE){
  mm_groups <- get_groups(design_matrix)
  n_groups <- length(unique(mm_groups))
  phi <- matrix(0, nrow = n_groups, ncol = ncol(embedding))
  phi[mm_groups + n_groups * (seq_along(mm_groups)-1)] <- 1
  phi <- as(phi, "sparseMatrix")

  N_b <- MatrixGenerics::rowSums2(phi)
  theta <- rep_len(theta, n_groups)
  theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))
  lambda_vec <- c(0, rep_len(lambda, n_groups))
  sigma <- rep_len(sigma, nclust)

  harmonyObj <- methods::new("Rcpp_harmony")
  harmonyObj$setup(embedding, phi, sigma, theta, lambda_vec, 0.2, max.iter.cluster, epsilon.cluster,
                    epsilon.harmony, nclust, block.size, rep(1, n_groups), verbose)
  ## NOTE: lemur's old wrapper additionally did
  ##   harmonyObj$Y <- t(stats::kmeans(t(harmonyObj$Z_cos), centers=K, ...)$centers)
  ## right before this call - but init_cluster_cpp() (C++) immediately
  ## overwrites Y with the result of its OWN internal kmeans_centers() call
  ## (a second, separate call to stats::kmeans via Rcpp). That manual R-side
  ## call is fully discarded except that it shifts the RNG stream consumed by
  ## the real (second) kmeans call - a redundant-computation artifact, not
  ## part of the clustering algorithm, so it's not replicated here.
  harmonyObj$init_cluster_cpp()
  harmonyObj
}

harmony_cluster_ref <- function(harmonyObj){
  err_status <- harmonyObj$cluster_cpp()
  stopifnot(err_status == 0)
  harmonyObj
}

## Z_cos/Z_orig/N/B are omitted: they're fully determined by the top-level
## `embedding`/`design_matrix` saved alongside each snapshot.
snapshot <- function(obj){
  list(
    Y = obj$Y, R = obj$R, E = obj$E, O = obj$O,
    K = obj$K,
    theta = obj$theta, sigma = obj$sigma, Pr_b = obj$Pr_b,
    objective_kmeans = obj$objective_kmeans,
    kmeans_rounds = obj$kmeans_rounds
  )
}

make_embedding <- function(n_emb, n_cells, n_centers, seed, spread = 1){
  set.seed(seed)
  centers <- matrix(rnorm(n_emb * n_centers, sd = 3), nrow = n_emb)
  assignment <- sample.int(n_centers, n_cells, replace = TRUE)
  centers[, assignment, drop = FALSE] + matrix(rnorm(n_emb * n_cells, sd = spread), nrow = n_emb)
}

## Generates one fixture: init (k-means + first soft assignment - exact,
## deterministic, no block_size or shuffle involved) and final (after running
## cluster_cpp() to convergence with an injected fixed cell-processing order,
## so the result is bit-exactly reproducible from the pure-R port).
##
## block_size is always kept < 1: block_size = 1 (a single block covering all
## cells) is numerically degenerate for this algorithm - removing the only
## block leaves E = O = 0, a genuine 0/0 - so the finite value that happens to
## come out is floating-point-summation-order noise, not a meaningful
## cross-implementation target.
save_case <- function(name, embedding, design_matrix, ..., seed, block_size = 0.1){
  set.seed(seed)
  init <- suppressWarnings(harmony_init_ref(embedding, design_matrix, ..., block.size = block_size))
  init_snapshot <- snapshot(init)

  perm <- sample.int(ncol(embedding))
  init$fixed_update_order <- perm - 1L  # Armadillo uvec is 0-indexed

  clustered <- harmony_cluster_ref(init)
  final_snapshot <- snapshot(clustered)
  saveRDS(
    list(embedding = embedding, design_matrix = design_matrix,
         params = c(list(...), list(block_size = block_size)),
         seed = seed, fixed_update_order = perm,
         init = init_snapshot, final = final_snapshot),
    file.path(out_dir, paste0(name, ".rds"))
  )
  cat("saved", name, "- kmeans_rounds:", clustered$kmeans_rounds, "\n")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 2. Synthetic scenarios ------------------------------------------------
n_emb <- 5; n_cells <- 300
embedding <- make_embedding(n_emb, n_cells, n_centers = 4, seed = 100)
batch <- sample(c("a", "b", "c"), n_cells, replace = TRUE, prob = c(0.6, 0.3, 0.1))
design_matrix <- model.matrix(~ batch)
save_case("multi_batch_unbalanced", embedding, design_matrix, nclust = 10, theta = 2, sigma = 0.1, seed = 1)
save_case("nclust_1", embedding, design_matrix, nclust = 1, theta = 2, sigma = 0.1, seed = 3)
save_case("tau_scaling", embedding, design_matrix, nclust = 10, theta = 2, sigma = 0.1, tau = 5, seed = 4)

n_cells_big <- 900
embedding_big <- make_embedding(3, n_cells_big, n_centers = 5, seed = 300)
patient <- sample(paste0("p", 1:4), n_cells_big, replace = TRUE)
condition <- sample(c("ctrl", "trt"), n_cells_big, replace = TRUE)
design_matrix_big <- model.matrix(~ patient + condition)
save_case("realistic_defaults", embedding_big, design_matrix_big,
          nclust = min(round(n_cells_big / 30), 100), theta = 2, sigma = 0.1, seed = 6)

n_emb2 <- 4; n_cells2 <- 240
embedding2 <- make_embedding(n_emb2, n_cells2, n_centers = 3, seed = 500)
batch2 <- sample(c("x", "y"), n_cells2, replace = TRUE, prob = c(0.7, 0.3))
design_matrix2 <- model.matrix(~ batch2)
save_case("small_block_size", embedding2, design_matrix2, nclust = 8, theta = 2, sigma = 0.1,
          seed = 11, block_size = 0.05)

## ---- 3. Guard conditions (N < 6 stop, N < 40 block.size override) --------
## Note: setup() forces block.size regardless of what's requested when N < 40,
## so these can't be deterministic golden snapshots (Armadillo shuffle kicks
## in). We only record the condition text here; the pure-R port is tested for
## equivalent guard behavior directly, not by comparing cluster values.
guard_conditions <- list()
embedding_tiny <- make_embedding(n_emb, 5, n_centers = 2, seed = 400)
design_tiny <- model.matrix(~ rep(c("a", "b"), length.out = 5))
guard_conditions$too_few_cells <- tryCatch({
  harmony_init_ref(embedding_tiny, design_tiny, nclust = 2)
  "no error"
}, error = function(e) conditionMessage(e))

embedding_36 <- make_embedding(n_emb, 36, n_centers = 3, seed = 401)
design_36 <- model.matrix(~ sample(c("a", "b"), 36, replace = TRUE))
guard_conditions$block_size_override <- tryCatch({
  tryCatch({
    harmony_init_ref(embedding_36, design_36, nclust = 4, block.size = 0.9)
    "no warning"
  }, warning = function(w) conditionMessage(w))
}, error = function(e) conditionMessage(e))
saveRDS(guard_conditions, file.path(out_dir, "guard_conditions.rds"))
cat("saved guard_conditions\n")

## ---- 4. Real-data fixture (Kang et al. 2018) ------------------------------
real_data_deps <- c("muscData", "ExperimentHub", "irlba", "sparseMatrixStats")
missing_deps <- real_data_deps[! vapply(real_data_deps, requireNamespace, logical(1), quietly = TRUE)]
if(length(missing_deps) > 0){
  message("Skipping the real-data (Kang et al.) fixture - missing packages: ",
          paste(missing_deps, collapse = ", "))
}else{
  suppressPackageStartupMessages(library(muscData))
  sce <- Kang18_8vs8()
  counts <- assay(sce, "counts")

  set.seed(42)
  n_sub <- 6000
  sub_idx <- sample.int(ncol(sce), n_sub)
  counts <- counts[, sub_idx]
  col_data <- as.data.frame(colData(sce)[sub_idx, ])

  lib_size <- Matrix::colSums(counts)
  norm <- log1p(Matrix::t(Matrix::t(counts) / lib_size * median(lib_size)))
  gene_var <- sparseMatrixStats::rowVars(norm)
  hvg <- order(gene_var, decreasing = TRUE)[1:1000]
  norm_hvg <- as.matrix(norm[hvg, ])
  centered <- norm_hvg - rowMeans(norm_hvg)

  pca <- irlba::irlba(t(centered), nv = 15)
  kang_embedding <- t(pca$u %*% diag(pca$d))

  kang_design_matrix <- model.matrix(~ ind, data = col_data)
  save_case("kang_real_data", kang_embedding, kang_design_matrix,
            nclust = 20, theta = 2, sigma = 0.1, seed = 42, block_size = 0.05)
}

cat("Done. Fixtures written to '", out_dir, "'.\n", sep = "")

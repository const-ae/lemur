## Golden-fixture tests for the pure-R maximum diversity clustering
## implementation (R/max_diversity_clustering.R), which replaces the
## `harmony` package dependency.
##
## Fixtures were generated from the real `harmony` R package, version 1.2.4
## (the last version before the breaking upstream changes that motivated this
## reimplementation), via data-raw/generate_max_diversity_clustering_fixtures.R.
## The fixture .rds files are NOT committed (see .gitignore) - they're
## sizable binary files, and CI/autobuild systems don't have (and shouldn't
## need) that script's dependencies. Every test below skips gracefully if its
## fixture is missing; run the generation script locally to exercise them.
## See the fixture files' `params`/`seed`/`fixed_update_order` for exact
## reproduction. Six scenarios are covered:
##  - multi_batch_unbalanced: 3 batches with unbalanced sizes (60/30/10%)
##  - nclust_1: degenerate single-cluster edge case
##  - tau_scaling: tau != 0, exercises the theta rescaling formula
##  - realistic_defaults: 900 cells, 8 patient x condition groups, nclust
##    computed the way align_harmony() computes it - the longest chain (75
##    iterations), used with a looser tolerance (see below)
##  - small_block_size: block_size = 0.05, close to align_harmony()'s default
##  - kang_real_data: real single-cell data (Kang et al. 2018 IFN-beta
##    stimulation dataset, muscData::Kang18_8vs8), 6000-cell subsample,
##    PCA embedding on top-1000 HVGs, design = ~ind (8 patients)
##
## For each scenario we compare:
##  (a) the initialization step (k-means centroids + first soft assignment) -
##      this only depends on stats::kmeans(), so with a matching seed it's
##      exactly reproducible (no RNG-mismatch or numerical-degeneracy issues).
##  (b) the iterative clustering loop - real harmony's block-shuffle uses
##      Armadillo's (R-unreachable) RNG, so to get an exact reference the
##      generation script patches a `fixed_update_order` test-only hook into
##      a fresh harmony 1.2.4 build that lets update_R() use an externally
##      supplied permutation. The exact same permutation is injected here via
##      `state$fixed_update_order` (a hook that only exists for testing; normal
##      use always draws a fresh `sample.int()` per call).
##
## Note: block_size = 1 (a single block covering all cells) is deliberately
## never used as a "final" reference - it's numerically degenerate (removing
## the only block leaves E = O = 0, a genuine 0/0), so the finite value that
## happens to come out is floating-point-summation-order noise, not a
## meaningful cross-implementation target.

fixture_dir <- function(name) test_path("fixtures", "max_diversity_clustering", name)

load_fixture <- function(name){
  path <- fixture_dir(paste0(name, ".rds"))
  skip_if_not(file.exists(path),
              paste0("fixture '", name, ".rds' not found - run ",
                     "data-raw/generate_max_diversity_clustering_fixtures.R to create it"))
  readRDS(path)
}

expect_state_close <- function(state, snapshot, tolerance){
  expect_equal(state$Y, snapshot$Y, tolerance = tolerance, ignore_attr = TRUE)
  expect_equal(state$R, snapshot$R, tolerance = tolerance, ignore_attr = TRUE)
  expect_equal(state$E, snapshot$E, tolerance = tolerance, ignore_attr = TRUE)
  expect_equal(state$O, snapshot$O, tolerance = tolerance, ignore_attr = TRUE)
}

tight_cases <- c("multi_batch_unbalanced", "nclust_1", "tau_scaling", "small_block_size", "kang_real_data")

for(case in tight_cases){
  test_that(paste0("matches harmony 1.2.4 golden fixture: ", case), {
    f <- load_fixture(case)
    p <- f$params

    set.seed(f$seed)
    st <- init_max_diversity_clustering(f$embedding, f$design_matrix,
                                         theta = p$theta, sigma = p$sigma, nclust = p$nclust,
                                         tau = if(is.null(p$tau)) 0 else p$tau,
                                         block_size = p$block_size, verbose = FALSE)
    expect_state_close(st, f$init, tolerance = 1e-6)
    expect_equal(as.numeric(st$objective_kmeans), as.numeric(f$init$objective_kmeans), tolerance = 1e-3)

    st$fixed_update_order <- f$fixed_update_order
    st <- run_max_diversity_clustering(st)

    expect_state_close(st, f$final, tolerance = 1e-4)
    expect_equal(st$kmeans_rounds, f$final$kmeans_rounds)
    # f$final$objective_kmeans is harmony's *persistent* trace: it includes the
    # init-step value the R port does not carry over between calls (see
    # run_max_diversity_clustering()'s docs), hence the [-1].
    expect_equal(as.numeric(st$objective_kmeans), as.numeric(f$final$objective_kmeans[-1]), tolerance = 1e-2)
  })
}

test_that("matches harmony 1.2.4 golden fixture: realistic_defaults (long chain, looser tolerance)", {
  # This is the longest chain (75 harmony iterations). Tiny per-iteration
  # floating-point differences between R's and Armadillo's BLAS accumulate
  # over that many rounds and can shift the exact convergence iteration by
  # one - expected numerical drift, not an algorithmic mismatch (the other
  # five, shorter-chain fixtures all match to floating-point precision).
  f <- load_fixture("realistic_defaults")
  p <- f$params

  set.seed(f$seed)
  st <- init_max_diversity_clustering(f$embedding, f$design_matrix,
                                       theta = p$theta, sigma = p$sigma, nclust = p$nclust,
                                       block_size = p$block_size, verbose = FALSE)
  expect_state_close(st, f$init, tolerance = 1e-6)

  st$fixed_update_order <- f$fixed_update_order
  st <- run_max_diversity_clustering(st)

  expect_lte(abs(st$kmeans_rounds - f$final$kmeans_rounds), 2)
  expect_state_close(st, f$final, tolerance = 1e-2)
})


test_that("guard conditions match harmony 1.2.4 behavior", {
  guards <- load_fixture("guard_conditions")
  expect_match(guards$too_few_cells, "less than 6 cells")
  expect_match(guards$block_size_override, "Too few cells")

  design_matrix <- model.matrix(~ rep(c("a", "b"), length.out = 5))
  expect_error(
    init_max_diversity_clustering(randn(3, 5), design_matrix, nclust = 2, verbose = FALSE),
    "less than 6 cells"
  )

  design_matrix_36 <- model.matrix(~ sample(c("a", "b"), 36, replace = TRUE))
  expect_warning(
    st <- init_max_diversity_clustering(randn(3, 36), design_matrix_36, nclust = 4,
                                         block_size = 0.9, verbose = FALSE),
    "Too few cells"
  )
  expect_equal(st$block_size, 0.2)
})


test_that("init_max_diversity_clustering rejects a design with a single group", {
  design_matrix <- matrix(1, nrow = 50, ncol = 1)
  expect_error(
    init_max_diversity_clustering(randn(3, 50), design_matrix, nclust = 5, verbose = FALSE),
    "at least two groups"
  )
})


test_that("R is always a valid soft assignment (columns sum to 1, non-negative)", {
  set.seed(1)
  for(case in tight_cases){
    f <- load_fixture(case)
    p <- f$params
    st <- init_max_diversity_clustering(f$embedding, f$design_matrix,
                                         theta = p$theta, sigma = p$sigma, nclust = p$nclust,
                                         tau = if(is.null(p$tau)) 0 else p$tau,
                                         block_size = p$block_size, verbose = FALSE)
    expect_true(all(st$R >= 0))
    expect_equal(unname(colSums(st$R)), rep(1, ncol(st$R)), tolerance = 1e-8)

    st <- run_max_diversity_clustering(st)  # regular (random shuffle) path
    expect_true(all(st$R >= 0))
    expect_equal(unname(colSums(st$R)), rep(1, ncol(st$R)), tolerance = 1e-8)
  }
})


test_that("run_max_diversity_clustering stops at max_iter_cluster if it doesn't converge", {
  f <- load_fixture("multi_batch_unbalanced")
  p <- f$params
  set.seed(f$seed)
  st <- init_max_diversity_clustering(f$embedding, f$design_matrix,
                                       theta = p$theta, sigma = p$sigma, nclust = p$nclust,
                                       block_size = p$block_size, max_iter_cluster = 2, verbose = FALSE)
  st <- run_max_diversity_clustering(st)
  expect_equal(st$kmeans_rounds, 2)
  expect_equal(length(st$objective_kmeans), 2)
})


test_that("run_max_diversity_clustering can be re-entered with an updated embedding", {
  # Mirrors how align_harmony()'s outer loop re-clusters after each alignment
  # correction: Z_cos is refreshed from a new embedding, R/Y carry over.
  f <- load_fixture("small_block_size")
  p <- f$params
  set.seed(f$seed)
  st <- init_max_diversity_clustering(f$embedding, f$design_matrix,
                                       theta = p$theta, sigma = p$sigma, nclust = p$nclust,
                                       block_size = p$block_size, max_iter_cluster = 5, verbose = FALSE)
  st <- run_max_diversity_clustering(st)

  perturbed_embedding <- f$embedding + randn(nrow(f$embedding), ncol(f$embedding), sd = 0.01)
  st2 <- run_max_diversity_clustering(st, embedding = perturbed_embedding)
  expect_equal(st2$Z_cos, l2_normalize_columns(perturbed_embedding))
  expect_true(all(st2$R >= 0))
  expect_equal(unname(colSums(st2$R)), rep(1, ncol(st2$R)), tolerance = 1e-8)
})

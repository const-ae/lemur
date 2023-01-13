sce <- readRDS("~/Documents/PhD_Projects/NB_WaVE/NB_WaVE-Experiments/data/haber.Rds")
hvg <- order(-MatrixGenerics::rowVars(logcounts(sce)))
set.seed(1)
sce <- sce[hvg[1:500], sample.int(ncol(sce), 2000)]

fit <- differential_embedding(sce, ~ condition, n_embedding = 15, n_ambient = Inf)
DE <- test_differential_expression(fit, contrast = fact(condition = "Control") == fact(condition = "Salmonella"),
                                   return = "matrix", variance_est = "none")
fit <- add_knn_graph(fit, k = 5)

test_that("making knn graph works", {

  expect_s3_class(fit$knn_graph, "igraph")
  expect_equal(igraph::vcount(fit$knn_graph), 1000)
  expect_equal(MatrixGenerics::rowSums2(fit$knn_graph[]), rep(5L, 1000))


  fit_red <- fit[,1:10]
  expect_s3_class(fit_red$knn_graph, "igraph")
  expect_equal(igraph::vcount(fit_red$knn_graph), 10)
})


test_that("de_regions", {
  de_regions <- find_de_regions(fit, DE)
  expect_equal(colnames(de_regions), c("name", "indices", "n_cells", "mean", "sd", "z_statistic"))
  expect_equal(de_regions$z_statistic, de_regions$mean / de_regions$sd)
    expect_equal(de_regions$name, rownames(DE))
  expect_equal(sapply(seq_along(de_regions$indices), \(idx) mean(DE[idx,de_regions$indices[[idx]]])),
               de_regions$mean)
})

fit5 <- add_knn_graph(fit, k = 5)
system.time(
  de_regions <- find_de_regions(fit5, DE)
)

fit20 <- add_knn_graph(fit, k = 20)
system.time(
  de_regions <- find_de_regions(fit20, DE)
)

tmp <- find_de_regions(fit20, DE[4,,drop=FALSE])
de_regions_saved[4,]$indices
tmp$indices
expect_equal(de_regions, de_regions_saved)

find_mim_from_sel_fast(c(100, 1:10, 3), which = 2:12)

x <- rnorm(100000)
sel <- c(sample.int(length(x), size = 1000), rep(NA_integer_, 100))
find_max_from_sel_cmp <- compiler::cmpfun(find_max_from_sel)
find_max_from_sel_rcpp <- Rcpp::cppFunction(r'(

int find_max_from_sel_fast(NumericVector x, IntegerVector which){
  int val = which[0];
  double max = x[val-1];
  int which_length = which.length();
  for(int i = 1; i < which_length; i++){
    int w = which[i];
    if(IntegerVector::is_na(w)){;
     break;
    }else{
      if(x[w-1] > max){
        max = x[w-1];
        val = w;
      }
    }
  }
  return val;
}
)')
find_max_from_sel_rcpp(x, sel)
find_max_from_sel_cmp(x, sel)

bench::mark(
  sel[which.max(x[sel])],
  find_max_from_sel_fast(x, sel),
  find_max_from_sel(x, sel),
  find_max_from_sel_cmp(x, sel),
)



find_mim_from_sel_fast(c(-5, 100, -3, -5, 300),  which = c(2,3,4))
find_min_from_sel(c(-5, 100, -3, -5, 300),  which = c(2,3,4))


Rcpp::cppFunction('
// [[Rcpp::export(rng = false)]]
bool test(IntegerVector& x){
  // bool ret = false;
  // int x_size = x.size();
  // for(int i = 0; i < x_size; i++){
  //   ret = ret || IntegerVector::is_na(x[i]);
  // }
  // return ret;
  return IntegerVector::is_na(x[0]);
}')
Rcpp::sourceCpp("/tmp/test_func.cpp")

bench::mark(test(c(1:3, NA, 5)))


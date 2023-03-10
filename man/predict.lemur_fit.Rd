% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.lemur_fit}
\alias{predict.lemur_fit}
\title{Predict values from \code{lemur_fit} object}
\usage{
\method{predict}{lemur_fit}(
  object,
  newdata = NULL,
  newdesign = NULL,
  embedding = object$embedding,
  with_ambient_pca = TRUE,
  with_linear_model = TRUE,
  with_embedding = TRUE,
  with_alignment = TRUE,
  ...
)
}
\arguments{
\item{object}{an \code{lemur_fit} object}

\item{newdata}{a data.frame which passed to \code{\link{model.matrix}} with
\code{design} to make the \code{newdesign} matrix}

\item{newdesign}{a matrix with the covariates for which the output
is predicted. If \code{NULL}, the \code{object$design_matrix} is used. If
it is a vector it is repeated \code{ncol(embedding)} times to create
a design matrix with the same entry for each cell.}

\item{embedding}{the low-dimensional cell position for which the
output is predicted.}

\item{with_ambient_pca}{a boolean to indicate if the output is
the gene expression values or the value in the ambient PCA embedding.}

\item{with_linear_model}{a boolean to indicate if the linear regression
offset is included in the prediction.}

\item{with_embedding}{a boolean to indicate if the embedding contributes
to the output.}

\item{with_alignment}{a boolean to indicate if the alignment effect
is removed from the output.}

\item{...}{additional parameters passed to \code{predict_impl}.}
}
\value{
a matrix with the same dimension \code{nrow(object) * nrow(newdesign)}.
}
\description{
Predict values from \code{lemur_fit} object
}
\seealso{
\code{\link[=residuals,lemur_fit-method]{residuals}}
}

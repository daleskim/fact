% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ct_structure}
\alias{ct_structure}
\title{Correlation Thresholded Factor Analysis Structure}
\arguments{
\item{R;}{numeric matrix of correlations.}

\item{tau;}{numeric scalar of threshold \eqn{\tau} value.}

\item{k;}{optional numeric scalar integer indicating minimum clique size
to search for. Must be 2 or greater, default is 2.}

\item{check;}{optional logical scalar to indicate if \code{R} should be
checked to be a valid input.}
}
\value{
\verb{numeric matrix} of the support of the coefficients or support of
the factor loadings: \eqn{\mathcal{A}(\Lambda)}.
}
\description{
Generates a factor analysis structure through correlation
thresholding.
}
\references{
Kim, D. S., & Zhou, Q. (2023). Structure learning of latent
factors via clique search on correlation thresholded graphs. \emph{Proceedings
of the 40th International Conference on Machine Learning., 202},
16978–16996. \url{https://proceedings.mlr.press/v202/kim23aa.html}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ct_algorithm}
\alias{ct_algorithm}
\title{Factor Analysis by Correlation Thresholding}
\arguments{
\item{R;}{numeric matrix of correlations.}

\item{tau;}{optional numeric vector of threshold \eqn{\tau} values.}

\item{k;}{optional numeric scalar integer indicating minimum clique size
to search for. Must be 2 or greater, default is 2.}
}
\value{
\code{list} of \verb{numeric matrix}; each contains the support of the
coefficients or support of the factor loadings: \eqn{\mathcal{A}(\Lambda)}.
One per (non-redundant) threshold \eqn{\tau}.
}
\description{
Generates a set of factor analysis structures by thresholding
the correlation graph at a sequence of points.
}
\details{
If \code{tau} is not given, then a set is automatically generated as
follows:
\itemize{
\item If the number of observed variables is \eqn{\leq 10}, then every
possible threshold is checked. This is at most 45 thresholds.
\item If the number of observed variables is \eqn{> 10}, then 50
evenly spaced thresholds are generated from the interval \eqn{[0, 1]}.
}
Prior to clique analysis, redundant \eqn{\tau} values are removed for
computational efficiency. Therefore, if two \eqn{\tau} values result in
identical factor analysis structures, the structure will only be returned
once.
}
\references{
Kim, D. S., & Zhou, Q. (2023). Structure learning of latent
factors via clique search on correlation thresholded graphs. \emph{Proceedings
of the 40th International Conference on Machine Learning., 202},
16978–16996. \url{https://proceedings.mlr.press/v202/kim23aa.html}
}

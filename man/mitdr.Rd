% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mitdr.R
\name{mitdr}
\alias{mitdr}
\title{Integral Transformation Methods for SDR Subspaces in Multivariate Regression}
\usage{
mitdr(X,Y,d,m,method="FT-IRE",
                lambda=NA,noB = 5,noC = 20,noW = 2,sparse.cov = FALSE, x.scale = FALSE)
}
\arguments{
\item{X}{Design matrix with dimension n-by-p}

\item{Y}{Response matrix with dimension n-by-q}

\item{d}{Structure dimension (default 2).}

\item{m}{The number of omegas, i.e., 2m number of integral transforms}

\item{method}{(default ``FT-IRE'') Specify the method of dimension reduction. Other possible choices are ``FT-DIRE'',``FT-SIRE'',``FT-RIRE'', ``FT-DRIRE'', and ``admmft''.}

\item{lambda}{Tuning Parameter for ``admmft'' method. If it is not provided, the optimal lambda value is chosen by cross-validation of the Fourier transformation method.}

\item{noB}{(default 5) Iterations for updating B. Only required for the ``admmft'' method.}

\item{noC}{(default 20) Iterations for updating C. Only required for the ``admmft'' method.}

\item{noW}{(default 2) Iterations for updating weight. Only required for the ``admmft'' method.}

\item{sparse.cov}{(default FALSE) If TRUE, calculates the soft-threshold matrix. Only required for the ``admmft'' method.}

\item{x.scale}{(default FALSE) If TRUE, standardizes each variable for the soft-threshold matrix. Only required for the ``admmft'' method.}
}
\value{
The function output is a p-by-d matrix and the estimated covariance matrix.
\item{Beta_hat}{An estimator for the SDR subspace.}

\item{sigma_X}{Estimated covariance matrix only from the ``admmft'' method and a null matrix for other methods.}
}
\description{
The ``\emph{mitdr()}'' function implements transformation method for multivariate regression
}
\details{
The ``\emph{mitdr()}'' function selects the sufficient variables using Fourier transformation sparse inverse regression estimators.
}
\examples{
\dontrun{
data(prostate)
Y <- as.matrix(prostate[, 9])
X <- as.matrix(prostate[, -9])
fit.ftire <- mitdr(X, Y, d = 1, method = "FT-DRIRE")
fit.ftire$Beta_hat
}
}
\references{
Weng, J. (2022), Fourier Transform Sparse Inverse Regression Estimators for Sufficient Variable Selection,
\emph{ Computational Statistics & Data Analysis}, 168, 107380.

Weng, J., & Yin, X. (2022). A Minimum Discrepancy Approach with Fourier
Transform in Sufficient Dimension Reduction. \emph{Statistica Sinica}, 32.
}

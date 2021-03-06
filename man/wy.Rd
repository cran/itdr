% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wy.R
\name{wy}
\alias{wy}
\title{Bootstrap estimation for the tuning parameter of the response variable.}
\usage{
wy(y,x,d,wx=0.1,wy_seq=seq(0.1,1,by=0.1),wh=1.5,B=500,
                                                     xdensity="normal",method="FM")
}
\arguments{
\item{y}{The n-dimensional response vector.}

\item{x}{The design matrix of the predictors with dimension n-by-p.}

\item{d}{An integer specifying the dimension of the sufficient dimension reduction subspace.}

\item{wx}{(default 0.1). The tuning parameter for the predictor variables.}

\item{wy_seq}{(default 0.1,0.2,...,1). A sequence of the candidate tuning parameter for the response.}

\item{wh}{(default 1.5). The bandwidth of the kernel density estimation function.}

\item{B}{(default 500). Number of bootstrap samples.}

\item{xdensity}{(default ``normal''). Density function of the predictor variables. 
If ``normal'' then predictor variables are coming from a multivariate normal distribution function. 
If ``elliptic''  then predictors are coming from an elliptical contoured distribution function. 
If the distribution of the predictor variables is unknown, then use ``kernel`` to estimate the distribution 
function using kernel smoothing method.}

\item{method}{(default ``FM''). The integral transformation method. ``FM'' for Fourier trans-formation method (Zhu and Zeng 2006), and
``CM'' for convolution transformation method (Zeng and Zhu 2010).}
}
\value{
The outputs are a table of average bootstrap distances between two subspaceses for each candidate value of wy and estimate value for wy.
\item{dis_wy}{A table of average bootstrap  distances for each candidate value of wy.}

\item{wy.hat}{The estimated value for tunning parameter \eqn{\sigma_t^2}.}
}
\description{
\emph{wy()} estimates the turning parameter for the response variable which required in `FM` and `CM` methods only when estimating the central subspace.
}
\references{
Zeng P. and Zhu Y. (2010).
An Integral Transform Method for Estimating the Central Mean and Central Subspaces. \emph{Journal of Multivariate Analysis}. 101, 1, 271--290.

Zhu Y. and Zeng P. (2006).
Fourier Methods for Estimating the Central Subspace and Central Mean Subspace in Regression. \emph{Journal of the American Statistical Association}. 101, 476, 1638--1651.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsp.R
\name{dsp}
\alias{dsp}
\title{Distance Between Two Subspaces.}
\usage{
dsp(A, B)
}
\arguments{
\item{A}{A matrix with dimension p-by-d.}

\item{B}{A matrix with dimension p-by-d.}
}
\value{
Outputs are the following scale values.
\item{r}{One mines the trace correlation. That is, \eqn{r=1-\gamma} }

\item{q}{One mines the vector correlation. That is, \eqn{q=1-\theta}}
}
\description{
The ``\emph{dsp()}'' function calculates the distance between two subspaces, which are spanned by the columns of two matrices.
}
\details{
Let \bold{A} and \bold{B} be two full rank matrices of size
\eqn{p \times q}. Suppose \eqn{\mathcal{S}(\textbf{A})} and
\eqn{\mathcal{S}(\textbf{B})} are the column subspaces of matrices
\bold{A} and \bold{B}, respectively.
And, let \eqn{\lambda_i} 's with
\eqn{1 \geq \lambda_1^2 \geq \lambda_2^2 \geq,\cdots,\lambda_p^2\geq 0},
be the eigenvalues of the matrix \eqn{\textbf{B}^T\textbf{A}\textbf{A}^T\textbf{B}}.

1.Trace correlation, (Hotelling, 1936): \deqn{\gamma=\sqrt{\frac{1}{p}\sum_{i=1}^{p}\lambda_i^2}}

2.Vector correlation, (Hooper, 1959): \deqn{\theta=\sqrt{\prod_{i=1}^{p}\lambda_i^2}}
}
\references{
Hooper J. (1959). Simultaneous Equations and Canonical Correlation Theory. \emph{Econometrica} 27, 245-256.

Hotelling H. (1936). Relations Between Two Sets of Variates. \emph{Biometrika} 28, 321-377.
}

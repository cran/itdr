#' Distance between two subspaces.
#'
#' \emph{dsp()} returns the distance between two subspaces, which are spanned by the columns of two matrices.
#' @param A A matrix with dimension p-by-d.
#' @param B A matrix with dimension p-by-d.
#' @return
#' Outputs are the following scale values.
#' \item{r}{One mines the trace correlation. That is, \eqn{r=1-\gamma} }
#'
#' \item{q}{One mines the vector correlation. That is, \eqn{q=1-\theta}}
#'
#' @export
#' @details 
#' Let \bold{A} and \bold{B} be two full rank matrices of size 
#' \eqn{p \times q}. Suppose \eqn{\mathcal{S}(\textbf{A})} and 
#' \eqn{\mathcal{S}(\textbf{B})} are the column subspaces of matrices
#' \bold{A} and \bold{B}, respectively. 
#' And, let \eqn{\lambda_i} 's with 
#' \eqn{1 \geq \lambda_1^2 \geq \lambda_2^2 \geq,\cdots,\lambda_p^2\geq 0}, 
#' are the eigenvalues of the matrix \eqn{\textbf{B}^T\textbf{A}\textbf{A}^T\textbf{B}}.
#' 
#' (Trace correlation, Hotelling, 1936) \deqn{\gamma=\sqrt{\frac{1}{p}\sum_{i=1}^{p}\lambda_i^2}}
#' 
#' (Vector correlation, Hooper, 1959) \deqn{\theta=\sqrt{\prod_{i=1}^{p}\lambda_i^2}}
#' @references
#' Hooper J. (1959). Simultaneous Equations and Canonical Correlation Theory. \emph{Econometrica} 27, 245-256.
#' 
#' Hotelling H. (1936). Relations Between Two Sets of Variates. \emph{Biometrika} 28, 321-377.

dsp <- function(A,B)
{
  A.orth <- qr.Q(qr(A))
  B.orth <- qr.Q(qr(B))
  p=nrow(B.orth)
  d=ncol(B.orth)

  BAAB <- t(B.orth) %*% A.orth %*% t(A.orth) %*% B.orth
  BAAB.eig <- eigen(BAAB, only.values = T)$values

  rohat=sqrt(abs(det(t(A.orth) %*%B%*% t(B) %*% A.orth)))
  rohat=sqrt(abs(det(t(B) %*%A.orth%*% t(A.orth) %*% B)))

  h=(diag(p)-B.orth%*%t(B.orth))%*%A.orth
  msq=apply(h,2,function(x) sum(x^2))

  list(r = 1 - sqrt(mean((BAAB.eig))), q = 1 - sqrt(prod((BAAB.eig))),ro=rohat,msq=msq)

}

#' Dimension Selection Testing Methods for the Central Mean Subspace.
#'
#' @description
#' The ``\emph{d.test()}'' function provides p-values for the hypothesis tests for the dimension of the subpsace. It employs three test statistics: Cook's test, Scaled test, and Adjusted test, using Fourier transform approach for inverse dimension reduction method.
#' @usage d.test(y,x,m)
#' @param y The n-dimensional response vector.
#' @param x The design matrix of the predictors with dimension n-by-p.
#' @param m An integer specifying the dimension of the central mean reduction subspace to be tested.
#' @details
#' The null and alternative hypothesis are
#'
#' \deqn{H_0: d=m}   \deqn{H_a: d>m}
#'
#' 1. Weighted Chi-Square test statistics (Weng and Yin, 2018):
#' \deqn{\hat{\Lambda}=n\sum_{j=m+1}^{p}\hat{\lambda}_j,}
#' where \eqn{\lambda_j}'s are the eigenvalues of \eqn{\widehat{\textbf{V}}}, defined under the ``\emph{invFM()}'' function.
#'
#'
#' 2. Scaled test statistic (Bentler and Xie, 2000):
#' \deqn{\overline{T}_m=[trace(\hat{\Omega}_n)/p^{\star}]^{-1}n\sum_{j=m+1}^{p}\hat{\lambda}_j \sim \mathcal{X}^2_{p^{\star}},}
#' where \eqn{\hat{\Omega}_n} is a covariance matrix, and \eqn{p^{\star} = (p-m)(2t-m)}.
#'
#'
#' 3. Adjusted test statistic (Bentler and Xie, 2000):
#' \deqn{\tilde{T}_m=[trace(\hat{\Omega}_n)/d^{\star}]^{-1}n\sum_{j=m+1}^{p}\hat{\lambda}_j \sim \mathcal{X}^2_{d^{\star}},}
#' where \eqn{d^{\star} = [trace(\hat{\Omega}_n)]^{2}/trace(\hat{\Omega}_n^2)} .
#' @export
#' @return The \emph{d.test()} function returns a table of p-values for each test.
#'
#' @examples
#' data(pdb)
#' colnames(pdb) <- NULL
#' p <- 15
#' df <- pdb[, c(79, 73, 77, 103, 112, 115, 124, 130, 132, 145, 149, 151, 153, 155, 167, 169)]
#' dff <- as.matrix(df)
#' planingdb <- dff[complete.cases(dff), ]
#' y <- planingdb[, 1]
#' x <- planingdb[, c(2:(p + 1))]
#' x <- x + 0.5
#' xt <- cbind(
#'   x[, 1]^(.33), x[, 2]^(.33), x[, 3]^(.57), x[, 4]^(.33), x[, 5]^(.4),
#'   x[, 6]^(.5), x[, 7]^(.33), x[, 8]^(.16), x[, 9]^(.27), x[, 10]^(.5),
#'   x[, 11]^(.5), x[, 12]^(.33), x[, 13]^(.06), x[, 14]^(.15), x[, 15]^(.1)
#' )
#' m <- 1
#' W <- sapply(1, rnorm)
#' d.test(y, x, m)
#'
#' @references
#' Bentler P. M., and Xie, J. (2000). Corrections to Test Statistics in Principal Hessian Directions.
#' \emph{Statistics and Probability Letters}. 47, 381-389.
#'
#' Weng J., and Yin X. (2018). Fourier Transform Approach for Inverse Dimension Reduction Method. \emph{Journal of Nonparametric Statistics}. 30, 4, 1029-0311.

d.test <- function(y, x, m) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  w <- rnorm(m)
  pvalue <- c()
  pv <- rep(0, 3)
  k <- 0
  while (any(pv < 0.05)) {
    omega <- ome(x, y, k, w)
    stat <- n * sum(invFM(x, y, p, w)$eigenvalue[(k + 1):p])
    pv <- c(cooktest(omega, stat, 10000), scaletest(omega, stat, p, k, m), addtest(omega, stat))
    pvalue <- c(pvalue, pv)
    k <- k + 1
  }
  if (pvalue[2] == 1) {
    text <- cat(
      "\t Hypothesis Tests for selecting sufficient dimension (d)\n",
      "Null: d=m\t",
      "vs\t Alternative:", "d>m \n \n",
      "Test \t\t", "W.Ch.Sq \t\t Scaled  \t Adjusted  \n",
      "p-value \t", pvalue[1], "\t \t", pvalue[2], "\t\t", pvalue[3]
    )
  } else {
    text <- cat(
      "\t Hypothesis Tests for selecting sufficient dimension (d)\n",
      "Null: d=m\t",
      "vs\t Alternative:", "d>m \n \n",
      "Test \t\t", "W.Ch.Sq \t\t Scaled  \t Adjusted  \n",
      "p-value \t", pvalue[1], "\t", pvalue[2], "\t\t", pvalue[3]
    )
  }
  # return(text)
}

ome <- function(x, y, k, w) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  m <- length(w)
  Q <- matrix(sapply(w, function(w) c(cos(w * y), sin(w * y))), nrow = n)
  z <- stand(x)
  U <- apply(Q, 2, function(vec) apply(z * vec, 2, mean))
  svd <- svd(U, nu = p, nv = 2 * m)

  lam0 <- svd$u[, ((k + 1):p)]
  phi0 <- svd$v[, ((k + 1):(2 * m))]

  ######### Sigma
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- t(x1) %*% (x1) / n

  ######### Delta_xy
  A <- array(apply(Q, 2, function(v) {
    as.vector(apply(Q, 2, function(vec) {
      cov(x * v, x * vec)
    }))
  }), c(p, p, 4 * m^2))

  del_xy <- matrix(NA, 2 * m * p, 2 * m * p)
  for (i in 1:(2 * m)) {
    del_xy[((i - 1) * p + 1):(i * p), ] <- as.vector(A[, , ((i - 1) * (2 * m) + 1):(i * 2 * m)])
  }

  ######### Delta_y
  del_y <- apply(Q, 2, function(v) {
    as.vector(apply(Q, 2, function(vec) {
      cov(v, vec)
    }))
  })

  ######### Delta_x
  del_x <- sigma


  ######### Delta_xyy
  A <- array(apply(Q, 2, function(v) {
    as.vector(apply(Q, 2, function(vec) {
      cov(x * v, vec)
    }))
  }), c(p, 1, 4 * m^2))

  del_xyy <- matrix(NA, 2 * m * p, 2 * m * 1)
  for (i in 1:(2 * m)) {
    del_xyy[((i - 1) * p + 1):(i * p), ] <- as.vector(A[, , ((i - 1) * (2 * m) + 1):(i * 2 * m)])
  }


  ######### Delta_xyx
  A <- array(apply(Q, 2, function(v) cov(x * v, x)), c(p, p, 2 * m))

  del_xyx <- matrix(NA, 2 * m * p, p)
  for (i in 1:(2 * m)) {
    del_xyx[((i - 1) * p + 1):(i * p), ] <- as.vector(A[, , i])
  }


  ######### Delta_yx
  A <- array(apply(Q, 2, function(v) cov(v, x)), c(1, p, 2 * m))

  del_yx <- matrix(NA, 2 * m, p)
  for (i in 1:(2 * m)) {
    del_yx[i, ] <- as.vector(A[, , i])
  }


  ########## Delta #######
  del <- rbind(cbind(del_xy, del_xyy, del_xyx), cbind(t(del_xyy), del_y, del_yx), cbind(t(del_xyx), t(del_yx), del_x))

  ########## A
  mu <- apply(x, 2, mean)
  mean_Q <- apply(Q, 2, mean)
  Mu <- matrix(0, nrow = 2 * p * m, ncol = 2 * m) # middle
  L <- matrix(0, nrow = 2 * p * m, ncol = p) # last
  for (i in 1:(2 * m)) {
    Mu[((i - 1) * p + 1):(i * p), i] <- mu
    L[((i - 1) * p + 1):(i * p), ] <- mean_Q[i] * diag(1, p)
  }
  A <- cbind(diag(1, 2 * p * m), Mu, L)
  sigmainv <- matpower(sigma, -0.5)

  omega <- kronecker(t(phi0), t(lam0) %*% sigmainv) %*% A %*% del %*% t(A) %*% kronecker(phi0, sigmainv %*% lam0)
  eiv <- ome_ev <- eigen(omega)$values
  ome_tra <- sum(diag(omega))
  return(list(ome = omega, eigen = ome_ev, trace = ome_tra))
}

cooktest <- function(omega, stat, B) {
  eiv <- omega$eigen
  hist <- apply(matrix(rchisq(length(eiv) * B, 1), ncol = B) * eiv, 2, sum)
  pvalue <- mean(hist > stat)
  return(pvalue)
}

scaletest <- function(omega, stat, p, k, m) {
  p_sta <- ((p - k) * (2 * m - k))
  scale_stat <- 1 / (omega$trace / p_sta) * stat
  pchisq(scale_stat, p_sta, lower.tail = F)
}

addtest <- function(omega, stat) {
  d_sta <- (omega$trace)^2 / sum(diag(omega$ome %*% omega$ome))
  adj_stat <- 1 / (omega$trace / d_sta) * stat
  pchisq(adj_stat, d_sta, lower.tail = F)
}

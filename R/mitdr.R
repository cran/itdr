#' Integral Transformation Methods for SDR Subspaces in Multivariate Regression
#' @description
#' The ``\emph{mitdr()}'' function implements transformation method for multivariate regression
#' 
#' @details
#' The ``\emph{mitdr()}'' function selects the sufficient variables using Fourier transformation sparse inverse regression estimators.
#' @usage mitdr(X,Y,d,m,method="FT-IRE",
#'                 lambda=NA,noB = 5,noC = 20,noW = 2,sparse.cov = FALSE, x.scale = FALSE)
#'
#' @param X Design matrix with dimension n-by-p
#' @param Y Response matrix with dimension n-by-q
#' @param d Structure dimension (default 2).
#' @param m The number of omegas, i.e., 2m number of integral transforms
#' @param method (default ``FT-IRE'') Specify the method of dimension reduction. Other possible choices are ``FT-DIRE'',``FT-SIRE'',``FT-RIRE'', ``FT-DRIRE'', and ``admmft''.
#' @param lambda Tuning Parameter for ``admmft'' method. If it is not provided, the optimal lambda value is chosen by cross-validation of the Fourier transformation method.
#' @param noB (default 5) Iterations for updating B. Only required for the ``admmft'' method.
#' @param noC (default 20) Iterations for updating C. Only required for the ``admmft'' method.
#' @param noW (default 2) Iterations for updating weight. Only required for the ``admmft'' method.
#' @param sparse.cov (default FALSE) If TRUE, calculates the soft-threshold matrix. Only required for the ``admmft'' method.
#' @param x.scale (default FALSE) If TRUE, standardizes each variable for the soft-threshold matrix. Only required for the ``admmft'' method.

#'
#'
#' @return The function output is a p-by-d matrix and the estimated covariance matrix.
#' \item{Beta_hat}{An estimator for the SDR subspace.}
#'
#' \item{sigma_X}{Estimated covariance matrix only from the ``admmft'' method and a null matrix for other methods.}
#'
#' @export
#' @examples
#' \dontrun{
#' data(prostate)
#' Y <- as.matrix(prostate[, 9])
#' X <- as.matrix(prostate[, -9])
#' fit.ftire <- mitdr(X, Y, d = 1, method = "FT-DRIRE")
#' fit.ftire$Beta_hat
#'}
#' @references
#' Weng, J. (2022), Fourier Transform Sparse Inverse Regression Estimators for Sufficient Variable Selection,
#' \emph{ Computational Statistics & Data Analysis}, 168, 107380.
#'
#' Weng, J., & Yin, X. (2022). A Minimum Discrepancy Approach with Fourier
#' Transform in Sufficient Dimension Reduction. \emph{Statistica Sinica}, 32.
#'
#'
mitdr <- function(X, Y, d, m = 30, method = "FT-IRE", lambda = NA,
                  noB = 5, noC = 20, noW = 2, sparse.cov = FALSE, x.scale = FALSE) {
  n <- length(Y) ## sample size
  p <- ncol(X) ## dimension of predictors
  Beta_hat <- matrix(NA, ncol = d, nrow = p)
  covxx <- NA
  lambda_hat <- NA

  if (method == "FT-IRE") {
    ## the initial value for iteration
    W <- Create_omega(Y, m) ## generate W using Section 2.1
    KM <- sdrkernelt(X, Y, W, "ftire")
    KernelX <- KM$Kernel %*% t(KM$Kernel) ## create the kernel matrix Kechi Section 2.2
    sigma <- KM$hatSigma ## by default sample covariance matrix
    svd_FT <- geigen(KernelX, sigma) ## generalize eigenvalue decomposition
    vectors_FT <- svd_FT$vectors[, order(abs(svd_FT$values), decreasing = T)]
    intB <- vectors_FT[, (1:d), drop = F] ## the first d eigenvectors is the initial estimate

    fire <- int_ft_ire(X, Y, m, sigma, sparse = T) ## soft-thresholding by default
    ans_fire <- ite2(X, Y, d = d, intB = intB, fire$zeta, fire$Gzhalf) ## QL decomposition
    # ans_fire.gi <- ite2.gi(X, Y, d=d, intB = intB, fire$zeta, fire$Gzinvhalf) ## Generalized Inverse
    Beta_hat <- ans_fire$beta
  } else if (method == "FT-DIRE") {
    ## the initial value for iteration
    W <- Create_omega(Y, m) ## generate W using Section 2.1
    KM <- sdrkernelt(X, Y, W, "ftire")
    KernelX <- KM$Kernel %*% t(KM$Kernel) ## create the kernel matrix Kechi Section 2.2
    sigma <- KM$hatSigma ## by default sample covariance matrix
    svd_FT <- geigen(KernelX, sigma) ## generalize eigenvalue decomposition
    vectors_FT <- svd_FT$vectors[, order(abs(svd_FT$values), decreasing = T)]
    intB <- vectors_FT[, (1:d), drop = F] ## the first d eigenvectors is the initial estimate

    dire <- int_ft_dire(X, Y, m, sigma, K = 2)
    ans_dire <- ite2(X, Y, d = d, intB = intB, dire$zeta, dire$Gzhalf)
    hbeta_xire <- ans_dire$beta
  } else if (method == "FT-SIRE") {
    ## the initial value for iteration
    W <- Create_omega(Y, m) ## generate W using Section 2.1
    KM <- sdrkernelt(X, Y, W, "ftire")
    KernelX <- KM$Kernel %*% t(KM$Kernel) ## create the kernel matrix Kechi Section 2.2
    sigma <- KM$hatSigma ## by default sample covariance matrix
    svd_FT <- geigen(KernelX, sigma) ## generalize eigenvalue decomposition
    vectors_FT <- svd_FT$vectors[, order(abs(svd_FT$values), decreasing = T)]
    intB <- vectors_FT[, (1:d), drop = F] ## the first d eigenvectors is the initial estimate

    sire <- int_ft(X, Y, m, sigma)
    ans_sire <- ite2(X, Y, d = d, intB = intB, sire$zeta, sire$Gzhalf)
    Beta_hat <- ans_sire$beta
    # criteria(Gam,hbeta_sire,d)
  } else if (method == "FT-RIRE") {
    ## the initial value for iteration
    W <- Create_omega(Y, m) ## generate W using Section 2.1
    KM <- sdrkernelt(X, Y, W, "ftire")
    KernelX <- KM$Kernel %*% t(KM$Kernel) ## create the kernel matrix Kechi Section 2.2
    sigma <- KM$hatSigma ## by default sample covariance matrix
    svd_FT <- geigen(KernelX, sigma) ## generalize eigenvalue decomposition
    vectors_FT <- svd_FT$vectors[, order(abs(svd_FT$values), decreasing = T)]
    intB <- vectors_FT[, (1:d), drop = F] ## the first d eigenvectors is the initial estimate

    rfire <- int_ft_rire(X, Y, m, sigma)
    ans_rire <- ite2(X, Y, d = d, intB = intB, rfire$zeta, rfire$Gzhalf)
    Beta_hat <- ans_rire$beta
    # criteria(Gam,hbeta_rire,d)
  } else if (method == "FT-DRIRE") {
    ## the initial value for iteration
    W <- Create_omega(Y, m) ## generate W using Section 2.1
    KM <- sdrkernelt(X, Y, W, "ftire")
    KernelX <- KM$Kernel %*% t(KM$Kernel) ## create the kernel matrix Kechi Section 2.2
    sigma <- KM$hatSigma ## by default sample covariance matrix
    svd_FT <- geigen(KernelX, sigma) ## generalize eigenvalue decomposition
    vectors_FT <- svd_FT$vectors[, order(abs(svd_FT$values), decreasing = T)]
    intB <- vectors_FT[, (1:d), drop = F] ## the first d eigenvectors is the initial estimate

    rdire <- int_ft_dire(X, Y, m, sigma, K = 2, fun = int_ft_rire)
    ans_drire <- ite2(X, Y, d = d, intB = intB, rdire$zeta, rdire$Gzhalf)
    Beta_hat <- ans_drire$beta
  } else if (method == "admmft") {
    ans_admmft <- admmft(X, Y, d, m, lambda, sparse.cov = sparse.cov, scale.X = x.scale)
    Beta_hat <- ans_admmft$B
    covxx <- ans_admmft$covxx
    lambda_hat <- ans_admmft$lamb_cv
  } else {
    stop("Wrong Method! Check the spelling and try again")
  }
  return(list(Beta_hat = Beta_hat, sigma_X = covxx, lambda_hat = lambda_hat))
}

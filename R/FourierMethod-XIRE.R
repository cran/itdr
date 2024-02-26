fm_xire <- function(Y, X, d, m = 30, method = c("FT-IRE")) {
  n <- length(Y) ## sample size
  p <- ncol(X) ## dimension of predictors

  ## the initial value for iteration
  W <- Create_omega(Y, m) ## generate W using Section 2.1
  KM <- sdrkernelt(X, Y, W, "ftire")
  KernelX <- KM$Kernel %*% t(KM$Kernel) ## create the kernel matrix Kechi Section 2.2
  sigma <- KM$hatSigma ## by default sample covariance matrix
  svd_FT <- geigen(KernelX, sigma) ## generalize eigenvalue decomposition
  vectors_FT <- svd_FT$vectors[, order(abs(svd_FT$values), decreasing = T)]
  intB <- vectors_FT[, (1:d), drop = F] ## the first d eigenvectors is the initial estimate

  hbeta_xire <- matrix(NA, ncol = d, nrow = p)

  if (method == "FT-IRE") {
    fire <- int_ft_ire(X, Y, m, sigma, sparse = T) ## soft-thresholding by default
    ans_fire <- ite2(X, Y, d = d, intB = intB, fire$zeta, fire$Gzhalf) ## QL decomposition
    # ans_fire.gi <- ite2.gi(X, Y, d=d, intB = intB, fire$zeta, fire$Gzinvhalf) ## Generalized Inverse
    hbeta_xire <- ans_fire$beta
  } else if (method == "FT-DIRE") {
    dire <- int_ft_dire(X, Y, m, sigma, K = 2)
    ans_dire <- ite2(X, Y, d = d, intB = intB, dire$zeta, dire$Gzhalf)
    hbeta_xire <- ans_dire$beta
  } else if (method == "FT-SIRE") {
    sire <- int_ft(X, Y, m, sigma)
    ans_sire <- ite2(X, Y, d = d, intB = intB, sire$zeta, sire$Gzhalf)
    hbeta_xire <- ans_sire$beta
    # criteria(Gam,hbeta_sire,d)
  } else if (method == "FT-RIRE") {
    rfire <- int_ft_rire(X, Y, m, sigma)
    ans_rire <- ite2(X, Y, d = d, intB = intB, rfire$zeta, rfire$Gzhalf)
    hbeta_xire <- ans_rire$beta
    # criteria(Gam,hbeta_rire,d)
  } else if (method == "FT-DRIRE") {
    rdire <- int_ft_dire(X, Y, m, sigma, K = 2, fun = int_ft_rire)
    ans_drire <- ite2(X, Y, d = d, intB = intB, rdire$zeta, rdire$Gzhalf)
    hbeta_xire <- ans_drire$beta
  } else {
    stop("Wrong Method! Check the spelling and try again")
  }
  return(list(hbeta_xire = hbeta_xire))
}

# ##### cross-validation
# #' Cross Validation of Fourier Transforms
# #'
# #'
# #' @param X Design matrix with dimension n-by-p
# #' @param Y Response matrix with dimension n-by-q
# #' @param d Structure dimension of the SDR subspace, default value is 2
# #' @param m The number of omegas. That is, 2m integral transforms are used.
# #' @param nolamb Number of candidate tuning parameters.
# #' @param noB Iteration for update B, default value is 5
# #' @param noC Iteration for update C, default value is 20
# #' @param noW Iteration for update weight, default value is 2
# #' @param Kfold Number of cross validations, default value is 10
# #' @param sparse.cov If TRUE, then, calculate the soft-threshold matrix.
# #' @param scale.X If TRUE, standardize each variable for soft-threshold.
# #'
# #' @return The function outputs are the followings
# #' \item{B}{ An estimator for the SDR subspace}
# #'
# #' \item{lambcv}{ The optimal lambda value}
# #' @export
#'
# #' @references
# #' Weng, J. (2022),Fourier transform sparse inverse regression estimators for sufficient variable selection,
# #' Computational Statistics & Data Analysis, 168,107380.
# #'
cvft <- function(X, Y, d, m, nolamb = 30, noB = 5,
                 noC = 20, noW = 2, Kfold = 10, sparse.cov = FALSE, scale.X = FALSE) {
  y <- Y
  # Input:
  # X: n times p design matrix
  # y: n times q response
  # d: structure dimension, default = 2
  # m: the number of omega; that is, 2m integral transforms
  # nolamb: number of candidate tuning parameters
  # noB: iteration for update B, default = 5
  # noC: iteration for update C, default = 20
  # noW: iteration for update weight, default = 2
  # Kfold: number of cross validation = 10
  # sparse.cov: calcualte the soft-thresholding matrix
  # scale.X: if TRUE, standardize each variables for soft-thresholding matrix

  # Output:
  # B: estimation
  # lambcv: the optimal lamb value.
  n <- dim(X)[1]
  p <- dim(X)[2]
  lambmax <- 1
  lambmin <- lambmax / 20
  lambseq <- exp(seq(log(lambmin), log(lambmax), length.out = nolamb))
  loss <- matrix(100, nrow = Kfold, ncol = nolamb)
  ## make the splits
  flen <- floor(n / Kfold)
  fstart <- numeric(Kfold)
  fend <- fstart
  for (k in 1:(Kfold - 1)) {
    fstart[k] <- (k - 1) * flen + 1
    fend[k] <- k * flen
  }
  fstart[Kfold] <- (Kfold - 1) * flen + 1
  fend[Kfold] <- n
  ind <- sample(1:n, n, F)

  for (k in 1:Kfold) {
    print(paste("Fold -", k))
    ## start CV
    # indval = ind[1:(n/5)]
    # indtr = setdiff(1:n, indval)
    indval <- ind[fstart[k]:fend[k]]
    indtr <- setdiff(1:n, indval)

    Xtr <- X[indtr, ]
    Xval <- X[indval, ]
    Ytr <- y[indtr, , drop = F]
    Yval <- y[indval, , drop = F]

    for (i in 1:nolamb) {
      lamb <- lambseq[i]
      out <- admmft(Xtr, Ytr, d, m, lamb, noB, noC, noW, sparse.cov, scale.X)
      estB <- out$B
      sigs <- out$covxx
      if (kappa(estB) > 1e10) {
        loss[k, i] <- 100
      } else {
        if (d > 1) {
          estB <- estB %*% solve(matpower(t(estB) %*% sigs %*% estB, 0.5))
        } else {
          estB <- estB / as.vector(sqrt(t(estB) %*% sigs %*% estB))
        }
        loss[k, i] <- 1 - dcor(Xval %*% estB, Yval) / 2
      }
    }
  }

  loss_mean <- apply(loss, 2, mean)
  lambcv <- loss_mean[which.min(loss_mean)]
  out <- admmft(X, y, d, m, lambcv, noB, noC, noW, sparse.cov, scale.X)
  estB <- out$B
  covxx <- out$covxx
  return(list(B = estB, covxx = covxx, lambcv = lambcv))
}

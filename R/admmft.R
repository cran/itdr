admmft <- function(X, Y, d, m = 30, lambda = NA,
                   noB = 5, noC = 20, noW = 2, sparse.cov = FALSE, scale.X = FALSE) {
  y <- Y
  if (is.matrix(X) == FALSE) {
    stop("X must be a n by p matrix")
  }
  if (is.matrix(Y) == FALSE) {
    stop("Y must be a n by q matrix")
  }
  if (is.na(lambda) == TRUE) {
    lamb.select <- cvft(X, Y, d,
      m = 30, nolamb = 30, noB = 5,
      noC = 20, noW = 2, Kfold = 10, sparse.cov = FALSE, scale.X = FALSE
    )
    lambda <- lamb.select$lambcv
  }
  lamb <- lambda
  # Input:
  # X: n times p design matrix
  # y: n times q response
  # d: structure dimension, default = 2
  # m: the number of omega; that is, 2m integral transforms
  # lamb: tuning parameter
  # noB: iteration for update B, default = 5
  # noC: iteration for update C, default = 20
  # noW: iteration for update weight, default = 2
  # sparse.cov: calcualte the soft-thresholding matrix
  # scale.X: if TRUE, standardize each variables for soft-thresholding matrix
  # Output:
  # B: estimation
  rho <- 1
  n <- dim(X)[1]
  p <- dim(X)[2]
  weight <- rep(1, p)
  out <- kernel.FT(X, y, m)
  M <- out$M
  Ups <- out$U
  if (sparse.cov) {
    covxx <- spcovCV(X, standard = scale.X, method = "soft")$sigma
  } else {
    covxx <- cov(X)
  }

  out <- eigen(M)
  ## initial value
  B <- out$vectors[, 1:d, drop = F]

  ## update C
  updateC <- function(B, Ups) {
    Ut_B <- t(Ups) %*% B
    out <- svd(Ut_B)
    C <- out$v %*% t(out$u)
    return(C)
  }
  C <- updateC(B, Ups)
  ## update B
  epsrel <- epsabs <- 1e-4
  updateB <- function(C, weight) {
    U <- Z <- matrix(0, nrow = p, ncol = d)
    lambj <- lamb * weight
    for (i in 1:noB) {
      B <- solve(covxx + rho * diag(p), Ups %*% t(C) + rho * Z - rho * U)
      Zold <- Z
      normBU <- sqrt(diag((B + U) %*% t(B + U)))
      K <- pmax(1 - lambj / (rho * normBU), rep(0, p))
      Z <- diag(K) %*% (B + U)
      U <- U + B - Z
      epspri <- sqrt(p) * epsabs + epsrel * max(norm(B, "f"), norm(Z, "f"))
      epsdual <- sqrt(p) * epsabs + epsrel * norm(U, "f")
      if (norm(B - Z, "f") < epspri & norm(rho * (Zold - Z), "f") < epsdual) {
        break
      }
    }
    return(Z)
  }

  ## ADMM
  normB <- sqrt(diag(B %*% t(B)))
  small <- normB < 1e-10
  for (j in 1:noW) {
    Bold <- B
    err <- sum(diag(-t(Ups) %*% B %*% C + t(C) %*% t(B) %*% covxx %*% B %*% C)) + lamb * sum(weight[!small] * normB[!small])
    for (k in 1:noC) {
      errold <- err
      B <- updateB(C, weight)
      C <- updateC(B, Ups)
      normB <- sqrt(diag(B %*% t(B)))
      small <- normB < 1e-10
      err <- sum(diag(-t(Ups) %*% B %*% C + t(C) %*% t(B) %*% covxx %*% B %*% C)) + lamb * sum(weight[!small] * normB[!small])
      if (abs(err) < 1e-10) {
        break
      }
      if (abs(err - errold) / abs(errold) < 1e-4) {
        break
      }
    }
    weight[!small] <- sqrt(1 / normB[!small])
    if (sum(!small) > 1) {
      weight[!small] <- weight[!small] / min(weight[!small])
    }
    if (norm(Bold, "f") < 1e-10) {
      break
    }
    if (kappa(Bold) < 1e10 & kappa(B) < 1e10) {
      projloss <- norm(Bold %*% solve(t(Bold) %*% Bold, t(Bold)) - B %*% solve(t(B) %*% B, t(B)), "f")
    } else {
      projloss <- 100
    }
    if (projloss / norm(Bold, "f") < 1e-4) {
      break
    }
  }

  return(list(B = B, covxx = covxx, lamb_cv = lamb))
}

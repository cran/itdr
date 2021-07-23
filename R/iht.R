

iht <- function(y, x,d)
{
  x <- as.matrix(x)
  y <- scale(y)

  n <- nrow(x)
  p <- ncol(x)

  x.e <- eigen(cov(x), symmetric=T)
  halfinv <- x.e$vectors %*% diag(1/sqrt(x.e$values)) %*% t(x.e$vectors)
  z <- scale(x, center=T, scale=F) %*% halfinv

  M <- M.dimreduc.iht(y, z)

  M.raw <- eigen(M, symmetric = T)
  M.evectors <- halfinv %*% M.raw$vectors
  M.evectors <- apply(M.evectors, 2, function(x) x/sqrt(sum(x^2)))
  eta_hat=M.evectors[,c(1:d)]
  list(eta_hat=eta_hat,M=M)
}

M.dimreduc.iht <- function(y, x)
{
  p <- ncol(x)
  n <- nrow(x)

  yxmatrix <- matrix(y, ncol=p, nrow=n) * x
  Sigmayzz <- t(yxmatrix) %*% x / n;
  Betayz <- matrix(apply(yxmatrix, 2, mean), ncol=1)

  B <- matrix(0, ncol=p, nrow=p)
  B[ , 1] <- Betayz / sqrt(sum(Betayz * Betayz))

  for(i in 2:p)
  {
    B[ , i] = Sigmayzz %*% B[ , i-1]
    B[ , i] = B[ , i] / sqrt(sum(B[ , i] * B[ , i]))
  }

  B %*% t(B)
}

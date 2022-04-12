########################################################################
## Integral Transform Method for Sufficient Dimension Reduction
##
## Reference: 
##   Peng Zeng and Yu Zhu (2010). An integral transform method for 
##      estimating the central mean and central Subspaces.
##      Journal of Multivariate Analysis 101, 271--290.
##
## Peng Zeng @ Auburn University
## 03-30-2010
########################################################################
#' @useDynLib itdr, .registration = TRUE, .fixes = "C_"
#' @importFrom stats cov pchisq quantile rchisq rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics plot
#' @import MASS
FTM = function(x, y, hx, hy, h,
               space = c("mean", "pdf"),
               xdensity = c("normal", "kernel", "elliptic"),
               out = FALSE, outden = FALSE)
{
   n = NROW(x); p = NCOL(x)
   if(length(y) != n)
      stop("Error:\n\tThe length of the response is not equal to the
           number of rows of the predictors!\n")
   
   x = as.matrix(x)
   y = scale(y)
   temp = x.standardize(x)
   Sinvsqrt = temp$Sinvsqrt
   xstd = temp$x.std
   
   space = match.arg(space)
   xdensity = match.arg(xdensity)
   
   if(xdensity == "normal")
   { denidx = 0; index = NULL; f0gk = NULL }
   else
   {
      denidx = 1
      f0gk = switch(xdensity,
                    kernel    = .C("Fdlogden1", as.double(xstd), as.integer(n), 
                                   as.integer(p), as.double(h), as.integer(outden),
                                   f0 = double(n), dlogf = double(n*p)),
                    elliptic  = .C("Fdlogden3", as.double(xstd), as.integer(n), 
                                   as.integer(p), as.double(h), as.integer(outden),
                                   f0 = double(n), rdlogf = double(n), dlogf = double(n*p)))
      
      cutpoint = quantile(f0gk$f0, 0.1)  ## trim 10% points.
      index = (f0gk$f0 >= cutpoint)
      gk = matrix(f0gk$dlogf, ncol = p)[index, ]
      xstd = xstd[index, ]; y = y[index]; n = sum(index)
   } 
   
   Mfunindex = paste(space, denidx, sep="")
   M = switch(Mfunindex,
              mean0 = .C("FM_mean_norm", as.double(xstd), as.double(y),
                         as.double(hx), as.integer(n), as.integer(p),
                         as.integer(out), M = double(p * p))$M,
              mean1 = .C("FM_mean", as.double(xstd), as.double(gk),
                         as.double(y), as.double(hx), as.integer(n),
                         as.integer(p), as.integer(out), M = double(p * p))$M,
              pdf0  = .C("FM_pdf_norm", as.double(xstd), as.double(y),
                         as.double(hx), as.double(hy), as.integer(n), as.integer(p),
                         as.integer(out), M = double(p * p))$M,
              pdf1  = .C("FM_pdf", as.double(xstd), as.double(gk),
                         as.double(y), as.double(hx), as.double(hy), as.integer(n),
                         as.integer(p), as.integer(out), M = double(p * p))$M)
   
   M = matrix(M, nrow = p, ncol = p)
   M.raw = eigen(M, symmetric = T)
   M.evectors = Sinvsqrt %*% M.raw$vectors
   M.evectors = apply(M.evectors, 2, function(a) {a / sqrt(sum(a^2))})
   
   info = paste("targeted space:", space,
                "\ndensity of predictors:", xdensity, "\n")
   list(evalues = M.raw$values, evectors = M.evectors, 
        M = M, index = index, density = f0gk$f0, info = info)
}


x.standardize = function(x)
{
   x.e = eigen(cov(x), symmetric = T)	
   Sinvsqrt = x.e$vectors %*% diag(1/sqrt(x.e$values)) %*% t(x.e$vectors)
   x.std = scale(x, center = T, scale = F) %*% Sinvsqrt
   list(x.std = x.std, Sinvsqrt = Sinvsqrt)
}

dist.space <- function(A, B)
{
   A.orth <- qr.Q(qr(A))
   B.orth <- qr.Q(qr(B))
   
   BAAB <- t(B.orth) %*% A.orth %*% t(A.orth) %*% B.orth
   BAAB.eig <- eigen(BAAB, only.values = T)$values
   
   list(r = 1 - sqrt(mean(BAAB.eig)), q = 1 - sqrt(prod(BAAB.eig)))
}
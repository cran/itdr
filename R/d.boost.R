#' Bootstrap Estimation for Dimension (d) of Sufficient Dimension Reduction Subspaces.
#'
#' The function ``\emph{d.boots()}'' estimates the dimension of the central mean subspace and the central subspaces in regression.
#' @usage d.boots(y,x,wx=0.1,wy=1,wh=1.5,B=500,var_plot=FALSE,space="mean"
#'                                         ,xdensity="normal",method="FM")
#' @param y The n-dimensional response vector.
#' @param x The design matrix of the predictors with dimension n-by-p.
#' @param wx (default 0.1). The tuning parameter for predictor variables.
#' @param wy (default 1). The tuning parameter for the response variable.
#' @param wh (default 1.5). The bandwidth of the kernel density estimation.
#' @param B (default 500). Number of bootstrap samples.
#' @param var_plot (default FALSE). If TRUE, it provides the dimension variability plot.
#' @param space (default ``mean''). The defalult is ``mean'' for the central mean subspace. Other option is ``pdf'' for estimating the central subspace.
#' @param xdensity (default ``normal''). Density function of predictor variables.
#' Options are ``normal'' for multivariate normal distribution, ``elliptic''  for elliptical contoured distribution function, or ``kernel'' for estimating the distribution
#' using kernel smoothing.
#' @param method (default ``FM''). The integral transformation method. ``FM'' for Fourier trans-formation method (Zhu and Zeng 2006), and
#' ``CM'' for convolution transformation method (Zeng and Zhu 2010).
#' @return The outputs includes a table of average bootstrap distances between two subspaceses for each candidate value of \emph{d} and the estimated value for \emph{d}.
#' \item{dis_d}{A table of average bootstrap  distances for each candidate value of \emph{d}.}
#'
#' \item{d.hat}{The estimated value for \eqn{d}.}
#'
#' \item{plot}{Provides the dimension variability plot if \emph{plot=TRUE}.}
#' @export
#' @examples
#' \donttest{
#' # Use dataset available in itdr package
#' data(automobile)
#' head(automobile)
#' automobile.na <- na.omit(automobile)
#' # prepare response and predictor variables
#' auto_y <- log(automobile.na[, 26])
#' auto_xx <- automobile.na[, c(10, 11, 12, 13, 14, 17, 19, 20, 21, 22, 23, 24, 25)]
#' auto_x <- scale(auto_xx) # Standardize the predictors
#' # call to the d.boots() function with required arguments
#' d_est <- d.boots(auto_y, auto_x, var_plot = TRUE, space = "pdf", xdensity = "normal", method = "FM")
#' auto_d <- d_est$d.hat
#' }
d.boots <- function(y, x, wx = .1, wy = 1, wh = 1.5, B = 500,
                    var_plot = FALSE, space = "mean",
                    xdensity = "normal", method = "FM") {
  p <- ncol(x)
  n <- nrow(x)
  dj <- matrix(0, nrow = B, ncol = 1)
  dist.r <- matrix(0, nrow = p, ncol = 1)
  y <- as.matrix(y)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = p, style = 3)

  for (j in 1:p) {
    dd <- j
    xy.dr <- itdr(y, x, dd, wx = wx, wy = wy, wh = wh, space = space, xdensity = xdensity, method = method)
    s_d <- xy.dr$eta_hat
    dataframe <- data.frame(y, x)
    boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
    for (jj in 1:B) {
      boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
      y.boostrap <- boost.df[, 1]
      x.boostrap <- boost.df[, -c(1)]
      xy.dr <- itdr(y.boostrap, x.boostrap, dd, wx = wx, wy = wy, wh = wh, space = space, xdensity = xdensity, method = method)
      s_dj <- xy.dr$eta_hat
      dist.dj <- dsp(s_dj, s_d)
      dj[jj] <- dist.dj$r
    }
    dist.r[j, 1] <- mean(dj)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  disttab <- data.frame(d = c(1:p), dbar = dist.r)
  dseq <- c(1:p)
  maxd <- dseq[which.max(dist.r)]
  if (maxd == 1) {
    maxde1 <- dseq[which.max(dist.r[-c(maxd)])]
    dist.rds <- dist.r[1:maxde1]
  } else {
    dist.rds <- dist.r[1:maxd]
  }
  # lastd=1

  dhat <- dseq[which.min(dist.rds)]
  if (var_plot == "TRUE") {
    p1 <- plot(dist.r, type = "l", xlab = "dimension", ylab = "distance")
  } else {
    p1 <- NA
  }
  list(dis_d = disttab, d.hat = dhat, Plot = p1)
}

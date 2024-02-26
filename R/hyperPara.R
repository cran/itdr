#' Bootstrap Estimation for Hyperparameters.
#'
#' The ``\emph{hyperPara()}'' function estimates the hyperparameters that required in the Fourier transformation method.
#' @usage hyperPara(y,x,d,wx=0.1,wy=1,wh=1.5,range=seq(0.1,1,by=.1),
#' xdensity="normal", B=500,space="mean", method="FM",hyper="wy")
#' @param y The n-dimensional response vector.
#' @param x The design matrix of the predictors with dimension n-by-p.
#' @param d An integer specifying the dimension of the sufficient dimension reduction subspace.
#' @param wx (default 0.1). Tuning parameter for the predictor variables.
#' @param wy (default 1). Tuning parameter for the response variable.
#' @param wh (default 1.5). Turning parameter for the bandwidth.
#' @param range (default 0.1,0.2,...,1). A sequence of candidate values for the hyperparameter.
#' @param xdensity Density function of predictor variables.
#' @param B (default 500). Number of bootstrap samples.
#' @param space (default ``mean''). Specifies whether to estimate the central mean subspace (``mean'') or the central subspace (``pdf'').
#' @param method (default ``FM''). Integral transformation method. ``FM'' for the Fourier trans-formation method (Zhu and Zeng 2006), and
#' ``CM'' for the convolution transformation method (Zeng and Zhu 2010).
#' @param hyper (default ``wy''). The hyperparameter to be estimated. Other choices are ``wx'' and ``wy''.
#' @return The outputs are a table of average bootstrap distances between two subspaces for each candidate value of the hyper parameter.
#' \item{dis_h}{A table of average bootstrap distances for each candidate value of the hyperparameter.}
#'
#' \item{h.hat}{The estimated hyperparameter.}
#' @export
#' @references
#' Zeng P. and Zhu Y. (2010).
#' An Integral Transform Method for Estimating the Central Mean and Central Subspaces. \emph{Journal of Multivariate Analysis}. 101, 1, 271--290.
#'
#' Zhu Y. and Zeng P. (2006).
#' Fourier Methods for Estimating the Central Subspace and Central Mean Subspace in Regression. \emph{Journal of the American Statistical Association}. 101, 476, 1638--1651.

hyperPara <- function(y, x, d, wx = 0.1, wy = 1, wh = 1.5, range = seq(0.1, 1, by = .1), xdensity = "normal",
                      B = 500, space = "mean", method = "FM", hyper = "wy") {
  if (hyper == "wy") {
    wy(y, x, d, wx = wx, wy_seq = range, wh = wh, B = 500, xdensity = xdensity, method = "FM")
  } else if (hyper == "wx") {
    wx(y, x, d, wx_seq = range, wy = wx, wh = wh, B = 500, space = "mean", xdensity = xdensity, method = "FM")
  } else if (hyper == "wh") {
    wh(y, x, d, wx = wx, wy = wy, wh_seq = range, B = 500, space = "mean", method = "FM")
  } else {
    stop("A wrong hyper parameter has been selected")
  }
}

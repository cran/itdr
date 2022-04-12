#' Bootstrap estimation for the bandwidth of the Gaussian kernel density estimation.
#' 
#' \emph{wh()} estimates the bandwidth of the Gaussian kernel density estimation function if the distribution of the predictor variables is unknown.
#' @usage wh(y,x,d,wx=0.1,wy=1,wh_seq=seq(0.1,3,by=.1),B=500,space="mean",
#'                                         method="FM")
#' @param y The n-dimensional response vector.
#' @param x The design matrix of the predictors with dimension n-by-p.
#' @param d An integer specifying the dimension of the sufficient dimension reduction subspace.
#' @param wx (default 0.1). The tuning parameter for the predictor variables.
#' @param wy (default 1). The tuning parameter for the response variable.
#' @param wh_seq (default 0.1,0.2,...,3). A sequence of candidate bandwidth for the kernel smoothing method.
#' @param B (default 500). Number of bootstrap samples.
#' @param space (default ``mean''). The defalult is ``mean'' for the central mean subspace. Other option is ``pdf'' for estimating the central subspace.
#' @param method (default ``FM''). The integral transformation method. ``FM''for Fourier trans-formation method (Zhu and Zeng 2006), and
#' ``CM''for convolution transformation method (Zeng and Zhu 2010). 
#' @return The outputs are a table of average bootstrap distances between two subspaces for each candidate value of the bandwidth and the estimated value for the bandwith. 
#' \item{dis_h}{A table of average boostrap distances for each candidate value of the bandwidth.}
#' 
#' \item{h.hat}{The estimated bandwidth parameter for the Gaussian kernel function.}
#' @export
#' @details 
#' The kernel density estimation of \eqn{f_{\textbf{x}}(\textbf{x})} at a fixed point \eqn{\textbf{x}_0} is defined as
#' \deqn{\widehat{f}_{\textbf{x}_0}(\textbf{x}_0)=(nh^p)^{-1}\sum_{\ell=1}^n G\left(\frac{\textbf{x}_0-\textbf{x}_{\ell}}{h}\right),}
#' where \eqn{G(\cdot)} is a Gaussian kernel function and `h' is the bandwidth of the 
#' kernel function. We denote this parameter as `wh` in all functions.
#' @references
#' Zeng P. and Zhu Y. (2010).
#' An Integral Transform Method for Estimating the Central Mean and Central Subspaces. \emph{Journal of Multivariate Analysis}. 101, 1, 271--290.
#' 
#' Zhu Y. and Zeng P. (2006).
#' Fourier Methods for Estimating the Central Subspace and Central Mean Subspace in Regression. \emph{Journal of the American Statistical Association}. 101, 476, 1638--1651.

wh=function(y,x,d,wx=0.1,wy=1,wh_seq=seq(0.1,3,by=.1),
                B=500,space="mean",method="FM"){
  if(wh_seq[1]==0){
    stop("Error!: h Sequence should not start zero")
  }
  if(method=="FM" || method=="CM"){
  xdensity="kernel"
  h=wh_seq
  dj=matrix(0,nrow=B,ncol=1)
  dist.r=matrix(0,nrow=length(h),ncol=1)
  y=as.matrix(y)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = length(h), style = 3)

  p=ncol(x)
  n=nrow(x)
  for(j in 1:length(h)){
    H=h[j]
    xy.dr =itdr(y, x, d,wx,wy,H,space,xdensity,method)
    s_d = xy.dr$eta_hat
    dataframe=data.frame(y,x)
    boost.df=dataframe[sample(nrow(x),n,replace = TRUE),]
    for(jj in 1:B){
      boost.df=dataframe[sample(nrow(x),n,replace = TRUE),]
      y.boostrap=boost.df[,1]
      x.boostrap=boost.df[,-c(1)]
      xy.dr =itdr(y.boostrap, x.boostrap, d,wx,wy,H,space,xdensity,method)
      s_dj = xy.dr$eta_hat
      dist.dj=dsp(s_dj,s_d)
      dj[jj]=dist.dj$r
    }
    dist.r[j,1]=mean(dj)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  disttab=data.frame(h=wh_seq,dbar=dist.r)
  h.hat=wh_seq[which.min(dist.r)]
  }else{
    stop("Error!:  method should be either 'FM' or 'CM'")
  }
  list(dis_h=disttab,h.hat=h.hat)
}

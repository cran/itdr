#' Bootstrap estimation for the tuning parameter of the response variable.
#'
#' \emph{wy()} estimates the turning parameter for the response variable which required in `FM` and `CM` methods only when estimating the central subspace.
#' @usage wy(y,x,d,wx=0.1,wy_seq=seq(0.1,1,by=0.1),wh=1.5,B=500,
#'                                                      xdensity="normal",method="FM")
#' @param y The n-dimensional response vector.
#' @param x The design matrix of the predictors with dimension n-by-p.
#' @param d An integer specifying the dimension of the sufficient dimension reduction subspace.
#' @param wx (default 0.1). The tuning parameter for the predictor variables.
#' @param wy_seq (default 0.1,0.2,...,1). A sequence of the candidate tuning parameter for the response.
#' @param wh (default 1.5). The bandwidth of the kernel density estimation function.
#' @param B (default 500). Number of bootstrap samples.
#' @param xdensity (default ``normal''). Density function of the predictor variables. 
#' If ``normal'' then predictor variables are coming from a multivariate normal distribution function. 
#' If ``elliptic''  then predictors are coming from an elliptical contoured distribution function. 
#' If the distribution of the predictor variables is unknown, then use ``kernel`` to estimate the distribution 
#' function using kernel smoothing method.
#' @param method (default ``FM''). The integral transformation method. ``FM'' for Fourier trans-formation method (Zhu and Zeng 2006), and
#' ``CM'' for convolution transformation method (Zeng and Zhu 2010). 
#' @return The outputs are a table of average bootstrap distances between two subspaceses for each candidate value of wy and estimate value for wy.
#' \item{dis_wy}{A table of average bootstrap  distances for each candidate value of wy.}
#' 
#' \item{wy.hat}{The estimated value for tunning parameter \eqn{\sigma_t^2}.}
#' @export
#' @references
#' Zeng P. and Zhu Y. (2010).
#' An Integral Transform Method for Estimating the Central Mean and Central Subspaces. \emph{Journal of Multivariate Analysis}. 101, 1, 271--290.
#' 
#' Zhu Y. and Zeng P. (2006).
#' Fourier Methods for Estimating the Central Subspace and Central Mean Subspace in Regression. \emph{Journal of the American Statistical Association}. 101, 476, 1638--1651.
wy=function(y,x,d,wx=0.1,wy_seq=seq(0.1,1,by=.1),wh=1.5,
                B=500,xdensity="normal",method="FM"){
  space="pdf"
  if(wy_seq[1]==0){
    stop("Error!: h Sequence should not start zero")
  }
  if(method=="FM" || method=="CM"){
  hy=wy_seq
  dj=matrix(0,nrow=B,ncol=1)
  dist.r=matrix(0,nrow=length(hy),ncol=1)
  y=as.matrix(y)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = length(hy), style = 3)

  p=ncol(x)
  n=nrow(x)
  for(j in 1:length(hy)){
    Hy=hy[j]
    xy.dr =itdr(y, x, d,wx,Hy,wh,space="pdf",xdensity,method)
    s_d = xy.dr$eta_hat
    dataframe=data.frame(y,x)
    boost.df=dataframe[sample(nrow(x),n,replace = TRUE),]
    for(jj in 1:B){
      boost.df=dataframe[sample(nrow(x),n,replace = TRUE),]
      y.boostrap=boost.df[,1]
      x.boostrap=boost.df[,-c(1)]
      xy.dr =itdr(y.boostrap, x.boostrap, d,wx,Hy,wh,space="pdf",xdensity,method)
      s_dj = xy.dr$eta_hat
      dist.dj=dsp(s_dj,s_d)
      dj[jj]=dist.dj$r
    }
    dist.r[j,1]=mean(dj)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  disttab=data.frame(hy=wy_seq,dbar=dist.r)
  hy.hat=wy_seq[which.min(dist.r)]
  }else{
    stop("Error!:  method should be either 'FM' or 'CM'")
  }
  list(dis_wy=disttab,wy.hat=hy.hat)
}

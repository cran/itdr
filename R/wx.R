#' Bootstrap estimation for the tuning parameter for the predictor variables.
#'
#' \emph{wx()} estimates the turning parameter for the predictors which required in both `FM` and `CM` methods.
#' @usage wx(y,x,d,wx_seq=seq(0.1,5,by=.1),wy=1,wh=1.5,B=500,space="mean",
#'                                                xdensity="normal",method="FM")
#' @param y The n-dimensional response vector.
#' @param x The design matrix of the predictors with dimension n-by-p.
#' @param d An integer specifying the dimension of the sufficient dimension reduction subspace.
#' @param wx_seq (default 0.1,0.2,...,0.5). A sequence of candidiate tuning parameter for the predictor variables.
#' @param wy (default 1). The tuning parameter for the response variable.
#' @param wh (default 1.5). The bandwidth of the kernel density estimation function.
#' @param B (default 500). Number of boostrap samples.
#' @param space (default ``mean''). The defalult is ``mean'' for the central mean subspace. Other option is ``pdf'' for estimating the central subspace.
#' @param xdensity (default ``normal''). Density function of the predictor variables. 
#' If ``normal'' then predictor variables are coming from a multivariate normal distribution. 
#' If ``elliptic''  then predictors are coming from an elliptical contoured distribution. 
#' If the distribution of the predictor variables is unknown, then use ``kernel'' to estimate the distribution 
#' function using kernel smoothing method.
#' @param method (default ``FM''). The integral transformation method. ``FM'' for Fourier trans-formation method (Zhu and Zeng 2006), and
#' ``CM'' for convolution transformation method (Zeng and Zhu 2010). 
#' @return The outputs are a table of average boostrap distances between two subspaces for each candidate value of wx and estimated value of wx.
#' \item{dis_wx}{A table of average boostrap distances for each candidate value of wx.}
#' 
#' \item{wx.hat}{The estimated value for tunning parameter \eqn{\sigma_w^2}.}
#' @export
#' @references
#' Zeng P. and Zhu Y. (2010).
#' An Integral Transform Method for Estimating the Central Mean and Central Subspaces. \emph{Journal of Multivariate Analysis}. 101, 1, 271--290.
#' 
#' Zhu Y. and Zeng P. (2006).
#' Fourier Methods for Estimating the Central Subspace and Central Mean Subspace in Regression. \emph{Journal of the American Statistical Association}. 101, 476, 1638--1651.


wx=function(y,x,d,wx_seq=seq(0.1,5,by=.1),wy=1,wh=1.5,B=500,space="mean"
            ,xdensity="normal",method="FM"){
  if(wx_seq[1]==0){
    stop("Error!: hx/sw2 Sequence should not start with zero")
  }
  if(method=="FM" || method=="CM"){
  hx=wx_seq
  dj=matrix(0,nrow=B,ncol=1)
  dist.r=matrix(0,nrow=length(hx),ncol=1)
  y=as.matrix(y)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = length(hx), style = 3)
  
  p=ncol(x)
  n=nrow(x)
  for(j in 1:length(hx)){
    Hx=hx[j]
    xy.dr =itdr(y,x,d,Hx,wy,wh,space,
                         xdensity,method)
      
    s_d = xy.dr$eta_hat
    dataframe=data.frame(y,x)
    boost.df=dataframe[sample(nrow(x),n,replace = TRUE),]
    for(jj in 1:B){
      boost.df=dataframe[sample(nrow(x),n,replace = TRUE),]
      y.boostrap=boost.df[,1]
      x.boostrap=boost.df[,-c(1)]
      xy.dr =itdr(y.boostrap, x.boostrap, d,Hx,wy,wh,space,xdensity,method)
      s_dj = xy.dr$eta_hat
      dist.dj=dsp(s_dj,s_d)
      dj[jj]=dist.dj$r
    }
    dist.r[j,1]=mean(dj)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  disttab=data.frame(hx=wx_seq,dbar=dist.r)
  hx.hat=wx_seq[which.min(dist.r)]
  }else{
    stop("Error!:  method should be either 'FM' or 'CM'")
  }
  list(dis_wx=disttab,wx.hat=hx.hat)
}
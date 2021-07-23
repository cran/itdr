#'  Integral transformation Methods of Estimating Sufficient Dimension Reduction Subspaces in Regression.
#'
#' \emph{itdr()} function computes a basis for sufficient dimension reduction subspaces
#' in regression.
#'
#' @param y The n-dimensional response vector.
#' @param x The design matrix of the predictors with dimension n-by-p.
#' @param d An integer specifying the dimension of the sufficient dimension reduction subspace.
#' @param wx (default 0.1). The tuning parameter for the predictor variables.
#' @param wy (default 1). The tuning parameter for the response variable.
#' @param wh (default 1.5). The bandwidth of the kernel density estimation function.
#' @param space (default ``mean''). The defalult is ``mean'' for the central mean subspace. Other option is ``pdf'' for estimating the central subspace.
#' @param xdensity (default ``normal''). Density function of the predictor variables. 
#' If ``normal'' then predictor variables are coming from a multivariate normal distribution. 
#' If ``elliptic''  then predictors are coming from an elliptical contoured distribution. 
#' If the distribution of the predictor variables is unknown, then use ``kernel'' to estimate the distribution 
#' function using a kernel smoothing method.
#' @param method (default ``FM''). The integral transformation method. ``FM'' is for the Fourier trans-formation method (Zhu and Zeng 2006), 
#' ``CM'' for convolution transformation method (Zeng and Zhu 2010), 
#' and ``iht'' for the iterative Hessian transformation method (Cook and Li 2002).
#'
#' @return The outputs are a p-by-d matrix and a p-by-p matrix defined as follows.
#' 
#' \item{eta_hat}{The estimated p by d matrix, whose coloumns form a basis of the CMS/CS.}
#'
#' \item{M}{The estimated p by p candidate matrix.}
#'
#' @export
#' @examples 
#' library(itdr)
#' data(automobile)
#' head(automobile)
#' automobile.na=na.omit(automobile)
#' wx=.14; wy=.9;wh=1.5;d=2;p=13
#' df=cbind(automobile[,c(26,10,11,12,13,14,17,19,20,21,22,23,24,25)])
#' dff=as.matrix(df)
#' automobi=dff[complete.cases(dff),]
#' y=automobi[,1]
#' x=automobi[,c(2:14)]
#' xt=scale(x)
#' fit.F_CMS=itdr(y,xt,d,wx,wy,wh,space="pdf",xdensity = "normal",method="FM")
#' round(fit.F_CMS$eta_hat,2)
#'
#' @usage itdr(y,x,d,wx=0.1,wy=1,wh=1.5,space="mean",xdensity="normal",method="FM")
#' @details
#' Let m(\bold{x})=E[y|\bold{X}=\bold{x}]. Then, integral transformation of gradient of the mean function m(\bold{x}) is defined as
#' \deqn{\bold{\psi}(\bold{\omega}) =\int \frac{\partial}{\partial \bold{x}}m(\bold{x}) W(\bold{x},\bold{\omega})f(\bold{x})d\bold{x},}
#' where \eqn{W(\bold{x},\bold{\omega})} is said to be a nondegenerate kernel function. Set \eqn{W(\bold{x},\bold{\omega})=\exp(i\bold{\omega}^T\bold{x})} for 
#' Fourier transformation (FM) method  and \eqn{W(\bold{x},\bold{\omega})=H(\bold{x}-\bold{\omega})=(2\pi\sigma_w^2)^{-p/2}\exp(-(\bold{x}-\bold{\omega})^T(\bold{x}-\bold{\omega})/(2\sigma_w^2))} for convolution transformation (CM) method  
#' where \eqn{W(\bold{x},\bold{\omega})} is an absolutely integrable function. 
#' The candidate matrix to estimate the central mean subspace (CMS),
#' \deqn{\bold{M}_{CMS}=\int \bold{\psi}(\bold{\omega}) \bold{\psi}(\bold{\omega})^T K(\bold{\omega})d\bold{\omega}, }
#' where \eqn{K(\bold{\omega})=(2\pi \sigma_w^2)^{-p/2}\exp{(-||\bold{\omega}||}/2\sigma_w^2)} under 'FM', and \eqn{K(\bold{\omega})=1} under `CM`. 
#' Here, \eqn{\sigma_w^2} is a tuning parameter and it refers as "tuning parameter for the predictor variables" and denoted by `wx` in all functions.
#' 
#' 
#' Let \eqn{\{T_v(y)=H(y,v),~ for~~ y,v\in \mathcal{R}\}} be the family of transformations for the response variable. That is, \eqn{v \in \mathcal{R}}, the mean 
#' response of \eqn{T_v(y)} is \eqn{m(\bold{\omega},v)=E[H(y,v)| \bold{X}=\bold{x}]}. Then, integral transformation for the gradient of
#' \eqn{m(\bold{\omega},v)} is defined as
#' \deqn{\bold{\psi}(\bold{\omega},v)=\int \frac{\partial}{\partial \bold{x}}m(\bold{x},v) W(\bold{x},\bold{\omega})f(\bold{x})d\bold{x},}
#' where \eqn{W(\bold{x},\bold{\omega})} is the define as above. Then, the 
#' candidate matrix for the central subspace (CS) is defined as
#' \deqn{\bold{M}_{CS}=\int H(y_1,v)H(y_2,v)dv \int \bold{\psi}(\bold{\omega}) \bold{\psi}(\bold{\omega})^T K(\bold{\omega})d\bold{\omega},}
#' where \eqn{K(\bold{\omega})} is the same as above, and \eqn{H(y,v)=(2\pi \sigma_t^2)^{-1/2}\exp(v^2/(2\sigma_t^2))} under  `FM`, 
#' and \eqn{H(y,v)=(2\pi \sigma_t^2)^{-1/2}\exp((y-v)^2/(2\sigma_t^2))} under `CM`.
#' Here \eqn{\sigma_t^2} is a tuning parameter and it refers as the "tuning parameter for the response variable" and is denote by `wy` in all functions.
#' 
#' \bold{Remark:} There is only one tuning parameter in the candidate matrix for the estimate of the CMS, and there are two tuning
#' parameters in the candidate matrix for the estimate of the CS.
#'  
#' @references
#' Cook R. D., and Li, B., (2002). 
#' Dimension Reduction for Conditional Mean in Regression. \emph{The Annals of Statistics}. 30, 455-474.
#' 
#' Zeng P. and Zhu Y. (2010).
#' An Integral Transform Method for Estimating the Central Mean and Central Subspaces. \emph{Journal of Multivariate Analysis}. 101, 1, 271--290.
#' 
#' Zhu Y. and Zeng P. (2006).
#' Fourier Methods for Estimating the Central Subspace and Central Mean Subspace in Regression. \emph{Journal of the American Statistical Association}. 101, 476, 1638--1651.

itdr=function(y,x,d,wx=0.1,wy=1,wh=1.5,space="mean",
              xdensity="normal",method="FM")
{
  if(method=="CM"){
    xy.dr=ITM(x, y, wx, wy,wh, space, xdensity)
    eta_hat=xy.dr$evectors[ , c(1:d)]
    
  }else if(method=="FM"){
    xy.dr=FTM(x,y,wx,wy,wh,space,xdensity)
    eta_hat=xy.dr$evectors[ , c(1:d)]
    
  }else if(method=="iht"){
    xy.dr=iht(y,x,d)
    eta_hat=xy.dr$evectors[ , c(1:d)]
  }else{
    stop("Error!. Wrong method provided.")
  }
  list(eta_hat=eta_hat,M=xy.dr$M)
}
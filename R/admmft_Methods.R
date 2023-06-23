###################################
#### R project: Integral Transform
#### Paper: FT-sparse IRE
#### Jiaying Weng 
###################################

#' @import energy
##### matrix power #####
matpower <- function(a,alpha){
  # a: matrix
  # alpha: the number of power of matrix a
  # return the power matrix
  small <- .000001
  p1<-nrow(a)
  eva<-eigen(a)$values
  eve<-eigen(a)$vectors
  eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
  index<-(1:p1)[Re(eva)>small]
  evai<-eva
  evai[index]<-(eva[index])^(alpha)
  ai<-eve%*%diag(evai)%*%t(eve)
  return(ai)
}

##### create kernel matrix for FT #######
kernel.FT <- function(X,y,m=30){
  # Input:
  # X: n times p design matrix
  # y: n times q response
  # m: the number of omega; that is, 2m integral transforms
  # Return:
  # M: kernel matrix = UU^T
  # U: p times 2m matrix.
  n = dim(X)[1]
  p = dim(X)[2]
  if(is.matrix(y)){
    q = dim(y)[2]
  }else{
    y = matrix(y, ncol = 1)
    q = 1
  }
  # create omega
  sig  = 0.1 * pi^2 / median(diag(y%*%t(y)))
  W = matrix(rnorm(m*q, 0, sqrt(sig)), ncol = q)
  X0 = scale(X, scale=FALSE)
  My = matrix(0, nrow = n, ncol = 2*m)
  for(i in 1:m){
    w = W[i,]
    if(q > 1){
      w = w/(w %*% t(w))^(0.5)
    }
    yw = y %*% t(w)
    My[, 2*i - 1] = scale(cos(yw))
    My[, 2*i] = scale(sin(yw))
  }
  sigcov = cov(cbind(X0, My))
  sigxf = sigcov[1:p, (p+1):ncol(sigcov)]
  # sigff = sigcov[(p+1):ncol(sigcov), (p+1):ncol(sigcov)]
  M = sigxf %*% t(sigxf)
  return(list(M=M,U=sigxf))
}

##### sparse covariance and tune parameter #####
spcov <- function(s, lam, standard, method){
  # s: matrix
  # lam: para for tuning
  # stardard: if TRUE, standardize each varibles
  # method: if 'soft', calcualte the soft-thresholding matrix 
  # return: the sparse matrix if method is TRUE.
  if(standard){
    dhat = diag(sqrt(diag(s)))
    dhatinv = diag(1/sqrt(diag(s)))
    S = dhatinv %*% s %*% dhatinv 
    Ss = S
  }else{
    dhat = diag(dim(s)[1])
    S = s
    Ss = S
  }
  if(method == 'soft'){
    tmp = abs(S) - lam
    tmp = tmp * (tmp > 0)
    Ss = sign(S) * tmp
    Ss = Ss - diag(diag(Ss)) + diag(diag(S))
  }
  outCov = dhat %*% Ss %*% dhat
  return(outCov)
}

spcovCV <- function(X, lamseq=NULL, standard=F, method='none', nsplits=10){
  # X: predictors
  # lamseq: a sequence candidates for lambda for CV
  # stardard: if TRUE, standardize each variables
  # method: if 'soft', calcualte the soft-thresholding matrix 
  # nsplits: number of resampling
  ## return:
  # sigma: covariance matrix of x
  # bestlam: the best lambda that achieve the minimum error
  # cverr: average loss for all candidates lambda
  # ntr: the size of training set
  n = dim(X)[1]
  sampCov = cov(X) * (n-1)/n
  
  ntr = floor(n*(1-1/log(n)))
  nval = n - ntr
  
  if(is.null(lamseq) && standard){
    lamseq = seq(0, 0.95, 0.05)
  }else if(is.null(lamseq) && !standard){
    absCov = abs(sampCov)
    maxval = max(absCov - diag(diag(absCov)))
    lamseq = seq(from = 0, to = maxval, length.out = 20)
  }
  cvloss = matrix(NA, nrow = length(lamseq), ncol = nsplits)
  
  for(ns in 1:nsplits){
    ind = sample(x = 1:n, size = n, replace = F)
    indtr = ind[1:ntr]
    indval = ind[(ntr+1):n]
    sstr = cov(X[indtr,]) * (ntr-1)/ntr
    ssva = cov(X[indval,]) * (nval-1)/nval
    for(i in 1:length(lamseq)){
      outCov = spcov(sstr, lamseq[i], standard, method)
      cvloss[i,ns] = base::norm(outCov-ssva, 'F')^2
    }
  }
  
  cverr = rowMeans(cvloss)
  cvix = which.min(cverr)
  bestlam = lamseq[cvix]
  outCov1 = spcov(sampCov, bestlam, standard, method)
  
  ans <- list(sigma = outCov1,
              bestlam = bestlam,
              cverr = cverr, 
              lamseq = lamseq,
              ntr = ntr)
  return(ans)
}




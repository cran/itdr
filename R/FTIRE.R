
########################################################################################################################
# This R-file includes support functions for simulations
# please load this file before running example R files.
########################################################################################################################
#' @import geigen
#' @import stats
#' @import magic
#' @import tidyr
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

##### criteria r2 and norm #####
criteria <- function(beta,hbeta,d){
  # beta: a matrix
  # hheat: the estimated matrix
  # d: the structural dimension
  # return: R2 trace correlation
  
  # there is NA
  hbeta[is.na(hbeta)] = 0
  ## beta is zero vector.
  if(all(hbeta==0)){
    return(0)
  }
  ## standization 
  ifelse((is.vector(beta)==F&dim(beta)[2]>1), beta <- beta%*%matpower(t(beta)%*%beta,-0.5), beta <- beta/sqrt(sum(beta^2)))
  ifelse((is.vector(hbeta)==F&dim(hbeta)[2]>1), hbeta <- hbeta%*%matpower(t(hbeta)%*%hbeta,-0.5), hbeta <- hbeta/sqrt(sum(hbeta^2)))
  #standardizes hbeta
  if(is.vector(hbeta)==F&dim(hbeta)[2]>1){
    hbeta<-svd(hbeta)$u
  }
  #to calculate the different distances
  temp = svd(t(hbeta)%*%beta%*%t(beta)%*%hbeta)$d
  r2<-sqrt(mean(temp))
  if(d > 1){
    ans = base::norm(hbeta%*%t(hbeta)-beta%*%t(beta), 'F')
  }else{
    ans = base::norm(abs(hbeta) - abs(beta), '2')
  }
  
  return(c(r2,ans))
} 

##### generate w value #####
##### Section 2.2 the choice of omega on page 4
Create_omega <- function(Y,m,s = 0.1, ratio=100){
  # Y: response, $n \times q$ matrix
  # m: the number of omega
  # s: parameter for standard deviation of the normal distribution
  # ratio: cutoff for extreme values
  # return: W values, $m \times q$ matrix
  q = dim(Y)[2]
  sig1 = s * pi^2/mean(diag(Y%*%t(Y)))
  sig2 = s * pi^2/median(diag(Y%*%t(Y))) ## robust version
  (R = mean(diag(Y%*%t(Y))) / median(diag(Y%*%t(Y))))
  W1 = matrix(rnorm(m*q, 0, sqrt(sig1)), ncol = q)
  W2 = matrix(rnorm(m*q, 0, sqrt(sig2)), ncol = q)
  if(R>ratio){
    W = W2
  }else{
    W = W1
  }
  W = W/max(abs(W))
  return(W)
}

##### kernel matrix for FT #####
create_upsilon <- function(x,y,w){
  # x: predictors
  # y: response, $n \times q $
  # w: a value of W, $1 \times q $
  ## return: 
  # upsilon: the real and imaginary parts of Kechi
  # J: the cos and sin of y^Tw
  # mean_iwy: the mean values of the cos and sin of y^Tw
  n = dim(x)[1]
  p = dim(x)[2]
  if(is.matrix(y)){
    q = dim(y)[2]
    w = matrix(w, ncol = q)
  }else{
    y = matrix(y, ncol = 1)
    w = matrix(w, ncol = 1)
    q = 1
  }
  
  iwy=complex(real=cos(y %*% t(w)),imaginary=sin( y %*% t(w)))
  exp_iwy=matrix(iwy,nrow=n,ncol=p)
  mean_iwy = mean(iwy) 
  upsilon = colMeans(exp_iwy*x) - mean(iwy)* colMeans(x)
  #(apply(exp_iwy*x,2,mean) - mean(iwy)*apply(x,2,mean))
  
  return(list(upsilon = cbind(Re(upsilon),Im(upsilon)), J = cbind(Re(iwy), Im(iwy)), mean_iwy = c(Re(mean_iwy), Im(mean_iwy))))
  # Example:
  # n = 100
  # p = 500
  # tt = matrix(NA, p,p)
  # for (i in 1:p){
  #   for (j in 1:p){
  #     tt[i,j] = 0.5^(abs(i-j))
  #   }
  # }
  # X = mvrnorm(n, rep(0, p), tt)
  # B = c(rep(1,5), rep(0, p-5))
  # Y = X%*%B + sin(X%*%B) + rnorm(n, 0, 1)
  # w = rnorm(1, 0, 1)
  # ker = create_upsilon(X,Y,w)
  # str(ker)
}

##### sparse covariance #####
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

##### sparse covariance and tune parameter #####
spcovCV <- function(X, lamseq=NULL, standard=F, method='none', nsplits=10){
  # X: predictors
  # lamseq: a sequence candidates for lambda for CV
  # standard: if TRUE, calcualte the soft-thresholding matrix 
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

##### kernel matrix for FT #####
sdrkernelt <- function(X,Y,W,method,lambda0){
  # X: predictors
  # Y: response
  # W: Omega Value
  # method: SIR or FT or sparse version
  # lambda0: parameter for sparse matrix
  ## return:
  # Kernel: Kernel matrix
  # hatSigma: the estimate of covariance matrix of X
  n = dim(X)[1]
  p = dim(X)[2]
  
  if(method %in% c('sir','ssir')){
    h = max(Y)
    xbar = colMeans(X) # colMeans(X)
    Xybar = matrix(NA, nrow = p, ncol = h)
    for(k in 1:h){
      index = (Y == k)
      if(sum(index)!=0){
        Xybar[,k] = (colMeans(X[index,]) - xbar) * sqrt(sum(index)/n)##sqrt(p_i)
      }else{
        Xybar[,k] = rep(0, p)
      }
    }
  }else if (method %in% c('ftire','ftsire','sftire','sftsire')){
    Xybar = matrix(apply(W, 1, function(w) create_upsilon(X,Y,w)[[1]]),nrow=p)
  }
  hatSigma = cov(X) * (n-1)/n
  if (method %in% c('ssir','sftire','sftsire')){
    hatSigma = spcov(hatSigma, lambda0, TRUE, 'soft')
  }
  ans <- list(Kernel = Xybar,
              hatSigma = hatSigma)
  
}

##### intitialize for FT-IRE#####
#### prepare values for the optimization problem
int_ft_ire <- function(X,Y,m,sigma,sparse=T){
  # X: predictors
  # Y: response
  # m: the number of omegas
  # sigma: covariance matrix
  # sparse: If true, using the soft-thresholding covariance matrix of gamma matrix
  ## return:
  # zeta: zeta matrix in section 2
  # Gz: Gamma matrix the limiting covariance matrix of kechi
  # Gzhalf: the half power of Gz
  # Gzinvhalf: the inverse and half power of Gz
  n = nrow(X)
  p = ncol(X)
  W = Create_omega(Y,m)
  Xybar = solve(sigma) %*% matrix(apply(W, 1, function(w) create_upsilon(X,Y,w)[[1]]),nrow=p)
  J = matrix(apply(W, 1, function(w) create_upsilon(X,Y,w)[[2]]),nrow = n)
  meanJ = as.vector(apply(W, 1, function(w) create_upsilon(X,Y,w)[[3]]))
  
  vec = t(sapply(1:n, function(i) 
    as.vector( solve(sigma) %*%t (X[i,,drop=F]-apply(X,2,mean)) %*% 
                 (J[i,]-meanJ-t(X[i,]-apply(X,2,mean)) %*% Xybar ))))
  Gz = (n-1)/n*cov(vec)
  if(sparse==T){
    Gz = spcovCV(vec,standard=T,method = 'soft')$sigma
  }
  
  Gzinv = matpower(Gz, -1)
  ans <- NULL
  ans$zeta <- Xybar
  ans$Gz <- Gz
  ans$Gzhalf <- matpower(Gz, 0.5)
  ans$Gzinv <- Gzinv
  ans$Gzinvhalf <- matpower(Gzinv,0.5)
  return(ans)
}

##### iteration for FT using QL #####
#### The optimization requires \Gamma^{-1}, ite2 use QL decomposition
#### While ite2.gi employ the generalized inverse matrix.
ite2 <- function(X,Y,d,intB,zeta,Gz){
  # X: predictors
  # Y: response
  # d: structural dimension 
  # intB: the initial value of iteration
  # zeta: zeta matrix in Section 2
  # Gz: Gamma matrix
  ## return:
  # beta: the estimate B
  # nv: the estimate C
  # fn: the QDF evaluated at estimates
  # iter: the number of iteration
  n <- dim(X)[1]
  p <- dim(X)[2]
  m2 <- dim(zeta)[2]	# zeta is p by 2m
  Tt <- diag(rep(1,p))
  #B <- diag(rep(1,ncol(Tt)))[,1:d,drop=FALSE]
  B = intB
  fn <- function(B,C){
    if(d>0){
      return(n*sum( forwardsolve(t(Gz),as.vector(zeta)
                                 -as.vector(B%*%C))^2 ))
    }else{
      n*sum( forwardsolve(t(Gz),as.vector(zeta))^2 )
    }
    
  }#end of function "fn"
  
  updateC <- function(){
    C = matrix( qr.coef(qr(forwardsolve(t(Gz),kronecker(diag(rep(1,m2)),
                                                        B))),forwardsolve(t(Gz),as.vector(zeta))), nrow=d)
    C[is.na(C)] = 0
    return(C)
  }#end of function "updateC"
  
  updateB <- function() { 
    for (k in 1:d) { 
      alphak <- as.vector(zeta - 
                            B[,-k,drop=FALSE] %*% C[-k,,drop=FALSE])
      PBk <- qr(B[,-k])  
      bk <- qr.coef(
        qr(forwardsolve(t(Gz),
                        t(qr.resid(PBk,t(kronecker(C[k,],Tt)))))),
        forwardsolve(t(Gz ),as.vector(alphak)))
      bk[is.na(bk)] <- 0  # can use any OLS estimate; eg, set NA to 0
      bk <- qr.resid(PBk,bk)
      B[,k] <- bk/sqrt(sum(bk^2))
    }#end of for
    B
  }#end of function "updateB"
  
  if(!is.null(B)){
    C <- updateC()
    err <- fn(B,C)
  }
  
  iter <- 0
  
  repeat{
    if(is.null(B)){
      C <- NULL
      f <- fn(0,0)
      break
    }
    iter <- iter + 1
    B <- updateB()
    C <- updateC()
    errold <- err
    err <- f <- fn(B,C)
    if(abs(err-errold)<1e-7||iter>1000){
      #print(iter)
      break
    } 
  }#end of repeat
  
  ans <- NULL
  ans$beta <- B
  ans$nu <- C
  ans$fn <- f
  ans$iter <- iter
  return(ans)
  
}

##### iteration for FT using generalized inverse #####
ite2.gi <- function(X,Y,d,intB,zeta,Gzinvhalf){
  # X: predictors
  # Y: response
  # d: structural dimension 
  # intB: the initial value of iteration
  # zeta: zeta matrix in Section 2
  # Gzinvhalf: the inverse and half power of Gamma matrix
  ## return:
  # beta: the estimate B
  # nv: the estimate C
  # fn: the QDF evaluated at estimates
  # iter: the number of iteration
  n <- dim(X)[1]
  p <- dim(X)[2]#p <- dim(zeta)[1]
  m2 <- dim(zeta)[2]	# zeta is p by 2m
  Tt <- diag(rep(1,p))
  #B <- diag(rep(1,ncol(Tt)))[,1:d,drop=FALSE]
  B = intB
  fn <- function(B,C){
    if(d>0){
      return(n*sum( ( t(Gzinvhalf) %*% (as.vector(zeta)
                                 -as.vector(B%*%C)) )^2 ))
    }else{
      n*sum( (t(Gzinvhalf) %*% as.vector(zeta))^2 )
    }
    
  }#end of function "fn"
  
  updateC <- function(){
    C = matrix( qr.coef(qr(( t(Gzinvhalf) %*% kronecker(diag(rep(1,m2)),B) )),
                        ( t(Gzinvhalf) %*% as.vector(zeta) )), nrow=d)
    C[is.na(C)] = 0
    return(C)
  }#end of function "updateC"
  
  updateB <- function() { 
    for (k in 1:d) { 
      alphak <- as.vector(zeta - 
                            B[,-k,drop=FALSE] %*% C[-k,,drop=FALSE])
      PBk <- qr(B[,-k])  
      bk <- qr.coef(
        qr( (t(Gzinvhalf) %*%
                        t(qr.resid(PBk,t(kronecker(C[k,],Tt)))))
           ),
        (t(Gzinvhalf) %*% as.vector(alphak))
        )
      bk[is.na(bk)] <- 0  # can use any OLS estimate; eg, set NA to 0
      bk <- qr.resid(PBk,bk)
      B[,k] <- bk/sqrt(sum(bk^2))
    }#end of for
    B
  }#end of function "updateB"
  
  if(!is.null(B)){
    C <- updateC()
    err <- fn(B,C)
  }
  
  iter <- 0
  
  repeat{
    if(is.null(B)){
      C <- NULL
      f <- fn(0,0)
      break
    }
    iter <- iter + 1
    B <- updateB()
    C <- updateC()
    errold <- err
    err <- f <- fn(B,C)
    if(abs(err-errold)<1e-7||iter>1000){
      #print(iter)
      break
    } 
  }#end of repeat
  
  ans <- NULL
  ans$beta <- B
  ans$nu <- C
  ans$fn <- f
  ans$iter <- iter
  return(ans)
  
}

##### iteration for FT-DIRE #####
#### preparation for FT-DIRE
int_ft_dire <- function(X,Y,m,sigma,K,fun=int_ft_ire,...){
  # X: predictors
  # Y: response
  # m: the number of omegas
  # sigma: covariance matrix
  # K: the number of partition
  # fun: the method to use for each partition
  ## return:
  # zeta: zeta matrix in section 2
  # Gz: Gamma matrix the limiting covariance matrix of kechi
  # Gzhalf: the half power of Gz
  # Gzinvhalf: the inverse and half power of Gz
  ans = lapply(1:K, function(i) fun(X,Y,m/K,sigma,...))
  
  zeta = c()
  Gz <- Gzhalf <- Gzinv <- 0 
  for(i in 1:K){
    zeta = cbind(zeta,ans[[i]]$zeta)
    Gz = adiag(Gz, ans[[i]]$Gz)
    Gzhalf = adiag(Gzhalf, ans[[i]]$Gzhalf)
    Gzinv = adiag(Gzinv, ans[[i]]$Gzinv)
  } 
  
  
  Gz <- Gz[-1,-1]
  Gzhalf <- Gzhalf[-1,-1]
  Gzinv <- Gzinv[-1,-1]
  
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  ans$Gzhalf <- Gzhalf
  ans$Gzinv <- Gzinv
  ans$Gzinvhalf <- matpower(Gzinv, 0.5)
  return(ans)
}

##### iteration for FT-RIRE #####
#### preparation for FT-RIRE
int_ft_rire <- function(X,Y,m,sigma,sparse=T){
  # X: predictors
  # Y: response
  # m: the number of omegas
  # sigma: covariance matrix
  # sparse: If true, using the soft-thresholding covariance matrix of gamma matrix
  ## return:
  # zeta: zeta matrix in section 2
  # Gz: Gamma matrix the limiting covariance matrix of kechi
  # Gzhalf: the half power of Gz
  # Gzinvhalf: the inverse and half power of Gz
  n = nrow(X)
  p = ncol(X)
  W = Create_omega(Y,m)
  Xybar = solve(sigma) %*% matrix(apply(W, 1, function(w) create_upsilon(X,Y,w)[[1]]),nrow=p)
  J = matrix(apply(W, 1, function(w) create_upsilon(X,Y,w)[[2]]),nrow = n)
  meanJ = as.vector(apply(W, 1, function(w) create_upsilon(X,Y,w)[[3]]))
  
  vec = t(sapply(1:n, function(i) 
    as.vector( solve(sigma) %*%t (X[i,,drop=F]-apply(X,2,mean)) %*% 
                 (J[i,]-meanJ))))
  Gz = (n-1)/n*cov(vec)
  if(sparse==T){
    Gz = spcovCV(vec,standard=T,method = 'soft')$sigma
  }
  Gzinv = matpower(Gz, -1)
  ans <- NULL
  ans$zeta <- Xybar
  ans$Gz <- Gz
  ans$Gzhalf <- matpower(Gz, 0.5)
  ans$Gzinv <- Gzinv
  ans$Gzinvhalf <- matpower(Gzinv, 0.5)
  return(ans)
}

##### iteration for FT-SIRE #####
#### preparation for FT-SIRE
int_ft <- function(X,Y,m,sigma){
  # X: predictors
  # Y: response
  # m: the number of omegas
  # sigma: covariance matrix
  ## return:
  # zeta: zeta matrix in section 2
  # Gz: Gamma matrix the limiting covariance matrix of kechi
  # Gzhalf: the half power of Gz
  # Gzinvhalf: the inverse and half power of Gz
  
  p = ncol(X)
  
  W = Create_omega(Y,m)
  Xybar = solve(sigma) %*% matrix(apply(W, 1, function(w) create_upsilon(X,Y,w)[[1]]),nrow=p)
  
  
  Gzinv <- 0
  for(i in 1:(2*m)) Gzinv = adiag(Gzinv, sigma)
  Gzinv <- Gzinv[-1,-1]
  
  ans <- NULL
  ans$zeta <- Xybar
  ans$Gz <- matpower(Gzinv,-1)
  ans$Gzhalf <- matpower(Gzinv, -0.5)
  ans$Gzinv <- Gzinv
  ans$Gzinvhalf <- matpower(Gzinv, 0.5)
  return(ans)
}

##### estimate FT #####
FT <- function(X,Y,d,m,method='ftire'){
  p = ncol(X)
  W = Create_omega(Y,m)
  KM = sdrkernelt(X,Y,W,method = method)
  KernelX = KM$Kernel %*% t(KM$Kernel)
  sigma = KM$hatSigma
  svd_FT = geigen(KernelX,sigma)
  vectors_FT = svd_FT$vectors[,order(abs(svd_FT$values),decreasing = T)]
  intB = vectors_FT[,(1:d),drop=F]
  if(d == p){
    intBc = NULL
  }else{
    intBc = vectors_FT[,((d+1):p),drop=F]
  }
  ans=list(intB=intB,sigma=sigma,intBc = intBc,
           values= sort(abs(svd_FT$values),decreasing = T))
  return(ans)
}

##### Testing d using FT-IRE with QL #####
testd_fire <-function(X,Y,m,fire){
  p=ncol(X)
  pvalue <- c()
  k <- 1
  Stat <- c()
  DF <- c()
  while(T){
    intB <- FT(X,Y,k,m)$intB
    ans <- ite2(X, Y, k, intB = intB, fire$zeta, fire$Gzhalf)
    Stat <- c(Stat, ans$fn)
    DF <- c(DF, (p-k)*(2*m-k))
    pv <- pchisq(q = ans$fn, df = (p-k)*(2*m-k),lower.tail = F)
    #plot(dchisq(1:round(ans$fn),df = (p-k)*(2*m-k)))
    #print(pv)
    pvalue <- c(pvalue,pv)
    k = k+1
    if(k>p | k>2*m){ break
    }else if(pv>=0.05) break
  }
  d = length(pvalue)
  return(list(pvalue=pvalue,stat=Stat,df=DF,d=d,k=k))
}

##### Testing d using FT-IRE with Generalize Inverse#####
testd_fire.gi <-function(X,Y,m,fire){
  p=ncol(X)
  pvalue <- c()
  k <- 1
  Stat <- c()
  DF <- c()
  while(T){
    intB <- FT(X,Y,k,m)$intB
    ans <- ite2.gi(X, Y, k, intB = intB, fire$zeta, fire$Gzhalf)
    Stat <- c(Stat, ans$fn)
    DF <- c(DF, (p-k)*(2*m-k))
    pv <- pchisq(q = ans$fn, df = (p-k)*(2*m-k),lower.tail = F)
    #plot(dchisq(1:round(ans$fn),df = (p-k)*(2*m-k)))
    # print(pv)
    pvalue <- c(pvalue,pv)
    k = k+1
    if(k>p | k>2*m){ break
    }else if(pv>=0.05) break
  }
  d = length(pvalue)
  return(list(pvalue=pvalue,stat=Stat,df=DF,d=d,k=k))
}

##### construct weights for chi-square dist using FT-IRE#####
ome2 <- function(X,Y,m,dire,ans){
  p = dim(X)[2]
  Gz = dire$Gz
  ## Calcuate the delta if you have \hat B and \hat C. 
  if(is.null(ans$beta)){
    Del <- diag(1, p*m*2)
  }else{Del = cbind(kronecker(t(ans$nu), diag(1,p)), kronecker(diag(1,m*2), ans$beta))}
  
  V = dire$Gzinv
  Phi = matpower(V,0.5)%*% Del
  PPhi = Phi %*% matpower(t(Phi)%*%Phi, -1) %*% t(Phi)
  QPhi = diag(1,p*m*2) - PPhi
  Om = matpower(V, 0.5) %*% Gz %*% matpower(V,0.5)
  omega = QPhi %*% Om %*% QPhi
  ome_ev = eigen(omega)$values
  ome_tra = sum(diag(omega))
  
  return(list(ome = omega,eigen = ome_ev,trace= ome_tra))
}

##### Testing d using weighted chi-square #####
cooktest <- function(omega,stat,B){
  eiv  = omega$eigen
  if(any(is.complex(eiv))) eiv = Re(eiv)
  hist = apply(matrix(rchisq(length(eiv)*B,1),ncol=B)*eiv,2,sum)
  pvalue = mean(hist>stat)
  return(pvalue)
}

##### Testing d using scaled chi-square #####
scaletest <-function(omega,stat,p,k,m2){
  p_sta = ((p-k)*(m2-k))
  scale_stat = 1/(omega$trace/p_sta) * stat
  pchisq(scale_stat,p_sta,lower.tail = F)
}

##### Testing d using adjusted chi-square #####
addtest <-function(omega,stat){
  d_sta = (omega$trace)^2/sum(diag(omega$ome%*%t(omega$ome)))
  adj_stat = 1/(omega$trace/d_sta) * stat
  pchisq(adj_stat,d_sta,lower.tail = F)
}

##### Testing d using FT-DIRE with QL #####
testd_dire <-function(X,Y,m,dire){
  p=ncol(X)
  pvalue <- c()
  k <- 1
  Stat <- c()
  while(T){
    ft <- FT(X,Y,k,m)
    intB = ft$intB
    sigma = ft$sigma
    ans <- ite2(X, Y, k, intB = intB, dire$zeta, dire$Gz)
    Stat <- c(Stat, ans$fn)
    omega = ome2(X,Y,m,dire,ans)
    stat = ans$fn
    pv = cooktest(omega,stat,10000)
    #print(pv)
    pvalue <- c(pvalue,pv)
    k = k+1
    if(k>p | k>2*m){ break
    }else if((pv>=0.05)) break   
  }
  
  d = length(pvalue)#apply(pvalue, 2, function(x) which(x>0.05))
  
  return(list(pvalue=pvalue,stat=Stat,d=d,k=k))
}

##### Testing d using FT-DIRE with Generalize Inverse#####
testd_dire.gi <-function(X,Y,m,dire){
  p=ncol(X)
  pvalue <- c()
  k <- 1
  Stat <- c()
  while(T){
    ft <- FT(X,Y,k,m)
    intB = ft$intB
    sigma = ft$sigma
    ans <- ite2.gi(X, Y, k, intB = intB, dire$zeta, dire$Gz)
    Stat <- c(Stat, ans$fn)
    omega = ome2(X,Y,m,dire,ans)
    stat = ans$fn
    pv = cooktest(omega,stat,10000)
    # print(pv)
    pvalue <- c(pvalue,pv)
    k = k+1
    if(k>p | k>2*m){ break
    }else if((pv>=0.05)) break   
  }
  
  d = length(pvalue)#apply(pvalue, 2, function(x) which(x>0.05))
  
  return(list(pvalue=pvalue,stat=Stat,d=d,k=k))
}

##### construct weights for chi-square dist using FT Weng Yin 2018 #####
ome_ed <- function(X, Y, m, k, sigma){
 
  n = dim(X)[1]
  p = dim(X)[2]
  W = Create_omega(Y,m)
  Q =  matrix(sapply(1:m, function(i) c(cos(Y %*% W[i,]),sin(Y %*% W[i,]))),nrow = n)
  
  ##
  z = t(matpower(sigma,-0.5) %*% (t(X)-apply(X,2,mean)%*%matrix(1,nrow = 1,ncol = n)))
  U = apply(Q,2, function(vec) apply(z * vec,2,mean))
  SVD = svd(U, nu = p, nv = m*2)
  
  lam0 = SVD$u[,((k+1):p)]
  phi0 = SVD$v[,((k+1):(m*2))]
  
  ######### Delta_xy
  A= array(apply(Q,2,function(v) as.vector(apply(Q,2,function(vec) 
    cov(X*v,X*vec)))),c(p,p,(2*m)^2))
  
  del_xy = matrix(NA,(m*2)*p,m*2*p)
  for(i in 1:(m*2)){
    del_xy[((i-1)*p+1):(i*p),] = as.vector(A[,,((i-1)*(m*2)+1):(i*m*2)])
  }
  
  ######### Delta_y
  del_y = apply(Q,2,function(v) as.vector(apply(Q,2,function(vec) 
    cov(v,vec))))
  
  ######### Delta_x
  del_x = sigma
  
  
  ######### Delta_xyy
  A= array(apply(Q,2,function(v) as.vector(apply(Q,2,function(vec) 
    cov(X*v,vec)))),c(p,1,4*m^2))
  
  del_xyy = matrix(NA,m*2*p,m*2*1)
  for(i in 1:(m*2)){
    del_xyy[((i-1)*p+1):(i*p),] = as.vector(A[,,((i-1)*(m*2)+1):(i*m*2)])
  }
  
  
  ######### Delta_xyx
  A= array(apply(Q,2,function(v) cov(X*v,X)),c(p,p,m*2))
  
  del_xyx = matrix(NA,m*2*p,p)
  for(i in 1:(m*2)){
    del_xyx[((i-1)*p+1):(i*p),] = as.vector(A[,,i])
  }
  
  
  ######### Delta_yx
  A = array(apply(Q,2,function(v) cov(v,X)),c(1,p,m*2))
  
  del_yx = matrix(NA,m*2,p)
  for(i in 1:(m*2)){
    del_yx[i,] = as.vector(A[,,i])
  }
  
  
  ########## Delta #######
  del = rbind(cbind(del_xy, del_xyy, del_xyx), cbind(t(del_xyy), del_y, del_yx), 
              cbind(t(del_xyx), t(del_yx), del_x))
  
  ########## A 
  xb <- apply(X, 2, mean)
  mean_Q = apply(Q, 2, mean)
  Mu = matrix(0, nrow=p*m*2, ncol=m*2)#middle
  L = matrix(0, nrow=p*m*2, ncol=p)#last
  for(i in 1:(m*2)){
    Mu[((i-1)*p+1):(i*p),i] = xb
    L[((i-1)*p+1):(i*p),] = mean_Q[i]* diag(1,p)
  }  
  A = cbind(diag(1, p*m*2), Mu, L)
  
  #sigmainv = matpower(sigma,-0.5)
  sigmamrt = matpower(sigma,0.5)
  omega = kronecker(t(phi0),t(lam0) %*% sigmamrt) %*% A %*% del %*% t(A) %*% kronecker(phi0, sigmamrt %*% lam0)
  eiv = ome_ev = eigen(omega)$values
  ome_tra = sum(diag(omega))
  return(list(ome = omega,eigen = ome_ev,trace= ome_tra))
}

##### testing d using FT Weng Yin 2018 #####
testd <- function(X,Y,m){
  n = dim(X)[1]
  p = dim(X)[2]
  ft = FT(X,Y,p,m)
  sigma = ft$sigma
  pvalue = c()
  k=1
  Stat = c()
  while(T){
    omega = ome_ed(X,Y,m,k,sigma)
    stat = n*sum(ft$values[(k+1):p])
    pv = scaletest(omega,stat,p,k,2*m)#c(cooktest(omega,stat,10000),scaletest(omega,stat,p,k,2*m),addtest(omega,stat))
    pv
    pvalue = rbind(pvalue,pv)
    k = k+1
    if(k>p | k>2*m){ break
    }else if(all(pv>=0.05)) break   
  }
  # d = length(pvalue)
  d = apply(pvalue, 2, function(x) which(x>0.05))
  return(list(pvalue=pvalue,stat=Stat,d=d,k=k))
}

##### projection matrix #####
projection <- function(alpha){
  if(dim(alpha)[2]>1){
    P=alpha%*%matpower(t(alpha)%*%alpha, -1)%*%t(alpha)
  }else{
    P=alpha%*%t(alpha)/as.vector(t(alpha)%*%alpha)
  }
  
  if(is.vector(alpha)){
    n=length(alpha)
  }else{
    n=dim(alpha)[1]
  }
  Q=diag(1,n)-P
  return(list(P=P,Q=Q))
}

##### marginal testing predictors of FT-IRE #####
testHm <- function(X, m, H, fire){
  #intB <- ft$beta$beta.z
  #ans2 <- ite2(x, y, d, intB = intB, ftire$zeta, ftire$Gzhalf)
  n = dim(X)[1]
  test_stat <- n*matrix( t(H)%*% fire$zeta, nrow=1 ) %*% 
    matpower(kronecker(diag(2*m), t(H))%*%fire$Gz%*%kronecker(diag(2*m), H), alpha = -1) %*% 
    matrix(t(H)%*%fire$zeta, ncol=1)
  ##chi-sqare: 2rm
  r = dim(H)[2]
  df = 2*r*m
  pvalue = pchisq(test_stat,df, lower.tail = F)
  return(list(statistic=test_stat, df = df, pvalue=pvalue))
}

##### iteration  #####
ite3 <- function(X,Y,H0,t,intB,zeta,Gz){
  n <- dim(X)[1]
  p <- dim(X)[2]# p <- dim(zeta)[1]
  d=ncol(intB)
  h1 <- dim(zeta)[2]  # zeta is p by (h-1)
  
  B = intB
  Tt <- diag(rep(1,p))
  #B <- diag(rep(1,ncol(Tt)))[,1:d,drop=FALSE]
  
  fn <- function(B,C){
    if(d>0){
      return(n*sum( forwardsolve(t(Gz),as.vector(zeta)
                                 -as.vector(H0%*%B%*%C))^2 ))
    }else{
      return(n*sum( forwardsolve(t(Gz),as.vector(zeta))^2 ))
    }
    
  }#end of function "fn"
  
  updateC <- function(){
    C = matrix( qr.coef(qr(forwardsolve(t(Gz),kronecker(diag(rep(1,h1)),
                                                        H0%*%B))),
                        forwardsolve(t(Gz),as.vector(zeta))), nrow=t)
    C[is.na(C)] = 0
    return(C)
  }#end of function "updateC"
  
  updateB <- function(){
    B = matrix( qr.coef(qr(forwardsolve(t(Gz),kronecker(t(C),
                                                        H0))),
                        forwardsolve(t(Gz),as.vector(zeta))), ncol=t)
    return(B)
  }#end of function "updateB"
  
  if(!is.null(B)){
    C <- updateC()
    err <- fn(B,C)
  }
  
  iter <- 0
  repeat{
    if(is.null(B)){
      C <- NULL
      f <- fn(0,0)
      break
    }
    iter <- iter + 1
    B <- updateB()
    C <- updateC()
    errold <- err
    err <- f <- fn(B,C)
    if(abs(err-errold)/errold<1e-7||iter>200) break
  }#end of repeat
  
  ans <- NULL
  ans$beta <- B
  ans$nu <- C
  ans$fn <- f
  ans$iter <- iter
  return(ans)
  
}#end of function "ite"

##### joint testing predictors of FT-IRE #####
testHj <- function(X, Y, d, m, H, fire){
  n = dim(X)[1]
  p = dim(X)[2]
  r = dim(H)[2]
  test_stat <- n*matrix(t(H)%*% fire$zeta, nrow=1)
  Qh <- projection(alpha = H)$Q
  Qh_ei <- eigen(Qh)
  H0 <- Qh_ei$vectors[,abs(Qh_ei$values-1)<10e-10,drop=F]
  
  Stat = df = pvalue = c()
  t = 1
  while(T){
    ft = FT(X,Y,t,m)
    intB <- ft$intB
    intB_H0 = t(H0) %*% Qh %*% intB
    ans1 <- ite3(X, Y, H0, t, intB = intB_H0, zeta = fire$zeta, Gz = fire$Gzhalf)## check code for d =0
    Stat = c(Stat,ans1$fn)
    df = c(df, (p-t)*(2*m-t)+t*r)
    pv = pchisq(test_stat[t],df[t], lower.tail = F)
    pv
    pvalue = c(pvalue, pv)
    t = t+1
    if(t>p | t>2*m){ break
    }else if((pv>=0.05)) break   
    
  }
  d = length(pvalue)
  return(list(statistic=Stat, df = df, pvalue=pvalue, d=d,t= t))
}

##### conditional testing predictors of FT-IRE #####
testHd <- function(X, Y, d, m, H, fire){
  p = dim(X)[2]
  Qh <- projection(alpha = H)$Q
  Qh_ei <- eigen(Qh)
  H0 <- Qh_ei$vectors[,abs(Qh_ei$values-1)<10e-10,drop=F]
  
  ft = FT(X,Y,d,m)
  intB <- ft$intB
  intB_H0 = t(H0) %*% Qh %*% intB
  
  ans1 <- ite3(X, Y, H0, d, intB = intB_H0, zeta = fire$zeta, Gz = fire$Gzhalf)
  ans2 <- ite2(X, Y, d, intB = intB, fire$zeta, fire$Gzhalf)
  stat <- ans1$fn-ans2$fn
  df <- d*dim(H)[2]
  pvalue <- pchisq(q = stat, df = df, lower.tail = F)
  return(list(statistic=stat, df = df, pvalue=pvalue))
}

##### marginal testing predictors of FT-RIRE #####
testHm_robust <- function(X, m, H, rire){
  r = dim(H)[2]
  n = dim(X)[1]
  test_stat <- n*matrix(t(H)%*% rire$zeta, nrow=1)%*% 
    matpower(kronecker(diag(2*m), t(H))%*% rire$Gz%*%kronecker(diag(2*m), H), alpha = -1) %*% 
    matrix(t(H)%*% rire$zeta, ncol=1)
  #H0
  Qh <- projection(alpha = H)$Q
  Qh_ei <- eigen(Qh)
  H0 <- Qh_ei$vectors[,abs(Qh_ei$values-1)<10e-10,drop=F]
  #Qm
  Gnh = matpower(rire$Gzinv,alpha = 0.5) # G^{-1/2}
  Qm = projection(alpha = Gnh%*%kronecker(diag(2*m),H0))$Q
  ##
  eig = eigen(Qm%*%Gnh%*%rire$Gz%*%Gnh%*%Qm)
  eiv = eig$values
  if(any(is.complex(eiv))) eiv = Re(eiv)
  B = 10000
  hist = apply(matrix(rchisq(length(eiv)*B,1),ncol=B)*eiv,2,sum)
  pvalue = mean(hist>as.vector(test_stat))
  return(list(statistic=test_stat, pvalue=pvalue))
}

##### conditional testing predictors of FT-RIRE #####
testHd_robust <- function(X, Y, d, m, H, rire){
  p = dim(X)[2]
  Qh <- projection(alpha = H)$Q
  Qh_ei <- eigen(Qh)
  H0 <- Qh_ei$vectors[,abs(Qh_ei$values-1)<10e-10,drop=F]
  
  ft = FT(X,Y,d,m)
  intB <- ft$intB
  intB_H0 = t(H0) %*% Qh %*% intB
  ans1 <- ite3(X, Y, H0, d, intB = intB_H0, zeta = rire$zeta, Gz = rire$Gzhalf)
  ans2 <- ite2(X, Y, d, intB = intB, rire$zeta, rire$Gzhalf)
  test_stat <- ans1$fn-ans2$fn
  
  #Qj
  Gnh = matpower(rire$Gzinv,alpha = 0.5) # G^{-1/2}
  Qj = projection(alpha = Gnh%*%cbind(kronecker(t(ans1$nu),H0), 
                                      kronecker(diag(2*m), H0%*%ans1$beta)))$Q
  #Q 
  Q = projection(alpha = Gnh%*%cbind(kronecker(t(ans2$nu), diag(p)), 
                                     kronecker(diag(2*m), ans2$beta)))$Q
  #Qc
  Qc = Qj - Q
  ##
  eig = eigen(Qc%*%Gnh%*% rire$Gz%*%Gnh%*%Qc)
  eiv = eig$values
  if(any(is.complex(eiv))) eiv = Re(eiv)
  B = 10000
  hist = apply(matrix(rchisq(length(eiv)*B,1),ncol=B)*eiv,2,sum)
  pvalue = mean(hist>as.vector(test_stat))
  return(list(statistic=test_stat, pvalue=pvalue))
}

##### standardize x #####
stand <- function(x){
  n<-nrow(x)
  p<-ncol(x)
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- t(x1) %*% (x1)/n
  eva <- eigen(sigma)$values
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1/sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}

##### kernel for sir #####
kernysir <- function(yc,nslices){
  n <- length(yc)
  nwithin <- n/nslices
  ksir <- matrix(0,n,nslices)
  for(i in 1:nslices){
    tmp <- order(yc)[(nwithin*(i-1)+1):(nwithin*i)]
    ksir[tmp,i] <- 1
  }
  return(ksir)
}

##### initialize #####
ini <- function(x,y,n,p,h){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  f <- rep(1/h,h)
  J <- kernysir(y,h) #n*h
  xi <- sapply(1:h,function(j) solve(sigma)%*%(t(x)%*%J[,j,drop=F]/sum(J[,j])-apply(x,2,mean)))
  An <- qr.Q(qr(contr.helmert(h)))		
  zeta <- xi%*%diag(f)%*%An
  
  vec = t(sapply(1:n, function(i) as.vector(matpower(sigma,-0.5)%*%t(z[i,,drop=F])%*%
                                              (J[i,]-f-(x[i,,drop=F]-apply(x,2,mean))%*%sweep(xi,2,FUN="*",f)))))
  Gz <- chol(kronecker(t(An),diag(rep(1,p))) %*% (((n-1)/n)*cov(vec)) %*%
               kronecker(An,diag(rep(1,p))))
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}

##### initialize 2 #####
ini2 <- function(x,y,n,p,H,w){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  nh <- length(H)
  
  slice_info <- function(x,y,n,p,h,w){
    f <- rep(1/h,h)
    J <- kernysir(y,h) #n*h
    xi <- sapply(1:h,function(j) solve(sigma)%*%(t(x)%*%J[,j,drop=F]/sum(J[,j])-apply(x,2,mean)))
    An <- qr.Q(qr(contr.helmert(h)))		
    zeta <- xi%*%diag(f)%*%An
    
    vec = t(sapply(1:n, function(i) 
      as.vector(matpower(sigma,-0.5)%*%t(z[i,,drop=F])%*%
                  (J[i,]-f-(x[i,,drop=F]-apply(x,2,mean))%*%sweep(xi,2,FUN="*",f)))))
    
    Gz.ih <- matpower(kronecker(t(An),diag(rep(1,p))) %*% (((n-1)/n)*cov(vec)) %*%
                        kronecker(An,diag(rep(1,p))),0.5)
    Gz.ih <- Gz.ih*w
    return(list(zeta = zeta, Gz.ih = Gz.ih))
  }
  
  ans = sapply(1:length(H),function(i)
    slice_info(x,y,n,p,H[i],w[i]))
  zeta = c()
  Gz = 0
  for(i in (1:length(H))*2-1) zeta = cbind(zeta,ans[[i]])
  for(i in (1:length(H))*2) Gz = adiag(Gz,ans[[i]])
  
  Gz <- Gz[-1,-1]
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}

##### initialize 3 #####
ini3 <- function(x,y,n,p,H,w){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  nh <- length(H)
  zeta <- NULL
  vect <- NULL
  AA <- 0
  
  slice_info <- function(x,y,n,p,h,w){
    f <- rep(1/h,h)
    J <- kernysir(y,h) #n*h
    xi <- sapply(1:h,function(j) solve(sigma)%*%(t(x)%*%J[,j,drop=F]/sum(J[,j])-apply(x,2,mean)))
    An <- qr.Q(qr(contr.helmert(h)))		
    zeta <- xi%*%diag(f)%*%An
    
    vec = t(sapply(1:n, function(i) 
      as.vector(matpower(sigma,-0.5)%*%t(z[i,,drop=F])%*%
                  (J[i,]-f-(x[i,,drop=F]-apply(x,2,mean))%*%sweep(xi,2,FUN="*",f)))))
    
    A <- kronecker(t(An),diag(rep(1,p)))
    
    return(list(zeta = zeta, vec = vec, A = A))
  }
  ans = sapply(1:length(H),function(i)
    slice_info(x,y,n,p,H[i],w[i]))
  
  for(i in 1:length(H)){
    zeta = cbind(zeta,ans[[i*3-2]])
    vect = cbind(vect,ans[[i*3-1]])
    AA <- adiag(AA,ans[[i*3]])
  }
  AA <- AA[-1,-1]
  #	Gz <- chol(AA %*% (((n-1)/n)*cov(vect)) %*% t(AA))
  Gz <- matpower(AA %*% (((n-1)/n)*cov(vect)) %*% t(AA),0.5)
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}

##### initialize for sir #####
ini_sir <- function(x,y,n,p,h){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  f <- rep(1/h,h)
  J <- kernysir(y,h) #n*h
  xi <- sapply(1:h,function(j) solve(sigma)%*%(t(x)%*%J[,j,drop=F]/sum(J[,j])-apply(x,2,mean)))
  #An <- qr.Q(qr(contr.helmert(h)))    
  zeta <- xi%*%diag(f)
  
  #     vec = t(sapply(1:n, function(i) as.vector(matpower(sigma,-0.5)%*%t(z[i,,drop=F])%*%
  #                                                (J[i,]-f-(x[i,,drop=F]-apply(x,2,mean))%*%sweep(xi,2,FUN="*",f)))))
  #     Gz <- chol(kronecker(t(An),diag(rep(1,p))) %*% (((n-1)/n)*cov(vec)) %*%
  #                  kronecker(An,diag(rep(1,p))))
  Gz <- 0
  for(ih in 1:h){
    Gz <- adiag(Gz, sigma/f[ih])
  }
  Gz <- Gz[-1,-1]    
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}

#SIR initial Ni Cook 2005
ini_sir_ni <- function(x,y,n,p,h){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  f <- rep(1/h,h)
  J <- kernysir(y,h) #n*h
  xi <- sapply(1:h,function(j) solve(sigma)%*%(t(x)%*%J[,j,drop=F]/sum(J[,j])-apply(x,2,mean)))
  #An <- qr.Q(qr(contr.helmert(h)))    
  zeta <- xi%*%diag(sqrt(f))
  
  #     vec = t(sapply(1:n, function(i) as.vector(matpower(sigma,-0.5)%*%t(z[i,,drop=F])%*%
  #                                                (J[i,]-f-(x[i,,drop=F]-apply(x,2,mean))%*%sweep(xi,2,FUN="*",f)))))
  #     Gz <- chol(kronecker(t(An),diag(rep(1,p))) %*% (((n-1)/n)*cov(vec)) %*%
  #                  kronecker(An,diag(rep(1,p))))
  Gz <- 0
  for(ih in 1:h){
    Gz <- adiag(Gz, sigma)
  }
  Gz <- Gz[-1,-1]    
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}

##### initialize #####
ini_r <- function(x,y,n,p,h){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  f <- rep(1/h,h)
  J <- kernysir(y,h) #n*h
  xi <- matrix(0,nrow=p,ncol=h)
  for(j in 1:h){
    xi[,j] <- solve(sigma)%*%(t(x)%*%J[,j,drop=F]/sum(J[,j])-apply(x,2,mean))
  }#end of for
  An <- qr.Q(qr(contr.helmert(h)))		
  zeta <- xi%*%diag(f)%*%An
  vec <- NULL			
  for(i in 1:n){
    epsilon <- J[i,]-f
    v <- as.vector(t(z[i,,drop=F])%*%epsilon)
    vec <- rbind(vec,v)
  }#end of for
  Gz <- matpower(kronecker(t(An),matpower(sigma,-0.5)) %*% (((n-1)/n)*cov(vec)) %*%
                   kronecker(An,matpower(sigma,-0.5)),0.5)
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}#end of function "ini"

##### initialize 2 #####
ini2_r <- function(x,y,n,p,H,w){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  nh <- length(H)
  zeta <- NULL
  Gz <- 0
  for(ih in 1:nh){
    h <- H[ih]
    f <- rep(1/h,h)
    J <- kernysir(y,h) #n*h
    xi <- matrix(0,nrow=p,ncol=h)
    for(j in 1:h){
      xi[,j] <- solve(sigma)%*%(t(x)%*%J[,j,drop=F]/
                                  sum(J[,j])-apply(x,2,mean))
    }#end of for
    An <- qr.Q(qr(contr.helmert(h)))		
    zeta <- cbind(zeta,xi%*%diag(f)%*%An)
    
    vec <- NULL			
    for(i in 1:n){
      epsilon <- J[i,]-f
      v <- as.vector(t(z[i,,drop=F])%*%epsilon)
      vec <- rbind(vec,v)
    }#end of for
    #Gz.ih <- chol(kronecker(t(An),matpower(sigma,-0.5)) %*% (((n-1)/n)*cov(vec)) %*%
    #                kronecker(An,matpower(sigma,-0.5)))
    Gz.ih <- matpower(kronecker(t(An),matpower(sigma,-0.5)) %*% (((n-1)/n)*cov(vec)) %*%
                        kronecker(An,matpower(sigma,-0.5)),0.5)
    Gz.ih <- Gz.ih
    Gz <- adiag(Gz,Gz.ih)
  }#end of for
  Gz <- Gz[-1,-1]
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}#end of function "ini2"

##### initialize 3 #####
ini3_r <- function(x,y,n,p,H,w){
  z <- stand(x)
  sigma <- t(x)%*%x/n
  nh <- length(H)
  zeta <- NULL
  vect <- NULL
  AA <- 0
  for(ih in 1:nh){
    h <- H[ih]
    f <- rep(1/h,h)
    J <- kernysir(y,h) #n*h
    xi <- matrix(0,nrow=p,ncol=h)
    for(j in 1:h){
      xi[,j] <- solve(sigma)%*%(t(x)%*%J[,j,drop=F]/
                                  sum(J[,j])-apply(x,2,mean))
    }#end of for
    An <- qr.Q(qr(contr.helmert(h)))		
    zeta <- cbind(zeta,xi%*%diag(f)%*%An)
    
    vec <- NULL			
    for(i in 1:n){
      epsilon <- J[i,]-f
      v <- as.vector(t(z[i,,drop=F])%*%epsilon)
      vec <- rbind(vec,v)
    }#end of for
    
    vect <- cbind(vect, vec)
    A <- kronecker(t(An),matpower(sigma,-0.5))
    AA <- adiag(AA,A)
  }#end of for
  AA <- AA[-1,-1]
  #	Gz <- chol(AA %*% (((n-1)/n)*cov(vect)) %*% t(AA))
  Gz <- matpower(AA %*% (((n-1)/n)*cov(vect)) %*% t(AA),0.5)
  ans <- NULL
  ans$zeta <- zeta
  ans$Gz <- Gz
  return(ans)
}#end of function "ini3"

##### iteration #####
ite <- function(n,d,zeta,Gz){
  p <- dim(zeta)[1]
  h1 <- dim(zeta)[2]	# zeta is p by (h-1)
  T <- diag(rep(1,p))
  if(d>=1){
    B <- diag(rep(1,ncol(T)))[,1:d,drop=FALSE]
  }else{B <- NULL}
  
  fn <- function(B,C){
    n*sum( forwardsolve(t(Gz),as.vector(zeta)
                        -as.vector(T%*%B%*%C))^2 )
  }#end of function "fn"
  
  updateC <- function(){
    matrix( qr.coef(qr(forwardsolve(t(Gz),kronecker(diag(rep(1,h1)),
                                                    T%*%B))),forwardsolve(t(Gz),as.vector(zeta))), nrow=d)
  }#end of function "updateC"
  
  updateB <- function() { 
    for (k in 1:d) { 
      alphak <- as.vector(zeta - 
                            T %*% B[,-k,drop=FALSE] %*% C[-k,])
      PBk <- qr(B[,-k])  
      bk <- qr.coef(
        qr(forwardsolve(t(Gz),
                        t(qr.resid(PBk,t(kronecker(C[k,],T)))))),
        forwardsolve(t(Gz),as.vector(alphak)))
      bk[is.na(bk)] <- 0  # can use any OLS estimate; eg, set NA to 0
      bk <- qr.resid(PBk,bk)
      B[,k] <- bk/sqrt(sum(bk^2))
    }#end of for
    B
  }#end of function "updateB"
  
  if(!is.null(B)){
    C <- updateC()
    err <- fn(B,C)
  }else{
    C <- 0
  }
  
  iter <- 0
  repeat{
    if(is.null(B)){
      B = matrix(0,nrow=p,ncol=1)
      C = 0
      break
    }
    iter <- iter + 1
    B <- updateB()
    C <- updateC()
    errold <- err
    err <- fn(B,C)
    if(abs(err-errold)/errold<1e-7||iter>200) break
  }#end of repeat
  
  ans <- NULL
  ans$beta <- B
  ans$nu <- C
  ans$fn <- fn(B,C)
  ans$iter <- iter
  return(ans)
  
}#end of function "ite"



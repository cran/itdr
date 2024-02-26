invFM <- function(x, y, d, w, x_scale = TRUE) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  phiW <- matrix(sapply(w, create_phi, x = x, y = y), nrow = p)
  ei <- eigen(phiW %*% t(phiW))
  # beta=ei$vectors[,1:d]
  if (x_scale) {
    beta <- matpower(cov(x) * (n - 1) / n, -0.5) %*% ei$vectors[, 1:d]
  } else {
    beta <- ei$vectors[, 1:d]
  }
  # beta_x = matpower(cov(x)*(n-1)/n,-0.5)%*%beta
  return(list(beta = beta, eigenvalue = ei$values, psi = phiW %*% t(phiW)))
}
create_phi <- function(x, y, w) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  q <- nrow(y)
  if (is.null(q) == FALSE) {
    iwy <- complex(real = cos(t(w) %*% y), imaginary = sin(t(w) %*% y))
  } else {
    iwy <- complex(real = cos(w * y), imaginary = sin(w * y))
  }
  im <- matrix(iwy, nrow = n, ncol = p)
  z <- stand(x)
  phi <- apply(im * z, 2, mean)
  return(cbind(Re(phi), Im(phi)))
}



criteria <- function(beta, hbeta, d) {
  # standardizes hbeta
  hbeta <- svd(hbeta)$u
  # to calculate the different distances
  r2 <- sqrt(mean(svd(t(hbeta) %*% beta %*% t(beta) %*% hbeta)$d))
  r1 <- sqrt(prod(svd(t(hbeta) %*% beta %*% t(beta) %*% hbeta)$d))
  deltam <- max(svd(hbeta %*% t(hbeta) - beta %*% t(beta))$d)
  deltaf <- sqrt(sum(diag((hbeta %*% t(hbeta) - beta %*% t(beta))
  %*% t(hbeta %*% t(hbeta) - beta %*% t(beta)))))

  m2 <- c()
  m <- c()
  for (i in 1:d) {
    m2 <- c(m2, t(hbeta[, i] - beta %*% t(beta) %*% hbeta[, i]) %*% (hbeta[, i] - beta %*% t(beta) %*% hbeta[, i]))
    m <- c(m, sqrt(m2[i]))
  }
  if (d == 1) beta <- beta / as.numeric(sqrt(t(beta) %*% beta))
  if (d > 1) beta <- beta %*% matpower(t(beta) %*% beta, -0.5)
  Proj <- beta %*% solve(t(beta) %*% beta) %*% t(beta)
  if (d == 1) ang <- acos(t(hbeta) %*% Proj %*% hbeta) * 180 / pi
  if (d > 1) ang <- apply(hbeta[, 1:d], 2, function(x) acos(t(x) %*% Proj %*% x) * 180 / pi)
  return(c(r2, r1, deltam, deltaf, m2, m, ang))
}
#   m12<-t(hbeta[,1]-beta%*%t(beta)%*%hbeta[,1])%*%(hbeta[,1]-beta%*%t(beta)%*%hbeta[,1])
#   m1<-sqrt(m12)
#
#
#   if(d>1){
#     m22<-t(hbeta[,2]-beta%*%t(beta)%*%hbeta[,2])%*%(hbeta[,2]-beta%*%t(beta)%*%hbeta[,2])
#     m2<-sqrt(m22)
#   }
#   if(d>2){
#     m32<-t(hbeta[,3]-beta%*%t(beta)%*%hbeta[,3])%*%(hbeta[,3]-beta%*%t(beta)%*%hbeta[,3])
#     m3<-sqrt(m22)
#   }
#   if(d>3){
#     m42<-t(hbeta[,4]-beta%*%t(beta)%*%hbeta[,4])%*%(hbeta[,4]-beta%*%t(beta)%*%hbeta[,4])
#     m4<-sqrt(m42)
#   }
#
#   if(d==1)  return(c(r2,r1,deltam,deltaf,m12,m1,ang))
#   else if(d==2) return(c(r2,r1,deltam,deltaf,m12,m1,m22,m2,ang))
#   else if(d==3) return(c(r2,r1,deltam,deltaf,m12,m1,m22,m2,m32,m3,ang))
#   else if(d==4) return(c(r2,r1,deltam,deltaf,m12,m1,m22,m2,m32,m3,m42,m4,ang))


# ome <- function(x,y,k,w){
#   n = dim(x)[1]
#   p = dim(x)[2]
#   m = length(w)
#   Q =  matrix(sapply(w, function(w) c(cos(w*y),sin(w*y))),nrow = n)
#   z = stand(x)
#   U = apply(Q,2, function(vec) apply(z * vec,2,mean))
#   svd = svd(U, nu = p, nv = 2*m)
#
#   lam0 = svd$u[,((k+1):p)]
#   phi0 = svd$v[,((k+1):(2*m))]
#
#   A = array(apply(Q,2,function(v) as.vector(apply(Q,2,function(vec)
#     t(lam0)%*%cov(z*v,z*vec)%*%lam0))),c((p-k),(p-k),4*m^2))
#
#   del= matrix(NA,2*m*(p-k),2*m*(p-k))
#   for(i in 1:(2*m)){
#     del[((i-1)*(p-k)+1):(i*(p-k)),] = as.vector(A[,,((i-1)*(2*m)+1):(i*2*m)])
#   }
#   omega = kronecker(t(phi0),diag(1,p-k))%*%del%*%kronecker(phi0,diag(1,p-k))
#   eiv = ome_ev = eigen(omega)$values
#   ome_tra = sum(diag(omega))
#   return(list(ome = omega,eigen = ome_ev,trace= ome_tra))
# }


# ome <- function(x,y,k,w){
#   n = dim(x)[1]
#   p = dim(x)[2]
#   m = length(w)
#   Q =  matrix(sapply(w, function(w) c(cos(w*y),sin(w*y))),nrow = n)
#   z = stand(x)
#   U = apply(Q,2, function(vec) apply(z * vec,2,mean))
#   svd = svd(U, nu = p, nv = 2*m)
#
#   lam0 = svd$u[,((k+1):p)]
#   phi0 = svd$v[,((k+1):(2*m))]
#
#   ######### Sigma
#   xb <- apply(x, 2, mean)
#   xb <- t(matrix(xb, p, n))
#   x1 <- x - xb
#   sigma <- t(x1) %*% (x1)/n
#
#   ######### Delta_xy
#   A= array(apply(Q,2,function(v) as.vector(apply(Q,2,function(vec)
#     cov(x*v,x*vec)))),c(p,p,4*m^2))
#
#   del_xy = matrix(NA,2*m*p,2*m*p)
#   for(i in 1:(2*m)){
#     del_xy[((i-1)*p+1):(i*p),] = as.vector(A[,,((i-1)*(2*m)+1):(i*2*m)])
#   }
#
#   ######### Delta_y
#   del_y = apply(Q,2,function(v) as.vector(apply(Q,2,function(vec)
#     cov(v,vec))))
#
#   ######### Delta_x
#   del_x = sigma
#
#
#   ######### Delta_xyy
#   A= array(apply(Q,2,function(v) as.vector(apply(Q,2,function(vec)
#     cov(x*v,vec)))),c(p,1,4*m^2))
#
#   del_xyy = matrix(NA,2*m*p,2*m*1)
#   for(i in 1:(2*m)){
#     del_xyy[((i-1)*p+1):(i*p),] = as.vector(A[,,((i-1)*(2*m)+1):(i*2*m)])
#   }
#
#
#   ######### Delta_xyx
#   A= array(apply(Q,2,function(v) cov(x*v,x)),c(p,p,2*m))
#
#   del_xyx = matrix(NA,2*m*p,p)
#   for(i in 1:(2*m)){
#     del_xyx[((i-1)*p+1):(i*p),] = as.vector(A[,,i])
#   }
#
#
#   ######### Delta_yx
#   A = array(apply(Q,2,function(v) cov(v,x)),c(1,p,2*m))
#
#   del_yx = matrix(NA,2*m,p)
#   for(i in 1:(2*m)){
#     del_yx[i,] = as.vector(A[,,i])
#   }
#
#
#   ########## Delta #######
#   del = rbind(cbind(del_xy, del_xyy, del_xyx), cbind(t(del_xyy), del_y, del_yx), cbind(t(del_xyx), t(del_yx), del_x))
#
#   ########## A
#   mu = apply(x, 2, mean)
#   mean_Q = apply(Q, 2, mean)
#   Mu = matrix(0, nrow=2*p*m, ncol=2*m)#middle
#   L = matrix(0, nrow=2*p*m, ncol=p)#last
#   for(i in 1:(2*m)){
#     Mu[((i-1)*p+1):(i*p),i] = mu
#     L[((i-1)*p+1):(i*p),] = mean_Q[i]* diag(1,p)
#   }
#   A = cbind(diag(1, 2*p*m), Mu, L)
#   sigmainv = matpower(sigma,-0.5)
#
#   omega = kronecker(t(phi0),t(lam0) %*% sigmainv) %*% A %*% del %*% t(A) %*% kronecker(phi0, sigmainv %*% lam0)
#   eiv = ome_ev = eigen(omega)$values
#   ome_tra = sum(diag(omega))
#   return(list(ome = omega,eigen = ome_ev,trace= ome_tra))
# }
#
# cooktest <- function(omega,stat,B){
#   eiv  = omega$eigen
#   hist = apply(matrix(rchisq(length(eiv)*B,1),ncol=B)*eiv,2,sum)
#   pvalue = mean(hist>stat)
#   return(pvalue)
# }
#
# scaletest <-function(omega,stat,p,k,m){
#   p_sta = ((p-k)*(2*m-k))
#   scale_stat = 1/(omega$trace/p_sta) * stat
#   pchisq(scale_stat,p_sta,lower.tail = F)
# }
#
# addtest <-function(omega,stat){
#   d_sta = (omega$trace)^2/sum(diag(omega$ome%*%omega$ome))
#   adj_stat = 1/(omega$trace/d_sta) * stat
#   pchisq(adj_stat,d_sta,lower.tail = F)
# }
#
# testd <- function(x,y,m){
#   n = dim(x)[1]
#   p = dim(x)[2]
#   w = rnorm(m)
#   pvalue = c()
#   pv=rep(0,3)
#   k=0
#   while(any(pv<0.05)){
#     omega = ome(x,y,k,w)
#     stat = n*sum(mybeta(x,y,p,w)$eigenvalue[(k+1):p])
#     pv = c(cooktest(omega,stat,10000),scaletest(omega,stat,p,k),addtest(omega,stat))
#     pvalue = c(pvalue,pv)
#     k=k+1
#   }
#   return(as.vector(pvalue))
# }

##### standardize x #####

stand <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- cov(x1) # t(x1) %*% (x1)/n
  eva <- eigen(sigma)$values
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1 / sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}
matpower <- function(a, alpha) {
  small <- .000001
  p1 <- nrow(a)
  eva <- eigen(a)$values
  eve <- eigen(a)$vectors
  eve <- eve / t(matrix((diag(t(eve) %*% eve)^0.5), p1, p1))
  index <- (1:p1)[Re(eva) > small]
  evai <- eva
  evai[index] <- (eva[index])^(alpha)
  ai <- eve %*% diag(evai) %*% t(eve)
  return(ai)
}

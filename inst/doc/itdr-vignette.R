## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  install.packages("itdr")
#  library(itdr)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("~/itdr.zip")
#  library(itdr)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(devtools)
#  install_github("TharinduPDeAlwis/itdr")
#  library(itdr)

## ----eval=TRUE, include=TRUE--------------------------------------------------
library(itdr)
data(automobile)
automobile.na <- na.omit(automobile)
# prepare response and predictor variables
auto_y <- log(automobile.na[, 26])
auto_xx <- automobile.na[, c(10, 11, 12, 13, 14, 17, 19, 20, 21, 22, 23, 24, 25)]
auto_x <- scale(auto_xx) # Standardize the predictors
# call to the d.boots() function with required #arguments
d_est <- d.boots(auto_y, auto_x,
  var_plot = TRUE, space = "pdf",
  xdensity = "normal", method = "FM"
)
auto_d <- d_est$d.hat
auto_d

## ----eval=TRUE, include=TRUE--------------------------------------------------
auto_d <- 2 # The estimated value from Section 2.1
candidate_list<-seq(0.05, 1, by = 0.01)
auto_sw2 <- hyperPara(auto_y, auto_x, auto_d, range = candidate_list,xdensity = "normal", B = 500, space = "pdf", method = "FM", hyper = "wx")
auto_sw2$wx.hat # we get the estimator for sw2 as 0.09

## ----eval=TRUE, include=TRUE--------------------------------------------------
auto_d <- 2 # Estimated value from Section 2.1
auto_st2 <- hyperPara(auto_y, auto_x, auto_d, wx = 0.1, range = seq(0.1, 1, by = 0.1), xdensity = "normal", method = "FM", hyper = "wy")
auto_st2$wy.hat # we get the estimator for st2=0.9

## ----eval=TRUE, include=TRUE--------------------------------------------------
h_hat <- hyperPara(auto_y, auto_x, auto_d, wx = 5, wy = 0.1, range = seq(0.1, 2, by = .1), space = "pdf", method = "FM", hyper = "wh")
# Bandwidth estimator for Gaussian kernel density estimation for central subspace
h_hat$h.hat # we have the estimator as h_hat=0.1

## ----eval=TRUE, include=TRUE--------------------------------------------------
library(itdr)
data(automobile)
head(automobile)
df <- cbind(automobile[, c(26, 10, 11, 12, 13, 14, 17, 19, 20, 21, 22, 23, 24, 25)])
dff <- as.matrix(df)
automobi <- dff[complete.cases(dff), ]
d <- 2
# Estimated value from Section 2.1
wx <- .14 # Estimated value from Section 2.2.1
wy <- .9 # Estimated value from Section 2.2.2
wh <- 1.5 # Estimated value from Section 2.2.3
p <- 13 # Estimated value from Section 2.3
y <- automobi[, 1]
x <- automobi[, c(2:14)]
xt <- scale(x)
# Distribution of the predictors is a normal distribution
fit.F_CMS <- itdr(y, xt, d, wx, wy, wh, space = "pdf", xdensity = "normal", method = "FM")
round(fit.F_CMS$eta_hat, 2)

# Distribution of the predictors is a unknown (using kernel method)
fit.F_CMS <- itdr(y, xt, d, wx, wy, wh, space = "pdf", xdensity = "kernel", method = "FM")
round(fit.F_CMS$eta_hat, 2)

## ----eval=TRUE, include=TRUE--------------------------------------------------
library(itdr)
data("recumbent")
recumbent.na <- na.omit(recumbent)
y <- recumbent.na$outcome
X1 <- log(recumbent.na$ast)
X2 <- log(recumbent.na$ck)
X3 <- log(recumbent.na$urea)
p <- 3
x <- matrix(c(X1, X2, X3), ncol = p)
d <- 2
fit.iht_CMS <- itdr(y, x, 2, method = "iht")
fit.iht_CMS$eta_hat

## ----eval=TRUE, include=TRUE--------------------------------------------------
library(itdr)
data(pdb)
colnames(pdb) <- NULL
p <- 15
# select predictor vecotr (y) and response variables (X) according to Weng and Weng and Yin, (2018).
df <- pdb[, c(79, 73, 77, 103, 112, 115, 124, 130, 132, 145, 149, 151, 153, 155, 167, 169)]
dff <- as.matrix(df)
# remove the NA rows
planingdb <- dff[complete.cases(dff), ]

y <- planingdb[, 1] # n-dimensionl response vector
x <- planingdb[, c(2:(p + 1))] # raw desing matrix
x <- x + 0.5
# desing matrix after tranformations
xt <- cbind(
  x[, 1]^(.33), x[, 2]^(.33), x[, 3]^(.57), x[, 4]^(.33), x[, 5]^(.4),
  x[, 6]^(.5), x[, 7]^(.33), x[, 8]^(.16), x[, 9]^(.27), x[, 10]^(.5),
  x[, 11]^(.5), x[, 12]^(.33), x[, 13]^(.06), x[, 14]^(.15), x[, 15]^(.1)
)
m <- 1
W <- sapply(50, rnorm)
# run the hypothsis tests
d.test(y, x, m)

## ----eval=TRUE, include=TRUE--------------------------------------------------
library(itdr)
data(pdb)
colnames(pdb) <- NULL
p <- 15
# select predictor vecotr (y) and response variables (X) according to Weng and Weng and Yin, (2018).
df <- pdb[, c(79, 73, 77, 103, 112, 115, 124, 130, 132, 145, 149, 151, 153, 155, 167, 169)]
dff <- as.matrix(df)
# remove the NA rows
planingdb <- dff[complete.cases(dff), ]

y <- planingdb[, 1] # n-dimensionl response vector
x <- planingdb[, c(2:(p + 1))] # raw desing matrix
x <- x + 0.5
# desing matrix after tranformations give in Weng and Yin, (2018).
xt <- cbind(
  x[, 1]^(.33), x[, 2]^(.33), x[, 3]^(.57), x[, 4]^(.33), x[, 5]^(.4),
  x[, 6]^(.5), x[, 7]^(.33), x[, 8]^(.16), x[, 9]^(.27), x[, 10]^(.5),
  x[, 11]^(.5), x[, 12]^(.33), x[, 13]^(.06), x[, 14]^(.15), x[, 15]^(.1)
)

d <- 1 # estimated dimension of the CS from Section 4.1
invFM.fit <- itdr(y, x, d, m = 50, method = "invFM", x.scale = FALSE) # estimated basis
betahat <- invFM.fit$eta_hat

plot(y ~ xt %*% betahat,
  xlab = "First reduced predictor",
  ylab = "Health insurance coverage"
)

## ----eval=TRUE, include=TRUE--------------------------------------------------
library(itdr)
data(prostate)
X <- as.matrix(prostate[, 1:8])
Y <- matrix(prostate[, 9], ncol = 1)
fit.ftire <- mitdr(X, Y, d = 2, m = 10, method = "FT-IRE")
betahat <- fit.ftire$Beta_hat
betahat
newx <- X %*% betahat
plot(Y ~ newx[, 1],
  xlab = "First reduced predictor",
  ylab = paste0(expression(log), "(antigen)", sep = "")
)
plot(Y ~ newx[, 2],
  xlab = "Second reduced predictor",
  ylab = paste0(expression(log), "(antigen)", sep = "")
)

## ----eval=TRUE, include=TRUE--------------------------------------------------
data(raman)
Y <- as.matrix(raman[, c(1100)]) ## percentage of total fat content
X <- as.matrix(raman[c(2:501)]) ## first 500 wavelength variables
out <- mitdr(X, Y, d = 1, m = 30, method = "admmft", lambda = 0.5, sparse.cov = TRUE, x.scale = TRUE)
estbeta <- out$Beta_hat
estbeta
plot(Y ~ X %*% estbeta,
  xlab = "First reduced predictor",
  ylab = "Percentage of total fat"
)


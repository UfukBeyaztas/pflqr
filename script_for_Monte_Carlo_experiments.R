source("auxilary_functions.R")
source("auxilary_functions_fpca.R")
source("dgp1.R")
source("dgp2.R")

gpx <- (1:50)/50
gpy <- (1:60)/60

# An example run under DGP-I

set.seed(123)
data <- dgp1(n = 50, gpy = gpy, gpx = gpx, sd.error = 0.01)#, out.p = 0.1) # No outlier
Y <- data$Y
X <- data$X

data.test <- dgp1(n = 100, gpy = gpy, gpx = gpx, sd.error = 0.01)
X.test <- data.test$X
Y.test <- data.test$Ytrue

# Proposed method
pffqr.model <- fpenqr(y=Y, x=X, tau=0.5, nby=15, nbx=15, nb0=15, gpy=gpy, gpx=gpx)
pffqr.model.025 <- fpenqr(y=Y, x=X, tau=0.025, nby=15, nbx=15, nb0=15, gpy=gpy, gpx=gpx)
pffqr.model.975 <- fpenqr(y=Y, x=X, tau=0.975, nby=15, nbx=15, nb0=15, gpy=gpy, gpx=gpx)

# FPCA-based function-on-function linear quantile regression model
fpca.model <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.5, gpy=gpy, gpx=list(gpx),
                   fnbasisY=15, fnbasisX=15, fncomp_Y=5, fncomp_X=5)
fpca.model.025 <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.025, gpy=gpy, gpx=list(gpx),
                       fnbasisY=15, fnbasisX=15, fncomp_Y=5, fncomp_X=5)
fpca.model.975 <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.975, gpy=gpy, gpx=list(gpx),
                       fnbasisY=15, fnbasisX=15, fncomp_Y=5, fncomp_X=5)

# pffr method
pffr.model <- pffr(Y ~ ff(X, xind = gpx), yind = gpy, data = data, 
                   bs.int = list(bs = "ps", k = 15, m = c(2, 2)),
                   bs.yindex = list(bs = "ps", k = 15, m = c(2, 2)))


# Predictions for a given test sample
preds.pffr <- predict(pffr.model, newdata = list(X = X.test))
preds.fpca <- fpca.model$pred
preds.pffqr <- predict.fpenqr(object = pffqr.model, xnew = X.test)

preds.pffqr.025 <- predict.fpenqr(object = pffqr.model.025, xnew = X.test)
preds.pffqr.975 <- predict.fpenqr(object = pffqr.model.975, xnew = X.test)

# MSPE values
# pffr
mean((Y.test-preds.pffr)^2) # 0.0006631121
# flqr
mean((Y.test-preds.fpca)^2) # 0.04839337
# pflqr (proposed method)
mean((Y.test-preds.pffqr)^2) # 7.352454e-05


# Regression parameter functions
alpha <- function(t) 2 * exp(-(t - 1)^2)
beta <- function(t, s) {
  4 * cos(2 * pi * t) * sin(pi * s)
}

alpha.t <- alpha(gpy)
beta.ts <- outer(gpy, gpx, beta)

# Estimates of the regression parameter function by pffr
pffrfit <- plot(pffr.model, select = 0)
pffralpha <- pffrfit[[1]]$fit + pffr.model$coefficients[1]
pffralpgrid <- pffrfit[[1]]$x
pffrbetafit <- pffrfit[[2]]$fit
pffrgrids <- pffrfit[[2]]$x
pffrgridt <- pffrfit[[2]]$y
pffrbeta <- t(matrix(pffrbetafit, length(pffrgrids), length(pffrgridt)))
pffralpha.t <- approx(pffralpgrid, pffralpha, gpy)$y 
pffrbeta.ts <- outer(pffrgridt, pffrgrids, beta)

pffqr.alpha <- pffqr.model$b0.hat
pffqr.beta <- pffqr.model$b.hat

# ISE(alpha)
# pffr
mean((alpha.t-pffralpha.t)^2) # 1.148243e-06
# flqr
mean((alpha.t-fpca.model$intercep)^2) # 0.004407539
# pflqr
mean((alpha.t-pffqr.alpha)^2) # 1.365669e-06

# ISE(beta)
# pffr
mean((pffrbeta-pffrbeta.ts)^2) # 0.002531333
# flqr
mean((beta.ts-fpca.model$coef[[1]])^2) # 0.0196923
# pflqr
mean((beta.ts-pffqr.beta)^2) # 0.0003305728

# score and CPD values
sco.cpd.fpca.sim <- matrix(, nrow = 100, ncol = 2)
sco.cpd.pffqr.sim <- matrix(, nrow = 100, ncol = 2)

for(i in 1:100){
  sco.cpd.fpca.sim[i,] <- interval_score(Y.test[i,],fpca.model.025$pred[i,], fpca.model.975$pred[i,], 0.05)
  sco.cpd.pffqr.sim[i,] <- interval_score(Y.test[i,],preds.pffqr.025[i,], preds.pffqr.975[i,], 0.05)
}
# flqr
apply(sco.cpd.fpca.sim, 2, mean) # 7.8014272 0.9246667
# pflqr
apply(sco.cpd.pffqr.sim, 2, mean) # 0.07421137 0.05100000








# An example run under DGP-II

set.seed(123)
data <- dgp2(n = 50, gpy = gpy, gpx = gpx, sd.error = 0.01)#, out.p = 0.1) # No outlier
Y <- data$Y
X <- data$X

data.test <- dgp2(n = 100, gpy = gpy, gpx = gpx, sd.error = 0.01)
X.test <- data.test$X
Y.test <- data.test$Ytrue

# Proposed method
pffqr.model <- fpenqr(y=Y, x=X, tau=0.5, nby=15, nbx=15, nb0=15, gpy=gpy, gpx=gpx)
pffqr.model.025 <- fpenqr(y=Y, x=X, tau=0.025, nby=15, nbx=15, nb0=15, gpy=gpy, gpx=gpx)
pffqr.model.975 <- fpenqr(y=Y, x=X, tau=0.975, nby=15, nbx=15, nb0=15, gpy=gpy, gpx=gpx)

# FPCA-based function-on-function linear quantile regression model
fpca.model <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.5, gpy=gpy, gpx=list(gpx),
                   fnbasisY=15, fnbasisX=15, fncomp_Y=5, fncomp_X=5)
fpca.model.025 <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.025, gpy=gpy, gpx=list(gpx),
                       fnbasisY=15, fnbasisX=15, fncomp_Y=5, fncomp_X=5)
fpca.model.975 <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.975, gpy=gpy, gpx=list(gpx),
                       fnbasisY=15, fnbasisX=15, fncomp_Y=5, fncomp_X=5)

# pffr method
pffr.model <- pffr(Y ~ ff(X, xind = gpx), yind = gpy, data = data, 
                   bs.int = list(bs = "ps", k = 15, m = c(2, 2)),
                   bs.yindex = list(bs = "ps", k = 15, m = c(2, 2)))


# Predictions for a given test sample
preds.pffr <- predict(pffr.model, newdata = list(X = X.test))
preds.fpca <- fpca.model$pred
preds.pffqr <- predict.fpenqr(object = pffqr.model, xnew = X.test)

preds.pffqr.025 <- predict.fpenqr(object = pffqr.model.025, xnew = X.test)
preds.pffqr.975 <- predict.fpenqr(object = pffqr.model.975, xnew = X.test)

# MSPE values
# pffr
mean((Y.test-preds.pffr)^2) # 0.1244713
# flqr
mean((Y.test-preds.fpca)^2) # 0.1378173
# pflqr (proposed method)
mean((Y.test-preds.pffqr)^2) # 0.01161558


# Regression parameter functions
alpha <- function(t) 2 * sin(4*pi*t)
beta <- function(t, s) {
  4 * cos(4 * pi * t) * sin(4*pi * s)
}

alpha.t <- alpha(gpy)
beta.ts <- outer(gpy, gpx, beta)

# Estimates of the regression parameter function by pffr
pffrfit <- plot(pffr.model, select = 0)
pffralpha <- pffrfit[[1]]$fit + pffr.model$coefficients[1]
pffralpgrid <- pffrfit[[1]]$x
pffrbetafit <- pffrfit[[2]]$fit
pffrgrids <- pffrfit[[2]]$x
pffrgridt <- pffrfit[[2]]$y
pffrbeta <- t(matrix(pffrbetafit, length(pffrgrids), length(pffrgridt)))
pffralpha.t <- approx(pffralpgrid, pffralpha, gpy)$y 
pffrbeta.ts <- outer(pffrgridt, pffrgrids, beta)

pffqr.alpha <- pffqr.model$b0.hat
pffqr.beta <- pffqr.model$b.hat

# ISE(alpha)
# pffr
mean((alpha.t-pffralpha.t)^2) # 0.01156349
# flqr
mean((alpha.t-fpca.model$intercep)^2) # 0.06153902
# pflqr
mean((alpha.t-pffqr.alpha)^2) # 0.0003737888

# ISE(beta)
# pffr
mean((pffrbeta-pffrbeta.ts)^2) # 3.817366
# flqr
mean((beta.ts-fpca.model$coef[[1]])^2) # 3.586012
# pflqr
mean((beta.ts-pffqr.beta)^2) # 1.403133

# score and CPD values
sco.cpd.fpca.sim <- matrix(, nrow = 100, ncol = 2)
sco.cpd.pffqr.sim <- matrix(, nrow = 100, ncol = 2)

for(i in 1:100){
  sco.cpd.fpca.sim[i,] <- interval_score(Y.test[i,],fpca.model.025$pred[i,], fpca.model.975$pred[i,], 0.05)
  sco.cpd.pffqr.sim[i,] <- interval_score(Y.test[i,],preds.pffqr.025[i,], preds.pffqr.975[i,], 0.05)
}
# flqr
apply(sco.cpd.fpca.sim, 2, mean) # 17.93953  0.94700
# pflqr
apply(sco.cpd.pffqr.sim, 2, mean) # 1.2187493 0.1718333


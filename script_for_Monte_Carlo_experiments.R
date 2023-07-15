library(here)

source(here("run_order.R"))
# Load the source files.
run_order()

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

# FPLS-based function-on-function linear quantile regression model
X <- matrix(X, ncol = ncol(X))
Y <- matrix(Y, ncol = ncol(Y))
X.test <- matrix(X.test, ncol = ncol(X.test))

# nbx, nby, and h were detedrimen via cross-validation
fpqr.model <- fpqr(y=Y, x=X, xnew=X.test, h=4, tau=0.5, nby=5, nbx=5, gpy=gpy, gpx=gpx)
fpqr.model.025 <- fpqr(y=Y, x=X, xnew=X.test, h=4, tau=0.025, nby=5, nbx=5, gpy=gpy, gpx=gpx)
fpqr.model.975 <- fpqr(y=Y, x=X, xnew=X.test, h=4, tau=0.975, nby=5, nbx=5, gpy=gpy, gpx=gpx)

# Functional partial least squares regression
FVE.x = .99
FVE.y = .99
ridge.par = c(0)
tune.method = 'CV'
p.max = pUpper.compu(X, gpx, FVE.x, basis.name = "bspline")

fplsr.model <- fAPLS(X, Y, gpx, gpy, X.test, ridge.par, p.max, tune = tune.method)

# Functional principal component regression
fpcr.model <- FPCR(X, Y, gpx, gpy, X.test, FVE.x, FVE.y, tune = T)


# Predictions for a given test sample
preds.pffr <- predict(pffr.model, newdata = list(X = X.test))
preds.fpca <- fpca.model$pred
preds.pffqr <- predict.fpenqr(object = pffqr.model, xnew = X.test)
preds.fpqr <- fpqr.model$predicted.values
preds.fplsr <- fplsr.model$y.pred
preds.fpcr <- fpcr.model$y.pred

preds.pffqr.025 <- predict.fpenqr(object = pffqr.model.025, xnew = X.test)
preds.pffqr.975 <- predict.fpenqr(object = pffqr.model.975, xnew = X.test)

preds.fpca.025 <- fpca.model.025$pred
preds.fpca.975 <- fpca.model.975$pred

preds.fpqr.025 <- fpqr.model.025$predicted.values
preds.fpqr.975 <- fpqr.model.975$predicted.values

# MSPE values
# pffr
mean((Y.test-preds.pffr)^2) # 0.0006631121
# flqr
mean((Y.test-preds.fpca)^2) # 0.04839337
# pflqr (proposed method)
mean((Y.test-preds.pffqr)^2) # 7.352454e-05
# fpqr
mean((Y.test-preds.fpqr)^2) # 0.004559364
# fplsr
mean((Y.test-preds.fplsr)^2) # 0.02180835
# fpcr
mean((Y.test-preds.fpcr)^2) # 0.02470341


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

fpqr.alpha <- fpqr.model$b0hat
fpqr.beta <- t(fpqr.model$bhat)

fplsr.beta <- t(fplsr.model$Beta)
fpcr.beta <- t(fpcr.model$Beta)

# ISE(alpha)
# pffr
mean((alpha.t-pffralpha.t)^2) # 1.148243e-06
# flqr
mean((alpha.t-fpca.model$intercep)^2) # 0.004407539
# pflqr
mean((alpha.t-pffqr.alpha)^2) # 1.365669e-06
# fpqr
mean((alpha.t-fpqr.alpha)^2) # 0.00387519

# ISE(beta)
# pffr
mean((pffrbeta-pffrbeta.ts)^2) # 0.002531333
# flqr
mean((beta.ts-fpca.model$coef[[1]])^2) # 0.0196923
# pflqr
mean((beta.ts-pffqr.beta)^2) # 0.0003305728
# fpqr
mean((beta.ts-fpqr.beta)^2) # 0.002433181
# fplsr
mean((beta.ts-fplsr.beta)^2) # 0.0008759661
# fpcr
mean((beta.ts-fpcr.beta)^2) # 0.02142581

# score and CPD values
sco.cpd.fpca.sim <- matrix(, nrow = 100, ncol = 2)
sco.cpd.pffqr.sim <- matrix(, nrow = 100, ncol = 2)
sco.cpd.pfqr.sim <- matrix(, nrow = 100, ncol = 2)

for(i in 1:100){
  sco.cpd.fpca.sim[i,] <- interval_score(Y.test[i,],fpca.model.025$pred[i,], fpca.model.975$pred[i,], 0.05)
  sco.cpd.pffqr.sim[i,] <- interval_score(Y.test[i,],preds.pffqr.025[i,], preds.pffqr.975[i,], 0.05)
  sco.cpd.pfqr.sim[i,] <- interval_score(Y.test[i,],preds.fpqr.025[i,], preds.fpqr.975[i,], 0.05)
}
# flqr
apply(sco.cpd.fpca.sim, 2, mean) # 7.8014272 0.9246667
# pflqr
apply(sco.cpd.pffqr.sim, 2, mean) # 0.07421137 0.05100000
# fpqr
apply(sco.cpd.pfqr.sim, 2, mean) # 1.9980595 0.9103333








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

# FPLS-based function-on-function linear quantile regression model
X <- matrix(X, ncol = ncol(X))
Y <- matrix(Y, ncol = ncol(Y))
X.test <- matrix(X.test, ncol = ncol(X.test))

# nbx, nby, and h were detedrimen via cross-validation
fpqr.model <- fpqr(y=Y, x=X, xnew=X.test, h=4, tau=0.5, nby=10, nbx=10, gpy=gpy, gpx=gpx)
fpqr.model.025 <- fpqr(y=Y, x=X, xnew=X.test, h=4, tau=0.025, nby=10, nbx=10, gpy=gpy, gpx=gpx)
fpqr.model.975 <- fpqr(y=Y, x=X, xnew=X.test, h=4, tau=0.975, nby=10, nbx=10, gpy=gpy, gpx=gpx)

# Functional partial least squares regression
FVE.x = .99
FVE.y = .99
ridge.par = c(0)
tune.method = 'CV'
p.max = pUpper.compu(X, gpx, FVE.x, basis.name = "bspline")

fplsr.model <- fAPLS(X, Y, gpx, gpy, X.test, ridge.par, p.max, tune = tune.method)

# Functional principal component regression
fpcr.model <- FPCR(X, Y, gpx, gpy, X.test, FVE.x, FVE.y, tune = T)


# Predictions for a given test sample
preds.pffr <- predict(pffr.model, newdata = list(X = X.test))
preds.fpca <- fpca.model$pred
preds.pffqr <- predict.fpenqr(object = pffqr.model, xnew = X.test)
preds.fpqr <- fpqr.model$predicted.values
preds.fplsr <- fplsr.model$y.pred
preds.fpcr <- fpcr.model$y.pred

preds.pffqr.025 <- predict.fpenqr(object = pffqr.model.025, xnew = X.test)
preds.pffqr.975 <- predict.fpenqr(object = pffqr.model.975, xnew = X.test)

preds.fpca.025 <- fpca.model.025$pred
preds.fpca.975 <- fpca.model.975$pred

preds.fpqr.025 <- fpqr.model.025$predicted.values
preds.fpqr.975 <- fpqr.model.975$predicted.values

# MSPE values
# pffr
mean((Y.test-preds.pffr)^2) # 0.1602355
# flqr
mean((Y.test-preds.fpca)^2) # 0.1668938
# pflqr (proposed method)
mean((Y.test-preds.pffqr)^2) # 0.05797657
# fpqr
mean((Y.test-preds.fpqr)^2) # 0.1080985
# fplsr
mean((Y.test-preds.fplsr)^2) # 0.1725928
# fpcr
mean((Y.test-preds.fpcr)^2) # 0.2075199


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

fpqr.alpha <- fpqr.model$b0hat
fpqr.beta <- t(fpqr.model$bhat)

fplsr.beta <- t(fplsr.model$Beta)
fpcr.beta <- t(fpcr.model$Beta)

# ISE(alpha)
# pffr
mean((alpha.t-pffralpha.t)^2) # 0.006997925
# flqr
mean((alpha.t-fpca.model$intercep)^2) # 0.06176412
# pflqr
mean((alpha.t-pffqr.alpha)^2) # 0.004339547
# fpqr
mean((alpha.t-fpqr.alpha)^2) # 0.04303604

# ISE(beta)
# pffr
mean((pffrbeta-pffrbeta.ts)^2) # 3.580498
# flqr
mean((beta.ts-fpca.model$coef[[1]])^2) # 3.633358
# pflqr
mean((beta.ts-pffqr.beta)^2) # 2.028303
# fpqr
mean((beta.ts-fpqr.beta)^2) # 2.048421
# fplsr
mean((beta.ts-fplsr.beta)^2) # 3.589995
# fpcr
mean((beta.ts-fpcr.beta)^2) # 3.607519

# score and CPD values
sco.cpd.fpca.sim <- matrix(, nrow = 100, ncol = 2)
sco.cpd.pffqr.sim <- matrix(, nrow = 100, ncol = 2)
sco.cpd.pfqr.sim <- matrix(, nrow = 100, ncol = 2)

for(i in 1:100){
  sco.cpd.fpca.sim[i,] <- interval_score(Y.test[i,],fpca.model.025$pred[i,], fpca.model.975$pred[i,], 0.05)
  sco.cpd.pffqr.sim[i,] <- interval_score(Y.test[i,],preds.pffqr.025[i,], preds.pffqr.975[i,], 0.05)
  sco.cpd.pfqr.sim[i,] <- interval_score(Y.test[i,],preds.fpqr.025[i,], preds.fpqr.975[i,], 0.05)
}
# flqr
apply(sco.cpd.fpca.sim, 2, mean) # 18.1760561  0.7111667
# pflqr
apply(sco.cpd.pffqr.sim, 2, mean) # 6.1430686 0.7383333
# fpqr
apply(sco.cpd.pfqr.sim, 2, mean) # 2.145448 0.097500


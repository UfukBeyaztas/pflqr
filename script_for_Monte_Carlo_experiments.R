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

# RMSPE values
# pffr
rmspe(Y.test, preds.pffr,  gpy)  # 0.7375234
# flqr
rmspe(Y.test, preds.fpca, gpy)  # 10.28026
# pflqr (proposed method)
rmspe(Y.test, preds.pffqr, gpy)  # 0.2452742
# fpqr
rmspe(Y.test, preds.fpqr, gpy)  # 3.062987
# fplsr
rmspe(Y.test, preds.fplsr, gpy)  # 6.923498
# fpcr
rmspe(Y.test, preds.fpcr, gpy)  # 7.175906


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

# RRISPEE(alpha)
# pffr
reispe.a(alpha.t, pffralpha.t, gpy)  # 0.06648273
# flqr
reispe.a(alpha.t, fpca.model$intercep, gpy)  # 4.226501
# pflqr
reispe.a(alpha.t, pffqr.alpha, gpy) # 0.07477094
# fpqr
reispe.a(alpha.t, fpqr.alpha, gpy)  # 3.964192

# RRISPEE(beta)
# pffr
rispe.b(pffrbeta, pffrbeta.ts, pffrgridt, pffrgrids) # 2.422286
# flqr
rispe.b(beta.ts, fpca.model$coef[[1]], gpx, gpy) # 6.841824
# pflqr
rispe.b(beta.ts, pffqr.beta, gpx, gpy) # 0.8054657
# fpqr
rispe.b(beta.ts, fpqr.beta, gpx, gpy) # 2.274095
# fplsr
rispe.b(beta.ts, fplsr.beta, gpx, gpy) # 1.453825
# fpcr
rispe.b(beta.ts, fpcr.beta, gpx, gpy) # 7.146398

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

# RMSPE values
# pffr
rmspe(Y.test, preds.pffr,  gpy)  # 21.60608
# flqr
rmspe(Y.test, preds.fpca, gpy)  # 22.2455
# pflqr (proposed method)
rmspe(Y.test, preds.pffqr, gpy)  # 13.69135
# fpqr
rmspe(Y.test, preds.fpqr, gpy)  # 18.85225
# fplsr
rmspe(Y.test, preds.fplsr, gpy)  # 24.14197
# fpcr
rmspe(Y.test, preds.fpcr, gpy)  # 25.0166


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

# RRISPEE(alpha)
# pffr
reispe.a(alpha.t, pffralpha.t, gpy)  # 5.855526
# flqr
reispe.a(alpha.t, fpca.model$intercep, gpy)  # 16.87415
# pflqr
reispe.a(alpha.t, pffqr.alpha, gpy) # 4.450031
# fpqr
reispe.a(alpha.t, fpqr.alpha, gpy)  # 14.15366

# RRISPEE(beta)
# pffr
rispe.b(pffrbeta, pffrbeta.ts, pffrgridt, pffrgrids) # 340.1204
# flqr
rispe.b(beta.ts, fpca.model$coef[[1]], gpx, gpy) # 95.26373
# pflqr
rispe.b(beta.ts, pffqr.beta, gpx, gpy) # 70.40651
# fpqr
rispe.b(beta.ts, fpqr.beta, gpx, gpy) # 71.1427
# fplsr
rispe.b(beta.ts, fplsr.beta, gpx, gpy) # 94.68226
# fpcr
rispe.b(beta.ts, fpcr.beta, gpx, gpy) # 94.91387

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






# An example run under DGP-III

set.seed(123)
data <- dgp3(n = 50, gpy = gpy, gpx = gpx)
Y <- data$Y
X <- data$X

data.test <- dgp3(n = 100, gpy = gpy, gpx = gpx)
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

# RMSPE values
# pffr
rmspe(Y.test, preds.pffr,  gpy)  # 37.60322
# flqr
rmspe(Y.test, preds.fpca, gpy)  # 42.14422
# pflqr (proposed method)
rmspe(Y.test, preds.pffqr, gpy)  # 32.55478
# fpqr
rmspe(Y.test, preds.fpqr, gpy)  # 35.00295
# fplsr
rmspe(Y.test, preds.fplsr, gpy)  # 68.37568
# fpcr
rmspe(Y.test, preds.fpcr, gpy)  # 69.03762


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

# RRISPEE(alpha)
# pffr
reispe.a(alpha.t, pffralpha.t, gpy)  # 32.8456
# flqr
reispe.a(alpha.t, fpca.model$intercep, gpy)  # 40.64746
# pflqr
reispe.a(alpha.t, pffqr.alpha, gpy) # 30.75653
# fpqr
reispe.a(alpha.t, fpqr.alpha, gpy)  # 33.56794

# RRISPEE(beta)
# pffr
rispe.b(pffrbeta, pffrbeta.ts, pffrgridt, pffrgrids) # 344.1322
# flqr
rispe.b(beta.ts, fpca.model$coef[[1]], gpx, gpy) # 95.23561
# pflqr
rispe.b(beta.ts, pffqr.beta, gpx, gpy) # 71.52271
# fpqr
rispe.b(beta.ts, fpqr.beta, gpx, gpy) # 77.46408
# fplsr
rispe.b(beta.ts, fplsr.beta, gpx, gpy) # 94.82672
# fpcr
rispe.b(beta.ts, fpcr.beta, gpx, gpy) # 94.92413

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
apply(sco.cpd.fpca.sim, 2, mean) # 27.56852  0.85300
# pflqr
apply(sco.cpd.pffqr.sim, 2, mean) # 16.01658 0.95000
# fpqr
apply(sco.cpd.pfqr.sim, 2, mean) # 4.9388578 0.1783333

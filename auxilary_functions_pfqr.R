# Packages
library(Matrix)
library(fda)
library(quantreg)
library(refund)
library(nloptr)
library(goffda)
library(expm)

# Function to estimate function-on-function penalized quantile regression
fpenqr <- function(y, x, tau, nby, nbx, nb0, gpy, gpx, alpha = 0.005)
{
  lambda.cand <- c(10^{-3}, 10^{-2}, 10^{-1}, 1, 5)
  
  qbic.ij <- matrix(NA, nrow = length(lambda.cand)^2, 3)
  qbic.ij[,1] <- rep(1:length(lambda.cand), each = length(lambda.cand))
  qbic.ij[,2] <- rep(1:length(lambda.cand), length(lambda.cand))
  
  ijk = 1
  for(i in 1:length(lambda.cand)){
    try({
      for(j in 1:length(lambda.cand)){
        qbic.ij[ijk, 3] <- qbic(y, x, tau, nby, nbx, nb0, gpy, gpx,
                                l0=lambda.cand[i], l1=lambda.cand[j], alpha)
        ijk = ijk + 1
      }
    }, silent = T)
  }
  
  l0 <- lambda.cand[qbic.ij[,1][which.min(qbic.ij[,3])]]
  l1 <- lambda.cand[qbic.ij[,2][which.min(qbic.ij[,3])]]
  
  final.model <- penqr.lam(y, x, tau, nby, nbx, nb0, gpy, gpx, l0, l1, alpha)
  fb0.hat <- final.model$b0.hat
  fb.hat <- final.model$b.hat
  fy.hat <- final.model$y.hat
  f.resids <- final.model$resids
  
  return(list(b0.hat = fb0.hat, b.hat = fb.hat, y.hat = fy.hat, resids = f.resids,
              gpy = gpy, gpx = gpx, l0=l0, l1=l1))
}

# Function to obtain predictions for a given new set of functional predictor
predict.fpenqr <- function(object, xnew)
{
  
  gpx <- object$gpx
  b0.hat <- object$b0.hat
  b.hat <- object$b.hat
  
  y.pred <- xnew %*% t(b.hat) * (gpx[2] - gpx[1])
  for(i in 1:dim(y.pred)[1])
    y.pred[i,] <- y.pred[i,] + c(b0.hat)
  
  return(y.pred)
}

# A sub-function used to estimate function-on-function penalized quantile regression 
penqr.lam <- function(y, x, tau, nby, nbx, nb0, gpy, gpx, l0, l1, alpha)
{
  
  #y <- smooth_fun(y, argvals = gpy)
  #x <- smooth_fun(x, argvals = gpx)
  
  n <- dim(y)[1]
  py <- dim(y)[2]
  px <- dim(x)[2]
  y0 <- y
  y <- c(t(y))
  
  bs_basis0 <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = nb0)
  bs_basisy <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = nby)
  bs_basisx <- create.bspline.basis(rangeval = c(gpx[1], gpx[px]), nbasis = nbx)
  evalbase0 = eval.basis(gpy, bs_basis0)
  evalbasey = eval.basis(gpy, bs_basisy)
  evalbasex = eval.basis(gpx, bs_basisx)
  diff.x <- gpx[2] - gpx[1]
  
  x.p <- (x %*% evalbasex) * diff.x
  x0.k <- kronecker(rep(1, n), evalbase0)
  x1.k <- kronecker(x.p, evalbasey)
  x.k <- cbind(x0.k, x1.k)
  
  pm.b0 <- bsplinepen(bs_basis0, Lfdobj = 2)
  pm.b.t0 <- bsplinepen(bs_basisy, Lfdobj = 0)
  pm.b.s0 <- bsplinepen(bs_basisx, Lfdobj = 0)
  pm.b.t2 <- bsplinepen(bs_basisy, Lfdobj = 2)
  pm.b.s2 <- bsplinepen(bs_basisx, Lfdobj = 2)
  
  pm.bt <- kronecker(pm.b.s0, pm.b.t2)
  pm.bs <- kronecker(pm.b.s2, pm.b.t0)
  pen.mat <- as.matrix(bdiag(l0*pm.b0, l1*(pm.bt+pm.bs)))
  
  bhats <- gqlm(y, x.k, tau, pen.mat, alpha, n)
  fb0.hat <- evalbase0 %*% bhats[1:nb0]
  fb.hat <- evalbasey %*% matrix(bhats[-(1:nb0)], nrow = nby, ncol = nbx) %*%
    t(evalbasex)
  fy.hat <- matrix(c(x.k %*% bhats), nrow = n, byrow = TRUE)
  f.resids <- y0 - fy.hat
  
  return(list(b0.hat = fb0.hat, b.hat = fb.hat, y.hat = fy.hat, resids = f.resids))
  
}

# Check-loss function
check.loss <- function(u, tau)
{
  u * (tau - (u < 0))
}

# Approximation of the check-loss function (Zheng, 2011)
mcheck.loss <- function(u, tau, alpha)
{
  u * tau + alpha * log(1 + exp(-u/alpha))
}

# Function to estimate parameters of penalized linear regression
qr.coef.lpmat <- function(y, x, pen.mat)
{
  xtx <- crossprod(x)
  xty <- crossprod(x, y)
  arg <- xtx + pen.mat
  tmp <- qr.coef(qr(arg), xty)
  return(tmp)
}

# Function to estimate parameters of penalized quantile regression
# via derivative-free bound-constrained optimization BOBYQA
gqlm <- function(y, x, tau, pen.mat, alpha, n)
{
  
  sobj_fun <- function(b, tau)
  {
    err <- y - x %*% b
    err2 <- matrix(err, nrow = n, byrow = T)
    buff <- apply(err2^2, 1, sum)
    trim_index <- sort.int(buff, decreasing = FALSE,
                          index.return = TRUE)
    ntrim <- round(0.8 * n)
    index_trunc <- trim_index$ix[1:ntrim]
    newerr <- err2[index_trunc,]
    newerr <- c(t(newerr))
    sum(mcheck.loss(newerr, tau, alpha)) + t(b) %*% pen.mat %*% b
  }
  
  binit <- qr.coef.lpmat(y, x, pen.mat)
  tmp <- bobyqa(x0 = binit, fn = sobj_fun, tau = tau)
  return(tmp$par)
}


# Function to compute BIC for selection of smoothing parameters
qbic <- function(y, x, tau, nby, nbx, nb0, gpy, gpx, l0, l1, alpha)
{
  model.bic <- penqr.lam(y, x, tau, nby, nbx, nb0, gpy, gpx, l0, l1, alpha)
  fit.bic <- model.bic$y.hat
  res.bic <- y - fit.bic
  
  res.sq <- apply(res.bic^2, 1, sum)
  trim_index <- sort.int(res.sq, decreasing = FALSE,
                        index.return = TRUE)
  ntrim <- round(0.8 * length(res.sq))
  index_trunc <- trim_index$ix[1:ntrim]
  newerr <- res.bic[index_trunc,]
  
  BIC_val <- log(sum(mcheck.loss(newerr, tau, alpha))) + log(nrow(newerr))
  
  return(BIC_val)
}

# Function to smooth raw data
smooth_fun <- function(data, argvals = NULL){
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  if(is.null(argvals))
    argvals <- seq(0, 1, length.out = p)
  nbasis <- min(40, round(p/5))
  basis.obs <- create.bspline.basis(range(argvals), nbasis)
  
  sdata <- matrix(, nrow = n, ncol = p)
  for(i in 1:n){
    xs <- smooth.basis(argvals = argvals, y= c(data[i,]), fdParobj = basis.obs)
    xfd <- xs$fd
    sdata[i,] <- eval.fd(argvals, xfd)
  }
  
  return(sdata)
}


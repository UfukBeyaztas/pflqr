
svd_flip <- function(x,y)
{
  abs_idx <- which.max(abs(x))
  sgn <- sign((x[abs_idx]))
  x <- x*sgn
  y <- y*sgn
  
  return(list(x = x, y = y))
}

# Scale function
scale_fun <- function(data)
{
  n <- dim(data)[1]
  p <- dim(data)[2]
  
    cen.data <- apply(data, 2, mean)
  
  scl.data <- (data - matrix(cen.data, nrow = n, ncol = p, byrow = TRUE))
  
  return(scl.data)
}

# Quantile covariance function
qcov <- function(y, x, tau)
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  cf <- numeric()
  for(i in 1:p){
    mod.i <- rq(y~x[,i], tau)
    cf[i] <- mod.i$coefficients[2]
  }
  
  output <- diag(cov(x)) * cf
  
  return(output)
}

# Function to obtain basis expansion coefficients
getAmat <- function(data, nbf, gp)
{
  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  bs_basis <- create.bspline.basis(rangeval = c(gp[1], gp[p]), nbasis = nbf)
  inp_mat <- inprod(bs_basis, bs_basis)
  sinp_mat <- sqrtm(inp_mat)
  evalbase = eval.basis(gp, bs_basis)
  fdobj <- fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  Amat <- t(pcaobj$coefs) %*% sinp_mat
  
  return(list(Amat = Amat, bs_basis = bs_basis, inp_mat = inp_mat,
              sinp_mat = sinp_mat, evalbase = evalbase))
}

# Partial quantile regression
pqr <- function(y, x, tau, h)
{
  x0 <- x
  y0 <- y
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  q <- dim(y)[2]
  
  x <- scale(x, scale = F)
  y <- scale(y, scale = F)
  
  xw <- matrix(, p, h)
  yw <- matrix(, q, h)
  xs <- matrix(, n, h)
  ys <- matrix(, n, h)
  xl <- matrix(, p, h)
  yl <- matrix(, q, h)
  
  for(i in 1:h)
  {
    qcv <- matrix(, p, q)
    for(qi in 1:q)
      qcv[,qi] <- qcov(y[,qi], x, tau)
    qsvd <- svd(qcv)
    
    xw. <- qsvd$u[,1]
    yw. <- qsvd$v[1,]
    
    buff <- svd_flip(xw., yw.)
    xw. <- buff$x
    yw. <- buff$y
    
    xs. <- x %*% xw.
    ys. <- y %*% yw.
    
    xl. <- c(t(x) %*% xs.) / c(t(xs.) %*% xs.)
    x <- x - xs. %*% xl.
    yl. <- c(t(xs.) %*% y) / c(t(xs.) %*% xs.)
    y <- y - xs. %*% yl.
    
    xw[,i] <- xw.
    yw[,i] <- yw.
    xs[,i] <- xs.
    ys[,i] <- ys.
    xl[,i] <- xl.
    yl[,i] <- yl.
  }
  
  xr <- xw %*% ginv(t(xl) %*% xw)
  yr <- yw %*% ginv(t(yl) %*% yw)
  
  fin.cf <- matrix(, (p+1), q)
  pqr.cf <- matrix(, (h+1), q)
  for(ic in 1:q){
    model.ic <- rq(y0[,ic]~xs, tau)
    temp <- xr %*% as.matrix(model.ic$coefficients[-1])
    pqr.cf[,ic] <- model.ic$coefficients
    fin.cf[,ic] <- c(model.ic$coefficients[1], temp)
  }
  
  return(list(T = xs, R = xr, P = xl, W = xw, d.coef = fin.cf,
              pqr.coef = pqr.cf))
}

# Functional partial quantile regression
fpqr <- function(y, x, xnew, h, tau, nby, nbx, gpy, gpx)
{
  BS.sol.x <- getAmat(data = x, nbf = nbx, gp = gpx)
  x1 <- BS.sol.x$Amat
  BS.sol.y <- getAmat(data = y, nbf = nby, gp = gpy)
  y1 <- BS.sol.y$Amat
  
  Bs.sol.xnew <- getAmat(data = xnew, nbf = nbx, gp = gpx)
  x1new <- Bs.sol.xnew$Amat
  
  m.pqr <- pqr(y = y1, x = x1, tau = tau, h = h)
  
  fits <- (cbind(1,x1) %*% m.pqr$d.coef) %*% solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)
  preds <- (cbind(1,x1new) %*% m.pqr$d.coef) %*% solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)
  
  b.hat <- BS.sol.x$evalbase %*% solve(BS.sol.x$sinp_mat) %*% m.pqr$d.coef[-1,] %*% 
    solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)
  
  b0.hat <- t(solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)) %*% as.matrix(m.pqr$d.coef[1,])
  
  return(list(fitted.values = fits, predicted.values = preds, bhat = b.hat, b0hat = b0.hat))
  
}

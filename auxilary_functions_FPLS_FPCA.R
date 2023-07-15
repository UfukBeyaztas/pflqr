# Auxilary functions to perform standard functional partial least squares
# and functional principal component regression models.
# Codes were downloaded from:
# https://github.com/ZhiyangGeeZhou/fAPLS
# Zhiyang Zhou (2021), Fast implementation of partial least squares for function-on-function regression,
# Journal of Multivariate Analysis, 185, 104769


#' type = 
#' 100: \int f(t)dt
#' 201: \int f(s,t)dt
#' 211: \int f(s,t)g(t)dt
#' 222: \int f(s,w)g(w,t)dw

options(warn=-1)

integral = function(f, g = NULL, domain, type){
  if (type == 100){
    f = as.vector(f)
    len.f = length(f)
    result = sum((f[-1] + f[-len.f]) * diff(domain))/2
    return(result)
  }
  
  if (type == 201){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = matrix((rowSums(f[, -1] * gap.mat) + rowSums(f[, -ncol.f] * gap.mat))/2, nrow = nrow.f, ncol = 1)
    return(result)
  }
  
  if (type == 211){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = ((f[, -1] * gap.mat) %*% as.matrix(g[-1]) + (f[, -ncol.f] * gap.mat) %*% as.matrix(g[-ncol.f]))/2
    return(as.matrix(result))
  }
  
  if (type == 222){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = ((f[, -1] * gap.mat) %*% g[-1, ] + (f[, -ncol.f] * gap.mat) %*% g[-ncol.f, ])/2
    return(result)
  }
}

pUpper.compu = function(x, domain.x, proportion, basis.name){
  mu.x = colMeans(x)
  expansion.x = expand.basis(x, domain.x, basis.name = basis.name)
  W = inprod(expansion.x$basis, expansion.x$basis)
  # x.fd = smooth.basisPar(argvals=1:ncol(x),
  #                        y = t(sweep(x, 2, mu.x)),
  #                        fdobj=expansion.x$basis,
  #                        Lfdobj=int2Lfd(2),
  #                        lambda=expansion.x$lambda)
  x.pca.obj = pca.fd(expansion.x$fdObj, nharm = min(dim(x) - 1), centerfns = T)
  pMax = which.max(cumsum(x.pca.obj$values)/sum(x.pca.obj$values) >= proportion)
  return(pMax)
}

expand.basis = function(x, domain.x, nOrder = 4, basis.name){
  if (!("fda" %in% rownames(installed.packages()))) install.packages("fda")
  if (!all(c(
    exists("create.bspline.basis", mode="function"),
    exists("int2Lfd", mode="function"),
    exists("fdPar", mode="function"),
    exists("Data2fd", mode="function"),
    exists("inprod", mode="function"),
    exists('deriv.fd', mode='function')
  ))) library("fda")
  
  timePts = domain.x
  K = min(nOrder+ncol(x)-2, nrow(x)-1) # number of basis functions
  if (basis.name == "bspline")
    basis = create.bspline.basis(rangeval=range(timePts), nbasis=K, norder = nOrder, names="bspl")
  if (basis.name == "fourier")
    basis = create.fourier.basis(rangeval = c(0,1), nbasis = K)
  D2Lfd = int2Lfd(m=2)
  
  # tune the value of lambda for smoothing x
  log10lambda = (-10):10
  gcvsave = NULL
  for (i in log10lambda) {
    lambda = 10^i
    D2fdPar = fdPar(basis, D2Lfd, lambda)
    smooth = smooth.basis(y=t(x), argvals=timePts, D2fdPar, dfscale=1.2)
    gcvsave = c(gcvsave, sum(smooth$gcv))
  }
  
  lambda.opt = 10^log10lambda[which.min(gcvsave)]
  D2fdPar = fdPar(basis, D2Lfd, lambda = lambda.opt)
  fd.smooth = smooth.basis(y=t(x), argvals=timePts, D2fdPar, dfscale=1.2)
  
  return(list(fdObj = fd.smooth$fd, coef = t(fd.smooth$fd$coefs), basis = basis, D2fdPar = D2fdPar, lambda = lambda))
} 

# output basis (as a list) of a p-dimensional Krylov subspace
KS = function(p, r.xx, r.xy, domain.x){
  basis = array(NA, dim = c(dim(r.xy), p))
  for (i in 1:p){
    if (i==1){
      basis[, , 1] = r.xy
    }else
      basis[, , i] = integral(r.xx, basis[, , i-1], domain.x, type = 222)
  }
  return(basis)
}

# Gram-Schmidt orthonormalization w.r.t. ker
# basis.origi is a list
GSortho = function(basis.origi, ker, domain.x, domain.y = NULL){
  
  basis.ortho = array(NA, dim = dim(basis.origi))
  
  if (is.null(domain.y)){
    p = dim(basis.origi)[2]
    for (i in 1:p){
      psi.tmp = as.matrix(basis.origi[, i])
      if (i > 1){
        for (j in 1:(i-1)){
          integ1 = integral(ker, psi.tmp, domain.x, type = 211)
          integ2 = integral(integ1 * as.matrix(basis.ortho[, j]), domain = domain.x, type = 100)
          psi.tmp = psi.tmp - as.matrix(basis.ortho[, j]) * integ2
        }
      }
      integ1 = integral(ker, psi.tmp, domain.x, type = 211)
      integ2 = integral(integ1 * psi.tmp, domain = domain.x, type = 100)
      basis.ortho[, i] = psi.tmp / integ2^.5
    }
  }
  
  if (!is.null(domain.y)){
    p = dim(basis.origi)[3]
    for (i in 1:p){
      psi.tmp = basis.origi[, , i]
      if (i > 1){
        for (j in 1:(i-1)){
          integ1 = integral(ker, psi.tmp, domain.x, type = 222)
          integ2 = integral(integ1 * basis.ortho[, , j], domain = domain.y, type = 201)
          integ3 = integral(integ2, domain = domain.x, type = 100)
          psi.tmp = psi.tmp - basis.ortho[, , j] * integ3
        }
      }
      integ1 = integral(ker, psi.tmp, domain.x, type = 222)
      integ2 = integral(integ1 * psi.tmp, domain = domain.y, type = 201)
      integ3 = integral(integ2, domain = domain.x, type = 100)
      basis.ortho[, , i] = psi.tmp / integ3^.5
    }
  }
  
  return(basis.ortho)
}

test.ortho = function(basis.ortho, ker, domain.x, domain.y = NULL){
  if (is.null(domain.y) & is.null(ker)){
    p = ncol(basis.ortho)
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        H[i, j] = integral(basis.ortho[, i] * basis.ortho[, j], domain = domain.x, type = 100)
      }
    }
  }
  
  if (is.null(domain.y) & !is.null(ker)){
    p = ncol(basis.ortho)
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ = integral(ker, basis.ortho[, j], domain.x, type = 211)
        H[i, j] = integral(basis.ortho[, i] * integ, domain = domain.x, type = 100)
      }
    }
  }
  
  if (!is.null(domain.y) & is.null(ker)){
    p = dim(basis.ortho)[3]
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ = integral(basis.ortho[, , i] * basis.ortho[, , j], domain = domain.y, type = 201)
        H[i, j] = integral(integ, domain = domain.x, type = 100)
      }
    }
  }
  
  if (!is.null(domain.y) & !is.null(ker)){
    p = dim(basis.ortho)[3]
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ1 = integral(ker, basis.ortho[, , j], domain.x, type = 222)
        integ2 = integral(basis.ortho[, , i] * integ1, domain = domain.y, type = 201)
        H[i, j] = integral(integ2, domain = domain.x, type = 100)
      }
    }
  }
  
  return(H)
}

# fAPLS
fAPLS = function(x.old, y.old, domain.x, domain.y, x.new = NULL, ridge.par, p.max, tune = 'CV'){
  
  if (length(ridge.par) > 1)
    hyper.par = cbind(ridge.par, sample.int(p.max, size = length(ridge.par), replace = T))
  else 
    hyper.par = cbind(ridge.par, p.max)
  
  if (!is.null(x.new)){
    mu.x = colMeans(rbind(x.old, x.new))
  }else{
    mu.x = colMeans(x.old) 
  }
  mu.y = colMeans(y.old)
  n = nrow(x.old)
  r.xx = cov(x.old)
  r.xy = cov(x.old, y.old)
  
  if (tune == 'GCV'){
    for (i in 1:length(ridge.par)){
      basis.origi = KS(p.max, r.xx + ridge.par[i] * diag(ncol(r.xx)), r.xy, domain.x)
      basis.ortho = GSortho(basis.origi, r.xx, domain.x, domain.y)
      gamma.vec = numeric(hyper.par[i, 2])
      for (p in 1:hyper.par[i, 2]){
        integ = integral(r.xy * basis.ortho[, , p], domain = domain.y, type = 201)
        gamma.vec[p] = integral(integ, domain = domain.x, type = 100)
        if (p == 1){
          Beta.tilde.curr = gamma.vec[p] * basis.ortho[, , p]
          y.old.tilde.curr = matrix(mu.y, nrow = n, ncol = length(mu.y), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.tilde.curr, domain.x, type = 222)
          Beta.tilde.opti = Beta.tilde.curr
          y.old.tilde.opti = y.old.tilde.curr
          p.opti = p
          lambda.opti = ridge.par[i]
        }else{
          Beta.tilde.curr = Beta.tilde.curr + gamma.vec[p] * basis.ortho[, , p]
          y.old.tilde.curr = matrix(mu.y, nrow = n, ncol = length(mu.y), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.tilde.curr, domain.x, type = 222)
          if (
            sum(integral((y.old.tilde.curr - y.old)^2, domain = domain.y, type = 201))/(n - p - 1)^2 <
            sum(integral((y.old.tilde.opti - y.old)^2, domain = domain.y, type = 201))/(n - p.opti - 1)^2
          ){
            Beta.tilde.opti = Beta.tilde.curr
            y.old.tilde.opti = y.old.tilde.curr
            p.opti = p
            lambda.opti = ridge.par[i]
          }
        }
      }
    }
  }
  
  if (tune == 'CV'){
    nfold = 5
    cv.array = array(NA, dim = c(p.max, length(ridge.par), nfold))
    groups <- split(sample(1:n), rep(1:nfold, length = n))
    for (k in 1:nfold){
      x.train = x.old[-groups[[k]], ]
      y.train = y.old[-groups[[k]], ]
      x.vali = x.old[groups[[k]], ]
      y.vali = y.old[groups[[k]], ]
      mu.x.train = colMeans(x.train)
      mu.y.train = colMeans(y.train)
      
      r.xx.train = cov(x.train)
      r.xy.train = cov(x.train, y.train)
      
      for (i in 1:length(ridge.par)){
        basis.origi = KS(p.max, r.xx.train + ridge.par[i]*diag(ncol(r.xx.train)), r.xy.train, domain.x)
        basis.ortho = GSortho(basis.origi, r.xx.train, domain.x, domain.y)
        gamma.vec = numeric(hyper.par[i, 2])
        for (p in 1:hyper.par[i, 2]){
          integ = integral(r.xy.train * basis.ortho[, , p], domain = domain.y, type = 201)
          gamma.vec[p] = integral(integ, domain = domain.x, type = 100)
          if (p == 1){
            Beta.tilde.curr = gamma.vec[p] * basis.ortho[, , p]
          }else{
            Beta.tilde.curr = Beta.tilde.curr + gamma.vec[p] * basis.ortho[, , p]
          }
          y.vali.tilde.curr = matrix(mu.y.train, nrow = nrow(x.vali), ncol = length(mu.y.train), byrow = T) +
            integral(sweep(x.vali, 2, mu.x.train), Beta.tilde.curr, domain.x, type = 222)
          cv.array[p, i, k] =
            sum(integral((y.vali.tilde.curr - y.vali)^2, domain = domain.y, type = 201))/
            sum(integral(sweep(y.vali, 2, mu.y.train)^2, domain = domain.y, type = 201))
          # cv.array[p, i, k] = 
          #   mean(integral((y.vali.tilde.curr - y.vali)^2, domain = domain.y, type = 201))
        }
      }
    }
    cv.mat = apply(cv.array, c(1, 2), mean)
    loc.opt = which(cv.mat == min(cv.mat, na.rm = T), arr.ind = T)
    lambda.opti = ridge.par[loc.opt[1, 2]]
    p.opti = loc.opt[1, 1]
    
    basis.origi = KS(p.max, r.xx + lambda.opti*diag(ncol(r.xx)), r.xy, domain.x)
    basis.ortho = GSortho(basis.origi, r.xx, domain.x, domain.y)
    gamma.vec = numeric(p.opti)
    for (p in 1:p.opti){
      integ = integral(r.xy * basis.ortho[, , p], domain = domain.y, type = 201)
      gamma.vec[p] = integral(integ, domain = domain.x, type = 100)
      if (p == 1){
        Beta.tilde.opti = gamma.vec[p] * basis.ortho[, , p]
      }else{
        Beta.tilde.opti = Beta.tilde.opti + gamma.vec[p] * basis.ortho[, , p]
      }
    }
  }
  
  y.new.tilde.opti = NULL
  if (!is.null(x.new)){
    y.new.tilde.opti = matrix(mu.y, nrow = nrow(x.new), ncol = length(mu.y), byrow = T) +
      integral(sweep(x.new, 2, mu.x), Beta.tilde.opti, domain.x, type = 222)
  }
  return(list(Beta = Beta.tilde.opti, lambda = lambda.opti, ncomp = p.opti, y.pred = y.new.tilde.opti))
}

# FPC regression
FPCR = function(x.old, y.old, domain.x, domain.y, x.new = NULL, FVE.x, FVE.y, tune = T){
  expansion.x = expand.basis(x.old, domain.x, basis.name = "bspline")
  x.pca.obj = pca.fd(expansion.x$fdObj, nharm = min(dim(x.old) - 1))
  p.max.x = which.max(cumsum(x.pca.obj$varprop) >= FVE.x) 
  eigen.fun.x = t(eval.basis(evalarg = domain.x, basisobj = x.pca.obj$harmonics$basis) %*% x.pca.obj$harmonics$coefs[, 1:p.max.x])
  
  expansion.y = expand.basis(y.old, domain.y, basis.name = "bspline")
  y.pca.obj = pca.fd(expansion.y$fdObj, nharm = min(dim(y.old) - 1))
  p.max.y = which.max(cumsum(y.pca.obj$varprop) >= FVE.y) 
  eigen.fun.y = t(eval.basis(evalarg = domain.y, basisobj = y.pca.obj$harmonics$basis) %*% y.pca.obj$harmonics$coefs[, 1:p.max.y])
  
  if (!is.null(x.new)){
    mu.x = colMeans(rbind(x.old, x.new))
  }else{
    mu.x = colMeans(x.old)
  }
  mu.y = colMeans(y.old)
  
  r.xx = cov(x.old)
  r.xy = cov(x.old, y.old)
  n = nrow(x.old)
  
  Coef.all = array(NA, dim = c(dim(r.xy), p.max.x, p.max.y))
  for (p.x in 1:p.max.x){
    for (p.y in 1:p.max.y){
      integ = integral(r.xy, eigen.fun.y[p.y, ], domain.y, type = 211)
      Coef.all[, , p.x, p.y] = integral(eigen.fun.x[p.x, ] * integ, domain = domain.x, type = 100)/
        x.pca.obj$values[p.x] *
        as.matrix(eigen.fun.x[p.x, ]) %*% t(as.matrix(eigen.fun.y[p.y, ]))
      if (tune == F & p.x == p.max.x & p.y == p.max.y){
        Beta.hat.opti = rowSums(Coef.all, dims = 2)
        p.opti.x = p.max.x
        p.opti.y = p.max.y
      }
      if (tune == T){
        if (p.x == 1 & p.y == 1){
          Beta.hat.opti = Coef.all[, , 1, 1]
          y.old.hat.opti = matrix(mu.y, ncol = length(mu.y), nrow = nrow(x.old), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.hat.opti, domain.x, type = 222)
          p.opti.x = p.x
          p.opti.y = p.y
        }else{
          Beta.hat.curr = rowSums(Coef.all[, , 1:p.x, 1:p.y], dims = 2)
          y.old.hat.curr = matrix(mu.y, ncol = length(mu.y), nrow = nrow(x.old), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.hat.curr, domain.x, type = 222)
          if (
            sum(integral((y.old.hat.curr - y.old)^2, domain = domain.y, type = 201))/(n - max(p.x, p.y) - 1)^2 <
            sum(integral((y.old.hat.opti - y.old)^2, domain = domain.y, type = 201))/(n - max(p.opti.x, p.opti.y) - 1)^2
          ){
            p.opti.x = p.x
            p.opti.y = p.y
            Beta.hat.opti = Beta.hat.curr
            y.old.hat.opti = y.old.hat.curr
          }
        }
      }
    }
  }
  
  y.new.hat.opti = NULL
  reispe.opti = NULL
  if (!is.null(x.new)){
    y.new.hat.opti = matrix(mu.y, ncol = length(mu.y), nrow = nrow(x.new), byrow = T) +
      integral(sweep(x.new, 2, mu.x), Beta.hat.opti, domain.x, type = 222)
  }
  return(list(Beta = Beta.hat.opti, ncomp = c(p.opti.x, p.opti.y), y.pred = y.new.hat.opti))
}

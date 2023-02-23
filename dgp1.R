X.phi1 <- function(nphi, gpx) {
  phi <- matrix(0, nphi, length(gpx))
  for (j in 1 : nphi)
    phi[j, ] <- (j^{-2}) * sqrt(2) * sin(j * pi * gpx)
  return(phi)
}

X.phi2 <- function(nphi, gpx) {
  phi <- matrix(0, nphi, length(gpx))
  for (j in 1 : nphi)
    phi[j, ] <- (j^{-2}) * sqrt(2) * cos(j * pi * gpx)
  return(phi)
}

rX.s <- function(nphi, gpx){
  xsi <- rnorm(2 * nphi, 0, 1)
  X <- xsi * rbind(X.phi1(nphi, gpx), X.phi2(nphi, gpx)) 
  Xs <- colSums(X)
  return(Xs)
}

dgp1 <- function(n, nphi = 10, gpy, gpx, sd.error, out.p = 0){
  
  alpha <- function(t) 2 * exp(-(t - 1)^2)
  beta <- function(t, s) {
    4 * cos(2 * pi * t) * sin(pi * s)
  }
  
  data <- list()
  coeffun <- list()
  s <- gpx
  t <- gpy
  ngpx <- length(gpx)
  ngpy <- length(gpy)
  alpha.t <- matrix(alpha(t), nrow = n, ncol = ngpy, byrow = TRUE)
  beta.ts <- outer(t, s, beta)
  data$X <- I(t(replicate(n, rX.s(nphi, s))))
  Xbeta.t <- (data$X) %*% t(beta.ts) * (s[2] - s[1])
  data$Xbeta <- I(Xbeta.t)
  eps <- matrix(rnorm(n * ngpy, 0, sd.error), n, ngpy)
  data$Ytrue <- I(alpha.t + Xbeta.t)
  data$Y <- I(alpha.t + Xbeta.t + eps)
  coeffun$alpha <- alpha.t
  coeffun$beta <- beta.ts
  
  if(out.p > 0){
    
    alpha.out <- function(t) 4 * exp(-(t)^2)
    beta.out <- function(t, s) {
      6 * sin(4 * pi * t) * sin(2*pi * s)
    }
    
    data2 <- list()
    
    alpha.t.out <- matrix(alpha.out(t), nrow = n, ncol = ngpy, byrow = TRUE)
    beta.ts.out <- outer(t, s, beta.out)
    
    X.out <- I(t(replicate(n, rX.s(nphi, s))))
    Xbeta.t.out <- (X.out) %*% t(beta.ts.out) * (s[2] - s[1])
    eps2 <- matrix(rnorm(n * ngpy, 0, sd.error), n, ngpy)
    Y.out <- I(alpha.t.out + Xbeta.t.out + eps2)
    
    nout <- round(n * out.p)
    out.indx <- sample(1:n, nout)
    data$Y[out.indx,] <- Y.out[out.indx,]
  }
  return(structure(as.data.frame(data, row.names = 1 : n), xindex = s, yindex = t,
                   truth = coeffun))
}

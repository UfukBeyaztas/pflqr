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

dgp3 <- function(n, nphi = 10, gpy, gpx, rho=0.8){
  
  alpha <- function(t) 2 * sin(4*pi*t)
  beta <- function(t, s) {
    4 * cos(4 * pi * t) * sin(4*pi * s)
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
  Sigma <- matrix(rho, nrow=n, ncol=n)
  diag(Sigma) <- 1
  eps <- t(mvrnorm(length(t), mu = rep(0, n), Sigma = Sigma))
  data$Ytrue <- I(alpha.t + Xbeta.t)
  data$Y <- I(alpha.t + Xbeta.t + eps)
  coeffun$alpha <- alpha.t
  coeffun$beta <- beta.ts
  
  data$X <- data$X + matrix(rnorm(n*ngpx), ncol = ngpx)
  data$Y <- data$Y + matrix(rnorm(n*ngpy), ncol = ngpy)
  
  return(structure(as.data.frame(data, row.names = 1 : n), xindex = s, yindex = t,
                   truth = coeffun))
}

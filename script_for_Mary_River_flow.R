source("auxilary_functions.R")
source("auxilary_functions_fpca.R")

gpx <- 1:24
gpy <- 1:24

library(robflreg)
data("MaryRiverFlow")
MaryRiverFlow <- as.matrix(MaryRiverFlow)

var_fun = function(data,order){
  n = dim(data)[1]
  y = data[((order+1):n),]
  x=list()
  a=1
  b=order
  for(i in 1:order){
    x[[i]] = data[(a:(n-b)),]
    a = a+1
    b = b-1
  }
  return(list(x=x,y=y))
}

starting_value <- 2139

mspe.pffr <- numeric()
mspe.fpca <- numeric()
mspe.pffqr <- numeric()

sco.cpd.fpca <- matrix(, ncol = 2, nrow = 50)
sco.cpd.pffqr <- matrix(, ncol = 2, nrow = 50)

for(sim in 1:50){
  try({
    
    data.i <- MaryRiverFlow[1:starting_value,]
    XY = var_fun(data=data.i, order=1)
    X = XY$x[[1]]
    Y = XY$y
    
    X.test <- matrix(Y[dim(Y)[1],], nrow = 1)
    pffqr.model <- fpenqr(y=Y, x=X, tau=0.5, nby=10, nbx=10, nb0=10, gpy=gpy, gpx=gpx)
    pffqr.model.025 <- fpenqr(y=Y, x=X, tau=0.025, nby=10, nbx=10, nb0=10, gpy=gpy, gpx=gpx)
    pffqr.model.975 <- fpenqr(y=Y, x=X, tau=0.975, nby=10, nbx=10, nb0=10, gpy=gpy, gpx=gpx)
    
    
    
    fpca.model <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.5, gpy=gpy, gpx=list(gpx),
                       fnbasisY=10, fnbasisX=10, fncomp_Y=5, fncomp_X=5)
    fpca.model.025 <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.025, gpy=gpy, gpx=list(gpx),
                           fnbasisY=10, fnbasisX=10, fncomp_Y=5, fncomp_X=5)
    fpca.model.975 <- fpca(fY=Y, fX=list(X), fX_test=list(X.test), tau=0.975, gpy=gpy, gpx=list(gpx),
                           fnbasisY=10, fnbasisX=10, fncomp_Y=5, fncomp_X=5)
    
    data <- data.frame(Y=Y, X=X)
    pffr.model <- pffr(Y ~ ff(X, xind = gpx), yind = gpy, data = data, 
                       bs.int = list(bs = "ps", k = 10, m = c(2, 2)),
                       bs.yindex = list(bs = "ps", k = 10, m = c(2, 2)))
    
    
    preds.pffr <- predict(pffr.model, newdata = list(X = X.test))
    preds.fpca <- fpca.model$pred
    preds.pffqr <- predict.fpenqr(object = pffqr.model, xnew = X.test)
    
    preds.pffqr.025 <- predict.fpenqr(object = pffqr.model.025, xnew = X.test)
    preds.pffqr.975 <- predict.fpenqr(object = pffqr.model.975, xnew = X.test)
    
    Y.test <- MaryRiverFlow[starting_value+1,]
    mspe.pffr[sim] <- mean((Y.test-preds.pffr)^2)
    mspe.fpca[sim] <- mean((Y.test-preds.fpca)^2)
    mspe.pffqr[sim] <- mean((Y.test-preds.pffqr)^2)

    sco.cpd.fpca[sim,] <- interval_score(Y.test,fpca.model.025$pred, fpca.model.975$pred, 0.05)
    sco.cpd.pffqr[sim,] <- interval_score(Y.test,preds.pffqr.025, preds.pffqr.975, 0.05)
    starting_value <- starting_value +1
  }, silent = TRUE)
}


mean(mspe.pffr);sd(mspe.pffr)
mean(mspe.fpca);sd(mspe.fpca)
mean(mspe.pffqr);sd(mspe.pffqr)

apply(sco.cpd.fpca,2,mean);apply(sco.cpd.fpca,2,sd)
apply(sco.cpd.pffqr,2,mean);apply(sco.cpd.pffqr,2,sd)

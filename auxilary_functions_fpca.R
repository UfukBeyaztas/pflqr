#_______________________________________________________
#_________________Estimation_function___________________
#_______________________________________________________
# This function is used to obtain parameter estimates
# sco_Y is the PC scores of the response variable
# sco_X is the PC scores of the predictors
# taus is the tau level
est_fun = function(sco_Y, sco_X, tau){
  Bhat = t(rqs.fit(x = cbind(1,sco_X), y = sco_Y,
                   tau = tau))
  return(Bhat)
}
#_______________________________________________________
#_______________________________________________________


#________________________________________
#__________Prediction_for_Y______________
#________________________________________
# This function is used to obtain fitted
# and predicted functions
# comp_Y is the PC of response variable
# sco_X is the PC scores of the predictors
# Bhat is the estimated parameter vector
pred_fun = function(comp_Y, sco_X, Bhat){
  ncomp = dim(comp_Y$coefs)[2]
  nest = t(sco_X %*% Bhat)
  if(ncomp == 1){
    nh = nest[1] * comp_Y[1,]
  }else{
    nh = nest[1] * comp_Y[1,]
    for (j in 2:ncomp){
      nh = nh + nest[j] * comp_Y[j,]  
    }
  }
  return(nh)
}
#________________________________________
#________________________________________


#_______________________________________________________________________________________________________
#____________________________________________FPCA_main_effect___________________________________________
#_______________________________________________________________________________________________________
# This function is used to obtain PC scores and components
# X_tr is a list of functional predictors in the training sample
# X_te is a list of functional predictors in the test sample
# nbasis is the number of basis functions used to approximate FPCA basis functions
# ncomp is the number of PCA components
# cen is a selection of TRUE/FALSE
getPCA = function(X_tr, X_te, nbasis, ncomp, gp, cen){
  data = rbind(X_tr, X_te)
  n_train = dim(X_tr)[1]
  n_test = dim(X_te)[1]
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = gp
  bs_basis = create.bspline.basis(c(gp[1],gp[p]), nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  dpca = pca.fd(pcaobj[1:n_train,], nharm = ncomp, fdobj, centerfns = cen)
  PCAscore = dpca$scores
  PCAcoef = dpca$harmonics
  mean_coef = dpca$meanfd
  fdXnew = pcaobj[-(1:n_train),]
  fdXnew$coefs = t(scale(t(fdXnew$coefs), scale = F))
  PCAscore_test = inprod(fdXnew, dpca$harmonics)
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, PCAscore_test = PCAscore_test,
              meanScore = mean_coef, evalbase = evalbase))
}
#_______________________________________________________________________________________________________
#_______________________________________________________________________________________________________

#___________________________________________
#______________Check_loss_function__________
#___________________________________________
clf = function(data, tau){
  data * (tau - (data < 0))
}
#___________________________________________
#___________________________________________


#______________________________________________________________________
#_________________________BIC_function_(nc)____________________________
#______________________________________________________________________
# This function is used to obtain BIC value
# Y is the functional response
# Yfit is the fitted functional response
# nc is the total number of PCs
# tau is the tau level
BIC_fun = function(Y, Yfit, nc, tau){
  n = dim(Y)[1]
  
  res = Y - Yfit
  
  BIC_val = sum(clf(res, tau)) + nc * log(n)
  
  return(BIC_val)
}
#__________________________________________________________________
#__________________________________________________________________


#___________________________________________________________________________________________________________
#_____________________________________________________BIC___________________________________________________
#___________________________________________________________________________________________________________
# This function is used to compute BIC value for each number of PC
# Y is the functional response
# X is a list of functional predictors
# ncmaxY is the maximum number of PC for the response variable
# ncmaxX is the maximum number of PC for the predictor variables
# rangeval_Y contains the interval values where the functional response is evaluated
# rangeval_X is a list containing the interval values where the functional predictors are evaluated
BIC = function(Y, X, tau, ncmaxY, ncmaxX, gpy, gpx, nbasisY, nbasisX){
  fnp = length(X)
  fn = dim(Y)[1]
  fp = dim(Y)[2]
  
  BIC_mat = matrix(, nrow = ncmaxX, ncol = ncmaxY)
  for(i in 1:ncmaxY){
    for(j in 1:ncmaxX){
      fPCA_Y = getPCA(X_tr = Y, X_te = NULL, nbasis = nbasisY, ncomp = i,
                      gp = gpy, cen = FALSE)
      fsco_Y = fPCA_Y$PCAscore
      fcomp_Y = fPCA_Y$PCAcoef
      fmean_Y = fPCA_Y$meanScore
      evaly = fPCA_Y$evalbase
      
      fsco_X = list()
      evalx = list()
      fcomp_X = list()
      for(fij in 1:fnp){
        fPCA_X = getPCA(X_tr = X[[fij]], X_te = NULL, nbasis = nbasisX[[fij]], ncomp = j,
                        gp = gpx[[fij]], cen = TRUE)
        fsco_X[[fij]] = fPCA_X$PCAscore
        evalx[[fij]] = fPCA_X$evalbase
        fcomp_X[[fij]] = fPCA_X$PCAcoef
      }
      
      fBhat = est_fun(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X), tau = tau)
      B0 = fBhat[1,]
      fBhat = fBhat[-1,]
      fB0 = B0 %*% t(fcomp_Y$coefs) %*% t(evaly)
      
      
      fYfit = matrix(, nrow = fn, ncol = fp)
      for(fk in 1:fn){
        fXk = do.call(cbind, fsco_X)[fk,]
        fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat)
        fYfit[fk,] = eval.fd(gpy, fmodel_k) + t(fB0)
      }
      
      BIC_mat[(i), (j)] = BIC_fun(Y = Y, Yfit = fYfit, nc = (i+j), tau = tau)
    }
  }
  BIC_opt = which(BIC_mat == min(BIC_mat), arr.ind = TRUE)
  ncomp_opt_X = BIC_opt[1]
  ncomp_opt_Y = BIC_opt[2]
  return(list(ncompX = ncomp_opt_X, ncompY = ncomp_opt_Y, bic = BIC_mat[ncomp_opt_X, ncomp_opt_Y]))
}
#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________


#______________________________________________________________________
#_________________________BIC_function_(M)_____________________________
#______________________________________________________________________
# This function is used to compute BIC value in the 
# variable selection procedure
BICm_fun = function(Y, Yfit, nv, tau){
  n = dim(Y)[1]
  
  res = Y - Yfit
  
  BIC_val = sum(clf(res, tau)) + nv * log(n) / (2*n)
  
  return(BIC_val)
}
#______________________________________________________________________
#______________________________________________________________________



#_______________________________________________________________________________________________________________
#___________________________________________________FPCA_model__________________________________________________
#_______________________________________________________________________________________________________________

fpca = function(fY, fX, fX_test, tau, gpy, gpx,
                fnbasisY, fnbasisX, fncomp_Y, fncomp_X){
  
    BIC_res = BIC(Y=fY, X=fX, tau=tau, ncmaxY=fncomp_Y, ncmaxX=fncomp_X, gpy=gpy, gpx=gpx,
                  nbasisY=fnbasisY, nbasisX=fnbasisX)
    
    fnp = length(fX)
    fn = dim(fY)[1]
    fn_test = dim(fX_test[[1]])[1]
    fp = dim(fY)[2]
    
    ncY = BIC_res$ncompY
    ncX = BIC_res$ncompX
    
    fPCA_Y = getPCA(X_tr = fY, X_te = NULL, nbasis = fnbasisY, ncomp = ncY,
                    gp = gpy, cen = FALSE)
    fsco_Y = fPCA_Y$PCAscore
    fcomp_Y = fPCA_Y$PCAcoef
    fmean_Y = fPCA_Y$meanScore
    evaly = fPCA_Y$evalbase
    
    fsco_X = list()
    fsco_X_test = list()
    evalx = list()
    fcomp_X = list()
    for(fij in 1:fnp){
      fPCA_X = getPCA(X_tr = fX[[fij]], X_te = fX_test[[fij]], nbasis = fnbasisX[[fij]], ncomp = ncX,
                      gp = gpx[[fij]], cen = TRUE)
      fsco_X[[fij]] = fPCA_X$PCAscore
      fsco_X_test[[fij]] = fPCA_X$PCAscore_test
      evalx[[fij]] = fPCA_X$evalbase
      fcomp_X[[fij]] = fPCA_X$PCAcoef
    }
    
    fBhat = est_fun(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X), tau = tau)
    B0 = fBhat[1,]
    fc = dim(fBhat)[2]
    fBhat = fBhat[-1,]
    if(is.null(dim(fBhat)) == TRUE & fc >1){
      fBhat = t(as.matrix(fBhat))
    }else if(is.null(dim(fBhat)) == TRUE & fc == 1){
      fBhat = as.matrix(fBhat)
    }
    fB0 = B0 %*% t(fcomp_Y$coefs) %*% t(evaly)
    
    
    fYfit = matrix(, nrow = fn, ncol = fp)
    for(fk in 1:fn){
      fXk = do.call(cbind, fsco_X)[fk,]
      fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat)
      fYfit[fk,] = eval.fd(gpy, fmodel_k) + t(fB0)
    }
    
    fYpred = matrix(, nrow = fn_test, ncol = fp)
    for(fk in 1:fn_test){
      fXk = do.call(cbind, fsco_X_test)[fk,]
      fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat)
      fYpred[fk,] = eval.fd(gpy, fmodel_k) + t(fB0)
    }
    
    coef_main = list()
    km = 1
    for(im in 1:fnp){
      coef_main[[im]] = t(evalx[[im]] %*% (fcomp_X[[im]]$coefs %*% 
                                             fBhat[km: (km+ncX-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly))
      km = im*ncX+1
    }
    
  return(list(fit = fYfit,
              pred = fYpred,
              coef = coef_main,
              intercep = fB0))
}
#_______________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________



interval_score <- function(holdout, lb, ub, alpha){
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
  cpd = abs(cover - (1 - alpha))
  return(c(mean(score), cpd))
}

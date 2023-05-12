## load packages
library(knitr)

## load help functions
source("function.R")

#########################################
########## Relative Efficiency ##########
#########################################

# stage 1 working correlation structures considered
corstr.v <- c("unstructured", "exchangeable", "independence")

# relative efficiency based on variance
calctheta <- function(output){
  # estimated stage 1 heterogeneous treatment effect conditional on Z1 = 0
  dat1.psd <- cbind.data.frame("Z1" = rep(0, 3), "A1" = rep(1, 3), "time" = c(1, 2, 3)) # pseudo data
  dat0.psd <- cbind.data.frame("Z1" = rep(0, 3), "A1" = rep(-1, 3), "time" = c(1, 2, 3))
  delta.hat <- predict(output$mod.s1, newdata = dat1.psd) - predict(output$mod.s1, newdata = dat0.psd)
  return(delta.hat)
}

simVarTheta <- function(I = 1000, alpha1 = 0.7, alpha2 = 1, alpha3 = 1.2, beta2 = 0.5, beta3 = 0.8, gamma2 = 1, gamma3 = 1.5, 
                       sigma.v = 4, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1){
  # matrix to store estimates
  var.mat <- cbind.data.frame("corstr" = corstr.v, matrix(0, nrow = 3, ncol = 3*2))
  colnames(var.mat)[2:7] <- apply(expand.grid(c("Time1", "Time2", "Time3"), c("SQGEE", "MQGEE")), 1, paste, collapse = ".")
  for (j in 1:length(corstr.v)){
    est.SQGEE <- matrix(NA, nrow = I, ncol = 3)
    est.MQGEE <- matrix(NA, nrow = I, ncol = 3)
    for (i in 1:I){
      dat <- generateData.joint(N = 200, 
                                alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                                sigma.v = sigma.v, sigma.e = sigma.e, 
                                lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
      est.SQGEE[i, ] <- calctheta(QGEE(dat, method = "SQGEE", corstr = c(corstr.v[j], "exchangeable"), w = c(1/3, 1/3, 1/3), misspec = FALSE))
    }
    for (i in 1:I){
      dat <- generateData.joint(N = 200, 
                                alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                                sigma.v = sigma.v, sigma.e = sigma.e, 
                                lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
      est.MQGEE[i, ] <- calctheta(QGEE(dat, method = "MQGEE", corstr = c(corstr.v[j], "exchangeable"), w = c(1/3, 1/3, 1/3), misspec = FALSE))
    }
    var.mat[j, 2:7] <- c(apply(est.SQGEE, 2, var), apply(est.MQGEE, 2, var))
  }
  return(var.mat)
}

var.low.pos <- simVarTheta(I = 1000, 
                           alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                           sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = 0.2, lambda3 = 0.2)
var.low.neg23 <- simVarTheta(I = 1000, 
                             alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                             sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = -0.2, lambda3 = -0.2)

var.med.pos <- simVarTheta(I = 1000, 
                           alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                           sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = 1, lambda3 = 1)
var.med.ind <- simVarTheta(I = 1000, 
                           alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                           sigma.v = 3.2, sigma.e = 2.7, lambda1 = 0, lambda2 = 0, lambda3 = 0)
var.med.neg23 <- simVarTheta(I = 1000, 
                             alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                             sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = -1, lambda3 = -1)

var.hgh.pos <- simVarTheta(I = 1000, 
                           alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                           sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = 3, lambda3 = 1)
var.hgh.neg23 <- simVarTheta(I = 1000, 
                             alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                             sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = -3, lambda3 = -1)

calcREVar <- function(res.list){
  res.dat <- cbind.data.frame("corstr" = rep(c("unstructured/independence", "exchangeable/independence"), length(res.list)), 
                              "correlation" = rep(c("med.pos", "med.ind", "med.neg", "low.pos", "low.neg", "high.pos", "high.neg"), each = 2), 
                              matrix(NA, nrow = 2*length(res.list), ncol = 6))
  colnames(res.dat)[3:8] <- apply(expand.grid(c("Time1", "Time2", "Time3"), c("SQGEE", "MQGEE")), 1, paste, collapse = ".")
  for (i in 1:length(res.list)){
    res <- res.list[[i]]
    res.dat[2*i-1, 3:8] <- round(res[1, 2:7]/res[3, 2:7], 2)
    res.dat[2*i, 3:8] <- round(res[2, 2:7]/res[3, 2:7], 2)
  }
  return(res.dat)
}
res.list <- list(var.med.pos, var.med.ind, var.med.neg23, var.low.pos, var.low.neg23, var.hgh.pos, var.hgh.neg23)
kable(calcREVar(res.list), "latex", longtable = T, booktabs = T, caption = "SimRE")

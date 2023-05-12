### code for simulation studies

## load help functions
source("function.R")

## two ways to calculate the true optimal marginal trajectory
# 1. calculate true optimal rules directly
calcOptMargTraj <- function(dat, w = c(1/3, 1/3, 1/3), 
                            alpha1, alpha2, alpha3, gamma2, gamma3, 
                            sigma.v, sigma.e, lambda1, lambda2, lambda3){
  # optimal rule
  # d2opt <- ifelse((w[2]*gamma2+w[3]*gamma3) * dat$A1 > 0, -1, 1)
  # frac2 <- gamma2 + lambda1*lambda2*sigma.v^2/(lambda1*sigma.v^2+sigma.e^2)*(dat$Y1-alpha1*dat$Z1*dat$A1)
  # frac3 <- gamma3 + lambda1*lambda3*sigma.v^2/(lambda1*sigma.v^2+sigma.e^2)*(dat$Y1-alpha1*dat$Z1*dat$A1)
  # temp <- (w[2]*frac2+w[3]*frac3) * dat$A1
  # d2opt <- ifelse(temp > 0, -1, 1)
  d1 <- function(a1){
    v <- ((w[1]*alpha1+w[2]*alpha2+w[3]*alpha3) * dat$Z1 + 
            (gamma2+gamma3) * ifelse((w[2]*gamma2+w[3]*gamma3) * a1 > 0, -1, 1)) * a1
    return(v)
  }
  d1opt <- ifelse(d1(-1) < d1(1), -1, 1)
  d2opt <- ifelse((w[2]*gamma2+w[3]*gamma3) * d1opt > 0, -1, 1)
  # trajectory under the optimal rule
  Y1 <- alpha1*dat$Z1*d1opt
  Y2 <- alpha2*dat$Z1*d1opt + gamma2*d1opt*d2opt
  Y3 <- alpha3*dat$Z1*d1opt + gamma3*d1opt*d2opt
  # Y2 <- alpha2*dat$Z1*d1opt + gamma2*d1opt*d2opt + lambda1*lambda2*sigma.v^2/(lambda1*sigma.v^2+sigma.e^2)*(dat$Y1-alpha1*dat$Z1*d1opt)
  # Y3 <- alpha3*dat$Z1*d1opt + gamma3*d1opt*d2opt + lambda1*lambda3*sigma.v^2/(lambda1*sigma.v^2+sigma.e^2)*(dat$Y1-alpha1*dat$Z1*d1opt)
  traj.true <- c(mean(Y1), mean(Y2), mean(Y3))
  return(traj.true)
}

# 2. calculate the individual trajectory for all DTRs
calcIndValue <- function(dat, d1, d2, w = c(1/3, 1/3, 1/3), 
                         alpha1, alpha2, alpha3, gamma2, gamma3, 
                         sigma.v, sigma.e, lambda1, lambda2, lambda3){
  # trajectory under the optimal rule
  Y1 <- alpha1*dat$Z1*d1
  Y2 <- alpha2*dat$Z1*d1 + gamma2*d1*d2
  Y3 <- alpha3*dat$Z1*d1 + gamma3*d1*d2
  value.ind <- w[1]*Y1 + w[2]*Y2 + w[3]*Y3
  return(value.ind)
}

calcOptMargTraj <- function(dat, w = c(1/3, 1/3, 1/3), 
                            alpha1, alpha2, alpha3, gamma2, gamma3, 
                            sigma.v, sigma.e, lambda1, lambda2, lambda3){
  temp <- cbind.data.frame("Y.1.1" = calcIndValue(dat = dat, d1 = 1, d2 = 1, w = w, 
                                                  alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, gamma2 = gamma2, gamma3 = gamma3, 
                                                  sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3), 
                           "Y.1.n1" = calcIndValue(dat = dat, d1 = 1, d2 = -1, w = w, 
                                                   alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, gamma2 = gamma2, gamma3 = gamma3, 
                                                   sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3), 
                           "Y.n1.1" = calcIndValue(dat = dat, d1 = -1, d2 = 1, w = w, 
                                                   alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, gamma2 = gamma2, gamma3 = gamma3, 
                                                   sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3), 
                           "Y.n1.n1" = calcIndValue(dat = dat, d1 = -1, d2 = -1, w = w, 
                                                    alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, gamma2 = gamma2, gamma3 = gamma3, 
                                                    sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3))
  index <- apply(temp, 1, which.min)
  d1 <- ifelse(index == 1 | index == 2, 1, -1)
  d2 <- ifelse(index == 1 | index == 3, 1, -1)
  Y1 <- alpha1*dat$Z1*d1
  Y2 <- alpha2*dat$Z1*d1 + gamma2*d1*d2
  Y3 <- alpha3*dat$Z1*d1 + gamma3*d1*d2
  return(c(mean(Y1), mean(Y2), mean(Y3)))
}

# set sample size to 200
simulation <- function(I = 1000, w = c(1/3, 1/3, 1/3), 
                       alpha1 = 0.7, alpha2 = 1, alpha3 = 1.2, gamma2 = 1, gamma3 = 1.5, 
                       sigma.v = 4, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                       misspec = FALSE){
  bias.sqgee <- matrix(NA, nrow = I, ncol = 3)
  bias.mqgee <- matrix(NA, nrow = I, ncol = 3)
  bias.sq <- matrix(NA, nrow = I, ncol = 3)
  bias.mq <- matrix(NA, nrow = I, ncol = 3)

  for (i in 1:I){
    dat <- generateData.joint(N = 200, 
                              alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, gamma2 = gamma2, gamma3 = gamma3, 
                              sigma.v = sigma.v, sigma.e = sigma.e, 
                              lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
    sqgee <- QGEE(dat, method = "SQGEE", corstr = c("unstructured", "exchangeable"), w = w, misspec = misspec)
    mqgee <- QGEE(dat, method = "MQGEE", corstr = c("unstructured", "exchangeable"), w = w, misspec = misspec)
    sq <- Q(dat, method = "SQ", w = w, misspec = misspec)
    mq <- Q(dat, method = "MQ", w = w, misspec = misspec)
    traj.true <- calcOptMargTraj(dat = dat, w = w, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, gamma2 = gamma2, gamma3 = gamma3, 
                                 sigma.v, sigma.e, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
    bias.sqgee[i, ] <- sqgee$traj - traj.true
    bias.mqgee[i, ] <- mqgee$traj - traj.true
    bias.sq[i, ] <- sq$traj - traj.true
    bias.mq[i, ] <- mq$traj - traj.true
  }
  
  mySum <- function(x){
    paste(round(mean(x, na.rm = TRUE), 3), " (", round(sd(x, na.rm = TRUE), 3), ")", sep = "")
  }
  sum <- c(apply(bias.sqgee, 2, mySum), apply(bias.mqgee, 2, mySum), 
           apply(bias.sq, 2, mySum), apply(bias.mq, 2, mySum))
  bias <- cbind.data.frame("method" = rep(c("sqgee", "mqgee", "sq", "mq")), 
                           "bias1" = sum[c(1,4,7,10)],
                           "bias2" = sum[c(2,5,8,11)], 
                           "bias3" = sum[c(3,6,9,12)])
  return(bias)
}

set.seed(1001)
N.sim = 1000
sim.pos <- simulation(I = N.sim, 
                      alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, gamma2 = -1.2, gamma3 = -0.8, 
                      sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                      misspec = FALSE)
sim.ind <- simulation(I = N.sim, 
                      alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, gamma2 = -1.2, gamma3 = -0.8, 
                      sigma.v = 3.2, sigma.e = 2.7, lambda1 = 0, lambda2 = 0, lambda3 = 0, 
                      misspec = FALSE)
sim.neg23 <- simulation(I = N.sim, 
                        alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, gamma2 = -1.2, gamma3 = -0.8, 
                        sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = -1, lambda3 = -1, 
                        misspec = FALSE)
sim.pos.mis <- simulation(I = N.sim, 
                          alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, gamma2 = -1.2, gamma3 = -0.8, 
                          sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                          misspec = TRUE)
sim.ind.mis <- simulation(I = N.sim, 
                          alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, gamma2 = -1.2, gamma3 = -0.8, 
                          sigma.v = 3.2, sigma.e = 2.7, lambda1 = 0, lambda2 = 0, lambda3 = 0, 
                          misspec = TRUE)
sim.neg23.mis <- simulation(I = N.sim, 
                            alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, gamma2 = -1.2, gamma3 = -0.8, 
                            sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = -1, lambda3 = -1, 
                            misspec = TRUE)

### simulation study code

## load packages

## define metrics that assess the performance of a method

metric <- function(output, alpha1, alpha2, alpha3, beta2, beta3, gamma2, gamma3){
  # probability of correctly identifying the stage-specific optimal rule (PCI)
  d2opt <- ifelse((beta2+beta3) * output$dat$Y0 + (gamma2+gamma3) * output$dat$A1 < 0, -1, 1)
  pci.s2 <- sum(output$dat.opt$A2 == d2opt) / nrow(output$dat.opt)
  d1 <- function(a1){
    v <- ((alpha1+alpha2+alpha3) *  output$dat$Y0 + (gamma2+gamma3) * ifelse((beta2+beta3) * output$dat$Y0 + (gamma2+gamma3) * a1 < 0, -1, 1)) * a1
    return(v)
  }
  d1opt <- ifelse(d1(-1) > d1(1), -1, 1)
  pci.s1 <- sum(output$dat.opt$A1 == d1opt) / nrow(output$dat.opt)
  
  # root mean square error (RMSE) of heterogeneous causal effect estimators
  # true stage 1 heterogeneous treatment effect at time 1 (stage 1), 2 and 3 (stage 2)
  # temp <- function(a1){
  #   return(ifelse((beta2+beta3)*output$dat$Y0 + (gamma2+gamma3)*a1 < 0, -1, 1))
  # }
  delta1 <- 2*alpha1*output$dat.opt$Y0
  # delta2 <- 2*alpha2*output$dat.opt$Y0 + gamma2*(temp(1)+temp(-1))
  # delta3 <- 2*alpha3*output$dat.opt$Y0 + gamma3*(temp(1)+temp(-1))
  delta2 <- 2*alpha2*output$dat.opt$Y0 + 2*gamma2*output$dat.opt$A2
  delta3 <- 2*alpha3*output$dat.opt$Y0 + 2*gamma3*output$dat.opt$A2
  # estimated stage 1 heterogeneous treatment effect
  dat <- cbind.data.frame(do.call(rbind, replicate(3, output$dat.opt[ , c("ID", "Y0", "A1")], simplify = FALSE)), 
                          "time" = rep(c(1, 2, 3), each = nrow(output$dat.opt)))
  dat <- dat[order(dat$ID), ]
  dat1 <- dat
  dat1$A1 <- 1
  dat0 <- dat
  dat0$A1 <- -1
  delta.hat <- matrix(predict(output$mod.s1, newdata = dat1), ncol = 3, byrow = TRUE) - matrix(predict(output$mod.s1, newdata = dat0), ncol = 3, byrow = TRUE)
  # rmse
  rmse1 <- sqrt(mean((delta.hat[, 1] - delta1)^2))
  rmse2 <- sqrt(mean((delta.hat[, 2] - delta2)^2))
  rmse3 <- sqrt(mean((delta.hat[, 3] - delta3)^2))
  rmse <- c(rmse1, rmse2, rmse3)
  Y0.grid <- c(-3, -2, -1, 0, 1, 2, 3)
  
  # bias
  bias.mat <- matrix(NA, nrow = 3, ncol = length(Y0.grid)) # matrix to store bias
  # estimated stage 1 heterogeneous treatment effect
  delta1 <- 2*alpha1*Y0.grid
  delta2 <- 2*alpha2*Y0.grid + gamma2*(ifelse((beta2+beta3)*Y0.grid+gamma2+gamma3 < 0, -1, 1) + ifelse((beta2+beta3)*Y0.grid-gamma2-gamma3 < 0, -1, 1)) # calculate A2opt analytically
  delta3 <- 2*alpha3*Y0.grid + gamma3*(ifelse((beta2+beta3)*Y0.grid+gamma2+gamma3 < 0, -1, 1) + ifelse((beta2+beta3)*Y0.grid-gamma2-gamma3 < 0, -1, 1))
  dat1.psd <- cbind.data.frame("Y0" = rep(Y0.grid, each = 3), "A1" = rep(1, length(Y0.grid)*3), "time" = rep(c(1, 2, 3), length(Y0.grid))) # pseudo data
  dat0.psd <- cbind.data.frame("Y0" = rep(Y0.grid, each = 3), "A1" = rep(-1, length(Y0.grid)*3), "time" = rep(c(1, 2, 3), length(Y0.grid)))
  delta.hat <- matrix(predict(output$mod.s1, newdata = dat1.psd), ncol = 3, byrow = TRUE) - matrix(predict(output$mod.s1, newdata = dat0.psd), ncol = 3, byrow = TRUE)
  bias.mat[1, ] <- delta.hat[, 1] - delta1
  bias.mat[2, ] <- delta.hat[, 2] - delta2
  bias.mat[3, ] <- delta.hat[, 3] - delta3
  return(list("pci.s2" = pci.s2, "pci.s1" = pci.s1, "rmse" = rmse, "bias.mat" = bias.mat))
}

## Simulation
simulation <- function(I = 1000, 
                       alpha1 = 0.7, alpha2 = 1, alpha3 = 1.2, beta2 = 0.5, beta3 = 0.8, gamma2 = 1, gamma3 = 1.5, 
                       sigma.v = 4, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1){
  ss <- seq(40, 400, by = 40) # set sample size
  pci.sq.s1 <- matrix(NA, nrow = length(ss), ncol = I)
  pci.sq.s2 <- matrix(NA, nrow = length(ss), ncol = I)
  pci.mq.s1 <- matrix(NA, nrow = length(ss), ncol = I)
  pci.mq.s2 <- matrix(NA, nrow = length(ss), ncol = I)
  pci.sqgee.s1 <- matrix(NA, nrow = length(ss), ncol = I)
  pci.sqgee.s2 <- matrix(NA, nrow = length(ss), ncol = I)
  pci.mqgee.s1 <- matrix(NA, nrow = length(ss), ncol = I)
  pci.mqgee.s2 <- matrix(NA, nrow = length(ss), ncol = I)
  rmse1.sqgee <- matrix(NA, nrow = length(ss), ncol = I)
  rmse1.mqgee <- matrix(NA, nrow = length(ss), ncol = I)
  rmse2.sqgee <- matrix(NA, nrow = length(ss), ncol = I)
  rmse2.mqgee <- matrix(NA, nrow = length(ss), ncol = I)
  rmse3.sqgee <- matrix(NA, nrow = length(ss), ncol = I)
  rmse3.mqgee <- matrix(NA, nrow = length(ss), ncol = I)
  Y0.grid <- c(-3, -2, -1, 0, 1, 2, 3)
  bias.mat <- matrix(NA, nrow = length(Y0.grid), ncol = I)
  bias1.sqgee <- rep(list(bias.mat), length(ss))
  bias1.mqgee <- rep(list(bias.mat), length(ss))
  bias2.sqgee <- rep(list(bias.mat), length(ss))
  bias2.mqgee <- rep(list(bias.mat), length(ss))
  bias3.sqgee <- rep(list(bias.mat), length(ss))
  bias3.mqgee <- rep(list(bias.mat), length(ss))
  for (s in 1:length(ss)){
    for (i in 1:I){
      dat <- generateData.joint(N = ss[s], 
                                alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                                sigma.v = sigma.v, sigma.e = sigma.e, 
                                lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
      sq <- Q(dat, method = "SQ")
      mq <- Q(dat, method = "MQ")
      sqgee <- QGEE(dat, method = "SQGEE")
      mqgee <- QGEE(dat, method = "MQGEE")
      metrics.sq <- metric(sq, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.mq <- metric(mq, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.sqgee <- metric(sqgee, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.mqgee <- metric(mqgee, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      pci.sq.s1[s, i] <- metrics.sq$pci.s1
      pci.sq.s2[s, i] <- metrics.sq$pci.s2
      pci.mq.s1[s, i] <- metrics.mq$pci.s1
      pci.mq.s2[s, i] <- metrics.mq$pci.s2
      pci.sqgee.s1[s, i] <- metrics.sqgee$pci.s1
      pci.sqgee.s2[s, i] <- metrics.sqgee$pci.s2
      pci.mqgee.s1[s, i] <- metrics.mqgee$pci.s1
      pci.mqgee.s2[s, i] <- metrics.mqgee$pci.s2
      rmse1.sqgee[s, i] <- metrics.sqgee$rmse[1]
      rmse1.mqgee[s, i] <- metrics.mqgee$rmse[1]
      rmse2.sqgee[s, i] <- metrics.sqgee$rmse[2]
      rmse2.mqgee[s, i] <- metrics.mqgee$rmse[2]
      rmse3.sqgee[s, i] <- metrics.sqgee$rmse[3]
      rmse3.mqgee[s, i] <- metrics.mqgee$rmse[3]
      bias1.sqgee[[s]][, i] <- metrics.sqgee$bias.mat[1, ]
      bias1.mqgee[[s]][, i] <- metrics.mqgee$bias.mat[1, ]
      bias2.sqgee[[s]][, i] <- metrics.sqgee$bias.mat[2, ]
      bias2.mqgee[[s]][, i] <- metrics.mqgee$bias.mat[2, ]
      bias3.sqgee[[s]][, i] <- metrics.sqgee$bias.mat[3, ]
      bias3.mqgee[[s]][, i] <- metrics.mqgee$bias.mat[3, ]
    }
  }
  pci <- cbind.data.frame("sample.size" = rep(ss, 8), 
                          "stage" = rep(c(rep(1, length(ss)), rep(2, length(ss))), 4), 
                          "method" = rep(c("sq", "mq", "sqgee", "mqgee"), each = 2*length(ss)), 
                          "pci" = c(apply(pci.sq.s1, 1, mean), apply(pci.sq.s2, 1, mean), 
                                    apply(pci.mq.s1, 1, mean), apply(pci.mq.s2, 1, mean),
                                    apply(pci.sqgee.s1, 1, mean), apply(pci.sqgee.s2, 1, mean), 
                                    apply(pci.mqgee.s1, 1, mean), apply(pci.mqgee.s2, 1, mean)))
  rmse1.200 <- cbind.data.frame("method" = rep(c("sqgee", "mqgee"), each = I), 
                                "rmse" = c(rmse1.sqgee[5, ], rmse1.mqgee[5, ]))
  rmse2.200 <- cbind.data.frame("method" = rep(c("sqgee", "mqgee"), each = I), 
                                "rmse" = c(rmse2.sqgee[5, ], rmse2.mqgee[5, ]))
  rmse3.200 <- cbind.data.frame("method" = rep(c("sqgee", "mqgee"), each = I), 
                                "rmse" = c(rmse3.sqgee[5, ], rmse3.mqgee[5, ]))
  rmse1 <- cbind.data.frame("sample.size" = rep(ss, 2),
                            "method" = rep(c("sqgee", "mqgee"), each = length(ss)), 
                            "rmse.mean" = c(apply(rmse1.sqgee, 1, mean), apply(rmse1.mqgee, 1, mean)),
                            "rmse.sd" = c(apply(rmse1.sqgee, 1, sd), apply(rmse1.mqgee, 1, sd)),
                            "rmse.lower" = c(apply(rmse1.sqgee, 1, quantile, probs = 0.025), apply(rmse1.mqgee, 1, quantile, probs = 0.025)),
                            "rmse.upper" = c(apply(rmse1.sqgee, 1, quantile, probs = 0.975), apply(rmse1.mqgee, 1, quantile, probs = 0.975)))
  rmse2 <- cbind.data.frame("sample.size" = rep(ss, 2),
                            "method" = rep(c("sqgee", "mqgee"), each = length(ss)), 
                            "rmse.mean" = c(apply(rmse2.sqgee, 1, mean), apply(rmse2.mqgee, 1, mean)),
                            "rmse.sd" = c(apply(rmse2.sqgee, 1, sd), apply(rmse2.mqgee, 1, sd)),
                            "rmse.lower" = c(apply(rmse2.sqgee, 1, quantile, probs = 0.025), apply(rmse2.mqgee, 1, quantile, probs = 0.025)),
                            "rmse.upper" = c(apply(rmse2.sqgee, 1, quantile, probs = 0.975), apply(rmse2.mqgee, 1, quantile, probs = 0.975)))
  rmse3 <- cbind.data.frame("sample.size" = rep(ss, 2),
                            "method" = rep(c("sqgee", "mqgee"), each = length(ss)), 
                            "rmse.mean" = c(apply(rmse3.sqgee, 1, mean), apply(rmse3.mqgee, 1, mean)),
                            "rmse.sd" = c(apply(rmse3.sqgee, 1, sd), apply(rmse3.mqgee, 1, sd)),
                            "rmse.lower" = c(apply(rmse3.sqgee, 1, quantile, probs = 0.025), apply(rmse3.mqgee, 1, quantile, probs = 0.025)),
                            "rmse.upper" = c(apply(rmse3.sqgee, 1, quantile, probs = 0.975), apply(rmse3.mqgee, 1, quantile, probs = 0.975)))
  bias1.sqgee <- cbind.data.frame("Y0" = rep(Y0.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Y0.grid)), 
                                  "bias.mean" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias1.mqgee <- cbind.data.frame("Y0" = rep(Y0.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Y0.grid)), 
                                  "bias.mean" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias2.sqgee <- cbind.data.frame("Y0" = rep(Y0.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Y0.grid)), 
                                  "bias.mean" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias2.mqgee <- cbind.data.frame("Y0" = rep(Y0.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Y0.grid)), 
                                  "bias.mean" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias3.sqgee <- cbind.data.frame("Y0" = rep(Y0.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Y0.grid)), 
                                  "bias.mean" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias3.mqgee <- cbind.data.frame("Y0" = rep(Y0.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Y0.grid)), 
                                  "bias.mean" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  
  return(list(pci = pci, 
              rmse1.200 = rmse1.200, rmse2.200 = rmse2.200, rmse3.200 = rmse3.200, 
              rmse1 = rmse1, rmse2 = rmse2, rmse3 = rmse3,
              bias1.sqgee = bias1.sqgee, bias1.mqgee = bias1.mqgee,
              bias2.sqgee = bias2.sqgee, bias2.mqgee = bias2.mqgee,
              bias3.sqgee = bias3.sqgee, bias3.mqgee = bias3.mqgee))
}

set.seed(1)
# correlated Outcomes
sim.pos <- simulation(I = 1000, 
                      alpha1 = 2.0, alpha2 = 0, alpha3 = 0, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                      sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1)
sim.ind <- simulation(I = 1000, 
                      alpha1 = 2.0, alpha2 = 0, alpha3 = 0, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                      sigma.v = 5, sigma.e = 3, lambda1 = 0, lambda2 = 0, lambda3 = 0)
sim.neg23 <- simulation(I = 1000, 
                        alpha1 = 2.0, alpha2 = 0, alpha3 = 0, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                        sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = -1, lambda3 = -1)

# model Misspecification
sim.pos.mis <- simulation(I = 1000, 
                          alpha1 = 2.0, alpha2 = 1.7, alpha3 = 1.6, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                          sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1)
sim.ind.mis <- simulation(I = 1000, 
                          alpha1 = 2.0, alpha2 = 1.7, alpha3 = 1.6, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                          sigma.v = 5, sigma.e = 3, lambda1 = 0, lambda2 = 0, lambda3 = 0)
sim.neg23.mis <- simulation(I = 1000, 
                            alpha1 = 2.0, alpha2 = 1.7, alpha3 = 1.6, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                            sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = -1, lambda3 = -1)

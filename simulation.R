### code for simulation studies

## load packages
library(knitr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(gtable)
library(lemon)
library(scales)

## define metrics to assess the performance of a method: PCI, RMSE, bias
# PCI is applicable to output of both functions "Q" and "QGEE", whereas RMSE and bias are only applicable to output of "QGEE"

metric <- function(output, alpha1, alpha2, alpha3, beta2, beta3, gamma2, gamma3){
  # probability of correctly identifying (PCI) the stage-specific optimal rule
  d2opt <- ifelse((beta2+beta3) * output$dat$Y0 + (gamma2+gamma3) * output$dat$A1 < 0, -1, 1) # true optimal rule at stage 2, as a function of Y0 and A1
  pci.s2 <- sum(output$dat.opt$A2 == d2opt) / nrow(output$dat.opt) # stage 2 PCI
  d1 <- function(a1){ # function to calculate the true value of average repeated-measures outcomes given baseline covariate Y0 and treatment A1 = a1, A2 = d2opt
    v <- ((alpha1+alpha2+alpha3) *  output$dat$Y0 + (gamma2+gamma3) * ifelse((beta2+beta3) * output$dat$Y0 + (gamma2+gamma3) * a1 < 0, -1, 1)) * a1
    return(v)
  }
  d1opt <- ifelse(d1(-1) > d1(1), -1, 1) # true optimal rule at stage 1, assuming that optimal rule is followed at stage 2
  pci.s1 <- sum(output$dat.opt$A1 == d1opt) / nrow(output$dat.opt) # stage 1 PCI
  
  # root mean square error (RMSE) of heterogeneous causal effect estimators
  delta1 <- 2*alpha1*output$dat.opt$Y0 # true stage 1 individual treatment effects at time 1
  delta2 <- 2*alpha2*output$dat.opt$Y0 + 2*gamma2*output$dat.opt$A2 # true stage 1 individual treatment effects at time 2, evaluated at ESTIMATED stage 2 optimal rule
  delta3 <- 2*alpha3*output$dat.opt$Y0 + 2*gamma3*output$dat.opt$A2 # true stage 1 individual treatment effects at time 3, evaluated at ESTIMATED stage 2 optimal rule
  # estimated stage 1 individual treatment effects
  dat <- cbind.data.frame(do.call(rbind, replicate(3, output$dat.opt[ , c("ID", "Y0", "A1")], simplify = FALSE)), 
                          "time" = rep(c(1, 2, 3), each = nrow(output$dat.opt)))
  dat <- dat[order(dat$ID), ]
  dat1 <- dat
  dat1$A1 <- 1
  dat0 <- dat
  dat0$A1 <- -1
  delta.hat <- matrix(predict(output$mod.s1, newdata = dat1), ncol = 3, byrow = TRUE) - matrix(predict(output$mod.s1, newdata = dat0), ncol = 3, byrow = TRUE)
  # calculate RMSE
  rmse1 <- sqrt(mean((delta.hat[, 1] - delta1)^2))
  rmse2 <- sqrt(mean((delta.hat[, 2] - delta2)^2))
  rmse3 <- sqrt(mean((delta.hat[, 3] - delta3)^2))
  rmse <- c(rmse1, rmse2, rmse3) # RMSE of heterogeneous stage 1 causal effects at times 1, 2, 3, respectively
  
  # bias of heterogeneous causal effect estimators
  Y0.grid <- c(-3, -2, -1, 0, 1, 2, 3) # a grid of Y0 values
  bias.mat <- matrix(NA, nrow = 3, ncol = length(Y0.grid)) # matrix to store bias
  # true stage 1 individual treatment effects at times 1, 2, 3, respectively, based on TRUE stage 2 optimal rule
  delta1 <- 2*alpha1*Y0.grid
  delta2 <- 2*alpha2*Y0.grid + gamma2*(ifelse((beta2+beta3)*Y0.grid+gamma2+gamma3 < 0, -1, 1) + ifelse((beta2+beta3)*Y0.grid-gamma2-gamma3 < 0, -1, 1)) # calculate A2opt analytically
  delta3 <- 2*alpha3*Y0.grid + gamma3*(ifelse((beta2+beta3)*Y0.grid+gamma2+gamma3 < 0, -1, 1) + ifelse((beta2+beta3)*Y0.grid-gamma2-gamma3 < 0, -1, 1))
  # estimated stage 1 individual treatment effects
  dat1.psd <- cbind.data.frame("Y0" = rep(Y0.grid, each = 3), "A1" = rep(1, length(Y0.grid)*3), "time" = rep(c(1, 2, 3), length(Y0.grid))) # pseudo data with A1 = 1
  dat0.psd <- cbind.data.frame("Y0" = rep(Y0.grid, each = 3), "A1" = rep(-1, length(Y0.grid)*3), "time" = rep(c(1, 2, 3), length(Y0.grid))) # pseudo data with A1 = -1
  delta.hat <- matrix(predict(output$mod.s1, newdata = dat1.psd), ncol = 3, byrow = TRUE) - matrix(predict(output$mod.s1, newdata = dat0.psd), ncol = 3, byrow = TRUE)
  # store bias values
  bias.mat[1, ] <- delta.hat[, 1] - delta1
  bias.mat[2, ] <- delta.hat[, 2] - delta2
  bias.mat[3, ] <- delta.hat[, 3] - delta3
  
  # return the calculated metrics
  return(list("pci.s2" = pci.s2, "pci.s1" = pci.s1, "rmse" = rmse, "bias.mat" = bias.mat))
}

## a function to carry out the simulation study

simulation <- function(I = 1000, # number of simulations
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

## simulation studies: based on 1000 simulations

set.seed(1)

# correlated outcomes

sim.pos <- simulation(I = 1000, 
                      alpha1 = 2.0, alpha2 = 0, alpha3 = 0, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                      sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1)
sim.ind <- simulation(I = 1000, 
                      alpha1 = 2.0, alpha2 = 0, alpha3 = 0, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                      sigma.v = 5, sigma.e = 3, lambda1 = 0, lambda2 = 0, lambda3 = 0)
sim.neg23 <- simulation(I = 1000, 
                        alpha1 = 2.0, alpha2 = 0, alpha3 = 0, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                        sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = -1, lambda3 = -1)

# model misspecification

sim.pos.mis <- simulation(I = 1000, 
                          alpha1 = 2.0, alpha2 = 1.7, alpha3 = 1.6, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                          sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1)
sim.ind.mis <- simulation(I = 1000, 
                          alpha1 = 2.0, alpha2 = 1.7, alpha3 = 1.6, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                          sigma.v = 5, sigma.e = 3, lambda1 = 0, lambda2 = 0, lambda3 = 0)
sim.neg23.mis <- simulation(I = 1000, 
                            alpha1 = 2.0, alpha2 = 1.7, alpha3 = 1.6, beta2 = 1.8, beta3 = 2.4, gamma2 = 0.8, gamma3 = 1.2, 
                            sigma.v = 5, sigma.e = 3, lambda1 = 1, lambda2 = -1, lambda3 = -1)

## results: tables and figures

# summary table: for sample size 200 and 400 only

mySum <- function(mean, sd){
  paste(round(mean, 2), " (", round(sd, 2), ")", sep = "")
}
tabRes <- function(sim){
  pci.200 <- sim$pci[sim$pci$sample.size == 200 & (sim$pci$method == "sqgee" | sim$pci$method == "mqgee") & sim$pci$stage == 1, ]
  rmse2.200 <- sim$rmse2[sim$rmse2$sample.size == 200, ]
  rmse3.200 <- sim$rmse3[sim$rmse3$sample.size == 200, ]
  pci.400 <- sim$pci[sim$pci$sample.size == 400 & (sim$pci$method == "sqgee" | sim$pci$method == "mqgee") & sim$pci$stage == 1, ]
  rmse2.400 <- sim$rmse2[sim$rmse2$sample.size == 400, ]
  rmse3.400 <- sim$rmse3[sim$rmse3$sample.size == 400, ]
  res.200 <- c(round(pci.200$pci[1], 3), mySum(rmse2.200$rmse.mean[1], rmse2.200$rmse.sd[1]), mySum(rmse3.200$rmse.mean[1], rmse3.200$rmse.sd[1]), 
               round(pci.200$pci[2], 3), mySum(rmse2.200$rmse.mean[2], rmse2.200$rmse.sd[2]), mySum(rmse3.200$rmse.mean[2], rmse3.200$rmse.sd[2]))
  res.400 <- c(round(pci.400$pci[1], 3), mySum(rmse2.400$rmse.mean[1], rmse2.400$rmse.sd[1]), mySum(rmse3.400$rmse.mean[1], rmse3.400$rmse.sd[1]), 
               round(pci.400$pci[2], 3), mySum(rmse2.400$rmse.mean[2], rmse2.400$rmse.sd[2]), mySum(rmse3.400$rmse.mean[2], rmse3.400$rmse.sd[2]))
  res <- rbind.data.frame(res.200, res.400)
  colnames(res) <- c("pci.sqgee", "rmse2.sqgee", "rmse3.sqgee",
                     "pci.mqgee", "rmse2.mqgee", "rmse3.mqgee")
  return(res)
}
res <- rbind.data.frame(tabRes(sim.pos), tabRes(sim.ind), tabRes(sim.neg23), 
                        tabRes(sim.pos.mis), tabRes(sim.ind.mis), tabRes(sim.neg23.mis))
kable(res, "latex", longtable = T, booktabs = T, caption = "SimRes") # need knitr package

# figure: plot of PCI at stage 1 and 2
pci <- cbind.data.frame(rbind.data.frame(sim.pos$pci, sim.ind$pci, sim.neg23$pci,
                                         sim.pos.mis$pci, sim.ind.mis$pci, sim.neg23.mis$pci),
                        "cor" = rep(c(rep("positive", 80), rep("independent", 80), rep("negative", 80)), 2),
                        "mis" = rep(c("(a)", "(b)"), each = 3*80))
pci$stage <- factor(pci$stage, levels = c(1, 2), labels = c("Stage 1", "Stage 2"))
pci$cor <- factor(pci$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))
png(filename = "pci.png", width = 30, height = 24, units = "cm", res = 300)
ggplot(pci, aes(x = sample.size, y = pci, shape = method, color = method)) +
  geom_line(aes(linetype = method)) + guides(linetype = FALSE) +
  geom_point() +
  labs(title = "", x = "Sample Size", y = "PCI") +
  ylim(0.5, 1) + scale_x_continuous(breaks = seq(0, 500, 100)) +
  scale_shape_discrete(name = "Method",
                       breaks = c("sq", "mq", "sqgee", "mqgee"),
                       labels = c("Standard Q-learning", "Modified Q-learning",
                                  "Standard Q-learning with GEE", "Modified Q-learning with GEE")) +
  scale_color_discrete(name = "Method",
                       breaks = c("sq", "mq", "sqgee", "mqgee"),
                       labels = c("Standard Q-learning", "Modified Q-learning",
                                  "Standard Q-learning with GEE", "Modified Q-learning with GEE")) + 
  facet_grid(cor ~ mis + stage) + 
  theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
dev.off()

# figure: plot of RMSE of stage 1 heterogeneous causal effects at times 1, 2, 3, conditional on estimated stage 2 optimal rules
rmse <- cbind.data.frame(rbind.data.frame(sim.pos$rmse1, sim.pos$rmse2, sim.pos$rmse3,
                                          sim.ind$rmse1, sim.ind$rmse2, sim.ind$rmse3, 
                                          sim.neg23$rmse1, sim.neg23$rmse2, sim.neg23$rmse3, 
                                          sim.pos.mis$rmse1, sim.pos.mis$rmse2, sim.pos.mis$rmse3, 
                                          sim.ind.mis$rmse1, sim.ind.mis$rmse2, sim.ind.mis$rmse3, 
                                          sim.neg23.mis$rmse1, sim.neg23.mis$rmse2, sim.neg23.mis$rmse3), 
                         "time" = rep(c(rep(1, 20), rep(2, 20), rep(3, 20)), 6), 
                         "cor" = rep(c(rep("positive", 60), rep("independent", 60), rep("negative", 60)), 2),
                         "mis" = rep(c("(a)", "(b)"), each = 20*3*3))
rmse$time <- factor(rmse$time, levels = c(1, 2, 3), labels = c("Time 1", "Time 2", "Time 3"))
rmse$cor <- factor(rmse$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))

png(filename = "rmse.png", width = 30, height = 20, units = "cm", res = 300)
pd <- position_dodge(10)
ggplot(rmse, aes(x = sample.size, y = rmse.mean, shape = method, color = method)) +
       geom_errorbar(aes(ymin = rmse.lower, ymax = rmse.upper), width = 3, position = pd)+
       geom_line(aes(linetype = method)) + guides(linetype = FALSE) +
       geom_point(position = pd) +
       labs(title = "", x = "Sample Size", y = "RMSE") +
       ylim(0, 10) + scale_x_continuous(breaks = seq(0, 400, 80)) +
       scale_shape_discrete(name = "Method",
                            breaks = c("sqgee", "mqgee"),
                            labels = c("Standard Q-learning with GEE", "Modified Q-learning with GEE")) +
       scale_color_discrete(name = "Method",
                            breaks = c("sqgee", "mqgee"),
                            labels = c("Standard Q-learning with GEE", "Modified Q-learning with GEE")) +
  facet_grid(cor ~ mis + time) + 
  theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
dev.off()

# figure: plot of bias of stage 1 heterogeneous causal effects at times 1, 2, 3, conditional on true stage 2 optimal rules
bias <- cbind.data.frame(rbind.data.frame(sim.pos$bias1.sqgee, sim.pos$bias1.mqgee, sim.pos$bias2.sqgee, sim.pos$bias2.mqgee, sim.pos$bias3.sqgee, sim.pos$bias3.mqgee, 
                                          sim.ind$bias1.sqgee, sim.ind$bias1.mqgee, sim.ind$bias2.sqgee, sim.ind$bias2.mqgee, sim.ind$bias3.sqgee, sim.ind$bias3.mqgee, 
                                          sim.neg23$bias1.sqgee, sim.neg23$bias1.mqgee, sim.neg23$bias2.sqgee, sim.neg23$bias2.mqgee, sim.neg23$bias3.sqgee, sim.neg23$bias3.mqgee, 
                                          sim.pos.mis$bias1.sqgee, sim.pos.mis$bias1.mqgee, sim.pos.mis$bias2.sqgee, sim.pos.mis$bias2.mqgee, sim.pos.mis$bias3.sqgee, sim.pos.mis$bias3.mqgee, 
                                          sim.ind.mis$bias1.sqgee, sim.ind.mis$bias1.mqgee, sim.ind.mis$bias2.sqgee, sim.ind.mis$bias2.mqgee, sim.ind.mis$bias3.sqgee, sim.ind.mis$bias3.mqgee, 
                                          sim.neg23.mis$bias1.sqgee, sim.neg23.mis$bias1.mqgee, sim.neg23.mis$bias2.sqgee, sim.neg23.mis$bias2.mqgee, sim.neg23.mis$bias3.sqgee, sim.neg23.mis$bias3.mqgee), 
                         "method" = rep(c(rep("sqgee", 70), rep("mqgee", 70)), 3*6), 
                         "time" = rep(c(rep(1, 70*2), rep(2, 70*2), rep(3, 70*2)), 6), 
                         "cor" = rep(c(rep("positive", 70*6), rep("independent", 70*6), rep("negative", 70*6)), 2),
                         "mis" = rep(c("(a)", "(b)"), each = 70*2*3*3))
bias <- bias[bias$sample.size == 200 | bias$sample.size == 400, ]
bias$sample.size <- factor(bias$sample.size, levels = c(200, 400))
bias$method <- factor(bias$method, levels = c("sqgee","mqgee"), labels = c("SQGEE", "MQGEE"))
bias$time <- factor(bias$time, levels = c(1, 2, 3), labels = c("Time 1", "Time 2", "Time 3"))
bias$cor <- factor(bias$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))


png(filename = "bias.png", width = 30, height = 30, units = "cm", res = 300)
pd <- position_dodge(0.1)
ggplot(bias, aes(x = Y0, y = bias.mean, shape = sample.size, color = sample.size)) + 
             geom_errorbar(aes(ymin = bias.lower, ymax = bias.upper), width = 0.1, position = pd) +
             geom_line(aes(linetype = sample.size)) + guides(linetype = FALSE) +
             geom_point(position = pd) +
             labs(title = "", x = "Y0", y = "Bias") +
             ylim(-20, 20) + scale_x_continuous(breaks = seq(-3, 3, 1)) +
             scale_shape_discrete(name = "Sample Size",
                                  breaks = c("200", "400"),
                                  labels = c("200", "400")) +
             scale_color_discrete(name = "Sample Size",
                                  breaks = c("200", "400"),
                                  labels = c("200", "400")) +
  facet_grid(mis + time ~ cor + method) + 
  theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
dev.off()

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

## load help functions
source("function.R")

###################################################
########## by sample size, equal weights ##########
###################################################

## a function to perform the simulation study

simulation <- function(I = 1000, # number of simulations
                       alpha1 = 0.7, alpha2 = 1, alpha3 = 1.2, beta2 = 0, beta3 = 0, gamma2 = 1, gamma3 = 1.5, 
                       sigma.v = 4, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1, ){
  ss <- seq(100, 500, by = 100) # set sample size
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
  Z1.grid <- c(-3, -2, -1, 0, 1, 2, 3)
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
      sq <- Q(dat, method = "SQ", w = w, misspec = misspec)
      mq <- Q(dat, method = "MQ", w = w, misspec = misspec)
      sqgee <- QGEE(dat, method = "SQGEE", corstr = c("unstructured", "exchangeable"), w = w, misspec = misspec)
      mqgee <- QGEE(dat, method = "MQGEE", corstr = c("unstructured", "exchangeable"), w = w, misspec = misspec)
      metrics.sq <- metric(output = sq, w = w, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.mq <- metric(output = mq, w = w, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.sqgee <- metric(output = sqgee, w = w, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.mqgee <- metric(output = mqgee, w = w, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
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
  bias1.sqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias1.mqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias2.sqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias2.mqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias3.sqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias3.mqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ss)), 
                                  "sample.size" = rep(ss, each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  
  return(list(pci = pci, 
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

sumSim <- function(N.sim = 1000, alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                   sigma.v = 3.2, sigma.e = 2.7){
  
  # correlated outcomes
  sim.pos <- simulation(I = N.sim, 
                        alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                        sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                        misspec = FALSE)
  sim.ind <- simulation(I = N.sim, 
                        alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                        sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = 0, lambda2 = 0, lambda3 = 0, 
                        misspec = FALSE)
  sim.neg23 <- simulation(I = N.sim, 
                          alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                          sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = 1, lambda2 = -1, lambda3 = -1, 
                          misspec = FALSE)
  
  # model misspecification
  sim.pos.mis <- simulation(I = N.sim, 
                            alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                            sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                            misspec = TRUE)
  sim.ind.mis <- simulation(I = N.sim, 
                            alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                            sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = 0, lambda2 = 0, lambda3 = 0, 
                            misspec = TRUE)
  sim.neg23.mis <- simulation(I = N.sim, 
                              alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                              sigma.v = sigma.v, sigma.e = sigma.e, lambda1 = 1, lambda2 = -1, lambda3 = -1, 
                              misspec = TRUE)
  
  # summary table (sample size 200)
  mySum <- function(vec){
    paste(round(vec[1], 2), " (", round(vec[2], 2), ")", sep = "")
  }
  tabRes <- function(sim){
    pci1 <- sim$pci[sim$pci$sample.size==200 & sim$pci$stage==1, ]
    pci2 <- sim$pci[sim$pci$sample.size==200 & sim$pci$stage==2, ]
    rmse1 <- sim$rmse1[sim$rmse1$sample.size==200, ]
    rmse2 <- sim$rmse2[sim$rmse2$sample.size==200, ]
    rmse3 <- sim$rmse3[sim$rmse3$sample.size==200, ]
    res <- cbind.data.frame("method" = c("Q", "mQ", "Q-GEE", "mQ-GEE"), 
                            "pci1" = pci1$pci, "pci2" = pci2$pci, 
                            "rmse1" = c(rep(NA, 2), apply(rmse1[, c("rmse.mean", "rmse.sd")], 1, mySum)), 
                            "rmse2" = c(rep(NA, 2), apply(rmse2[, c("rmse.mean", "rmse.sd")], 1, mySum)), 
                            "rmse3" = c(rep(NA, 2), apply(rmse3[, c("rmse.mean", "rmse.sd")], 1, mySum)))
    return(res)
  }
  res <- cbind.data.frame("correlation" = rep(c("positive", "independent", "negative"), each = 4), 
                          rbind.data.frame(tabRes(sim.pos), tabRes(sim.ind), tabRes(sim.neg23)))
  res.mis <- cbind.data.frame("correlation" = rep(c("positive", "independent", "negative"), each = 4), 
                              rbind.data.frame(tabRes(sim.pos.mis), tabRes(sim.ind.mis), tabRes(sim.neg23.mis)))
  
  # figure: plot of PCI at stage 1 and 2
  pci <- cbind.data.frame(rbind.data.frame(sim.pos$pci, sim.ind$pci, sim.neg23$pci,
                                           sim.pos.mis$pci, sim.ind.mis$pci, sim.neg23.mis$pci),
                          "cor" = rep(c(rep("positive", 40), rep("independent", 40), rep("negative", 40)), 2),
                          "mis" = rep(c("(a)", "(b)"), each = 3*40))
  pci$stage <- factor(pci$stage, levels = c(1, 2), labels = c("Stage 1", "Stage 2"))
  pci$cor <- factor(pci$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))
  pci.plot <- ggplot(pci, aes(x = sample.size, y = pci, shape = method, color = method)) +
    geom_line(aes(linetype = method)) + guides(linetype = FALSE) +
    geom_point() +
    labs(title = "", x = "Sample Size", y = "PCI") +
    ylim(0.5, 1) + scale_x_continuous(breaks = seq(0, 500, 100)) +
    scale_shape_discrete(name = "Method",
                         breaks = c("sq", "mq", "sqgee", "mqgee"),
                         labels = c("Composite Q-learning", "Modified composite Q-learning",
                                    "Q-learning with GEE", "Modified Q-learning with GEE")) +
    scale_color_discrete(name = "Method",
                         breaks = c("sq", "mq", "sqgee", "mqgee"),
                         labels = c("Composite Q-learning", "Modified composite Q-learning",
                                    "Q-learning with GEE", "Modified Q-learning with GEE")) + 
    facet_grid(cor ~ mis + stage) + 
    theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
  
  # figure: plot of RMSE of stage 1 heterogeneous causal effects at times 1, 2, 3, conditional on estimated stage 2 optimal rules
  rmse <- cbind.data.frame(rbind.data.frame(sim.pos$rmse1, sim.pos$rmse2, sim.pos$rmse3,
                                            sim.ind$rmse1, sim.ind$rmse2, sim.ind$rmse3, 
                                            sim.neg23$rmse1, sim.neg23$rmse2, sim.neg23$rmse3, 
                                            sim.pos.mis$rmse1, sim.pos.mis$rmse2, sim.pos.mis$rmse3, 
                                            sim.ind.mis$rmse1, sim.ind.mis$rmse2, sim.ind.mis$rmse3, 
                                            sim.neg23.mis$rmse1, sim.neg23.mis$rmse2, sim.neg23.mis$rmse3), 
                           "time" = rep(c(rep(1, 10), rep(2, 10), rep(3, 10)), 6), 
                           "cor" = rep(c(rep("positive", 30), rep("independent", 30), rep("negative", 30)), 2),
                           "mis" = rep(c("(a)", "(b)"), each = 10*3*3))
  rmse$time <- factor(rmse$time, levels = c(1, 2, 3), labels = c("Time 1", "Time 2", "Time 3"))
  rmse$cor <- factor(rmse$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))
  pd <- position_dodge(10)
  rmse.plot <- ggplot(rmse, aes(x = sample.size, y = rmse.mean, shape = method, color = method)) +
    geom_errorbar(aes(ymin = rmse.lower, ymax = rmse.upper), width = 3, position = pd)+
    geom_line(aes(linetype = method)) + guides(linetype = FALSE) +
    geom_point(position = pd) +
    labs(title = "", x = "Sample Size", y = "RMSE") +
    ylim(0, 10) + scale_x_continuous(breaks = seq(0, 500, 100)) +
    scale_shape_discrete(name = "Method",
                         breaks = c("sqgee", "mqgee"),
                         labels = c("Q-learning with GEE", "Modified Q-learning with GEE")) +
    scale_color_discrete(name = "Method",
                         breaks = c("sqgee", "mqgee"),
                         labels = c("Q-learning with GEE", "Modified Q-learning with GEE")) +
    facet_grid(cor ~ mis + time) + 
    theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
  
  # figure: plot of bias of stage 1 heterogeneous causal effects at times 1, 2, 3, conditional on true stage 2 optimal rules
  bias <- cbind.data.frame(rbind.data.frame(sim.pos$bias1.sqgee, sim.pos$bias1.mqgee, sim.pos$bias2.sqgee, sim.pos$bias2.mqgee, sim.pos$bias3.sqgee, sim.pos$bias3.mqgee, 
                                            sim.ind$bias1.sqgee, sim.ind$bias1.mqgee, sim.ind$bias2.sqgee, sim.ind$bias2.mqgee, sim.ind$bias3.sqgee, sim.ind$bias3.mqgee, 
                                            sim.neg23$bias1.sqgee, sim.neg23$bias1.mqgee, sim.neg23$bias2.sqgee, sim.neg23$bias2.mqgee, sim.neg23$bias3.sqgee, sim.neg23$bias3.mqgee, 
                                            sim.pos.mis$bias1.sqgee, sim.pos.mis$bias1.mqgee, sim.pos.mis$bias2.sqgee, sim.pos.mis$bias2.mqgee, sim.pos.mis$bias3.sqgee, sim.pos.mis$bias3.mqgee, 
                                            sim.ind.mis$bias1.sqgee, sim.ind.mis$bias1.mqgee, sim.ind.mis$bias2.sqgee, sim.ind.mis$bias2.mqgee, sim.ind.mis$bias3.sqgee, sim.ind.mis$bias3.mqgee, 
                                            sim.neg23.mis$bias1.sqgee, sim.neg23.mis$bias1.mqgee, sim.neg23.mis$bias2.sqgee, sim.neg23.mis$bias2.mqgee, sim.neg23.mis$bias3.sqgee, sim.neg23.mis$bias3.mqgee), 
                           "method" = rep(c(rep("sqgee", 35), rep("mqgee", 35)), 3*6), 
                           "time" = rep(c(rep(1, 35*2), rep(2, 35*2), rep(3, 35*2)), 6), 
                           "cor" = rep(c(rep("positive", 35*6), rep("independent", 35*6), rep("negative", 35*6)), 2),
                           "mis" = rep(c("(a)", "(b)"), each = 35*2*3*3))
  bias <- bias[bias$sample.size == 200 | bias$sample.size == 400, ]
  bias$sample.size <- factor(bias$sample.size, levels = c(200, 400))
  bias$method <- factor(bias$method, levels = c("sqgee","mqgee"), labels = c("Q-GEE", "mQ-GEE"))
  bias$time <- factor(bias$time, levels = c(1, 2, 3), labels = c("Time 1", "Time 2", "Time 3"))
  bias$cor <- factor(bias$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))
  pd <- position_dodge(0.1)
  bias.plot <- ggplot(bias, aes(x = Z1, y = bias.mean, shape = sample.size, color = sample.size)) + 
    geom_errorbar(aes(ymin = bias.lower, ymax = bias.upper), width = 0.1, position = pd) +
    geom_line(aes(linetype = sample.size)) + guides(linetype = FALSE) +
    geom_point(position = pd) +
    labs(title = "", x = "Z1", y = "Bias") +
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
  
  return(list("res" = res, "res.mis" = res.mis, "pci.plot" = pci.plot, "rmse.plot" = rmse.plot, "bias.plot" = bias.plot))
}

res.simulation <- sumSim(N.sim = 1000, alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                         sigma.v = 3.2, sigma.e = 2.7)

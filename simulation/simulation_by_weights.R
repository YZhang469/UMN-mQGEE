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

################################
########## by weights ##########
################################

## a function to perform the simulation study

simulation <- function(I = 1000, 
                       alpha1 = 0.7, alpha2 = 1, alpha3 = 1.2, beta2 = 0, beta3 = 0, gamma2 = 1, gamma3 = 1.5, 
                       sigma.v = 4, sigma.e = 3, lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                       misspec = FALSE){
  w0 <- c(1/3, 1/3, 1/3)
  w1 <- c(0, 1, 0)
  w2 <- c(0, 0.5, 0.5)
  w3 <- c(0, 0, 1)
  w4 <- c(0.25, 0.5, 0.25)
  w5 <- c(0.25, 0.25, 0.5)
  w6 <- c(0.5, 0.5, 0)
  w7 <- c(0.5, 0.25, 0.25)
  ws <- list(w0, w1, w2, w3, w4, w5, w6, w7)
  pci.sq.s1 <- matrix(NA, nrow = length(ws), ncol = I)
  pci.sq.s2 <- matrix(NA, nrow = length(ws), ncol = I)
  pci.mq.s1 <- matrix(NA, nrow = length(ws), ncol = I)
  pci.mq.s2 <- matrix(NA, nrow = length(ws), ncol = I)
  pci.sqgee.s1 <- matrix(NA, nrow = length(ws), ncol = I)
  pci.sqgee.s2 <- matrix(NA, nrow = length(ws), ncol = I)
  pci.mqgee.s1 <- matrix(NA, nrow = length(ws), ncol = I)
  pci.mqgee.s2 <- matrix(NA, nrow = length(ws), ncol = I)
  rmse1.sqgee <- matrix(NA, nrow = length(ws), ncol = I)
  rmse1.mqgee <- matrix(NA, nrow = length(ws), ncol = I)
  rmse2.sqgee <- matrix(NA, nrow = length(ws), ncol = I)
  rmse2.mqgee <- matrix(NA, nrow = length(ws), ncol = I)
  rmse3.sqgee <- matrix(NA, nrow = length(ws), ncol = I)
  rmse3.mqgee <- matrix(NA, nrow = length(ws), ncol = I)
  Z1.grid <- c(-3, -2, -1, 0, 1, 2, 3)
  bias.mat <- matrix(NA, nrow = length(Z1.grid), ncol = I)
  bias1.sqgee <- rep(list(bias.mat), length(ws))
  bias1.mqgee <- rep(list(bias.mat), length(ws))
  bias2.sqgee <- rep(list(bias.mat), length(ws))
  bias2.mqgee <- rep(list(bias.mat), length(ws))
  bias3.sqgee <- rep(list(bias.mat), length(ws))
  bias3.mqgee <- rep(list(bias.mat), length(ws))
  for (i in 1:I){
    dat <- generateData.joint(N = 200, 
                              alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3, 
                              sigma.v = sigma.v, sigma.e = sigma.e, 
                              lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
    for (s in 1:length(ws)){
      sq <- Q(dat, method = "SQ", w = ws[[s]], misspec = misspec)
      mq <- Q(dat, method = "MQ", w = ws[[s]], misspec = misspec)
      sqgee <- QGEE(dat, method = "SQGEE", corstr = c("unstructured", "exchangeable"), w = ws[[s]], misspec = misspec)
      mqgee <- QGEE(dat, method = "MQGEE", corstr = c("unstructured", "exchangeable"), w = ws[[s]], misspec = misspec)
      metrics.sq <- metric(sq, w = ws[[s]], alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.mq <- metric(mq, w = ws[[s]], alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.sqgee <- metric(sqgee, w = ws[[s]], alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
      metrics.mqgee <- metric(mqgee, w = ws[[s]], alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, beta2 = beta2, beta3 = beta3, gamma2 = gamma2, gamma3 = gamma3)
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
  pci <- cbind.data.frame("weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), 8), 
                          "stage" = rep(c(rep(1, length(ws)), rep(2, length(ws))), 4), 
                          "method" = rep(c("sq", "mq", "sqgee", "mqgee"), each = 2*length(ws)), 
                          "pci" = c(apply(pci.sq.s1, 1, mean), apply(pci.sq.s2, 1, mean), 
                                    apply(pci.mq.s1, 1, mean), apply(pci.mq.s2, 1, mean),
                                    apply(pci.sqgee.s1, 1, mean), apply(pci.sqgee.s2, 1, mean), 
                                    apply(pci.mqgee.s1, 1, mean), apply(pci.mqgee.s2, 1, mean)))
  rmse1 <- cbind.data.frame("weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), 2),
                            "method" = rep(c("sqgee", "mqgee"), each = length(ws)), 
                            "rmse.mean" = c(apply(rmse1.sqgee, 1, mean), apply(rmse1.mqgee, 1, mean)),
                            "rmse.sd" = c(apply(rmse1.sqgee, 1, sd), apply(rmse1.mqgee, 1, sd)),
                            "rmse.lower" = c(apply(rmse1.sqgee, 1, quantile, probs = 0.025), apply(rmse1.mqgee, 1, quantile, probs = 0.025)),
                            "rmse.upper" = c(apply(rmse1.sqgee, 1, quantile, probs = 0.975), apply(rmse1.mqgee, 1, quantile, probs = 0.975)))
  rmse2 <- cbind.data.frame("weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), 2),
                            "method" = rep(c("sqgee", "mqgee"), each = length(ws)), 
                            "rmse.mean" = c(apply(rmse2.sqgee, 1, mean), apply(rmse2.mqgee, 1, mean)),
                            "rmse.sd" = c(apply(rmse2.sqgee, 1, sd), apply(rmse2.mqgee, 1, sd)),
                            "rmse.lower" = c(apply(rmse2.sqgee, 1, quantile, probs = 0.025), apply(rmse2.mqgee, 1, quantile, probs = 0.025)),
                            "rmse.upper" = c(apply(rmse2.sqgee, 1, quantile, probs = 0.975), apply(rmse2.mqgee, 1, quantile, probs = 0.975)))
  rmse3 <- cbind.data.frame("weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), 2),
                            "method" = rep(c("sqgee", "mqgee"), each = length(ws)), 
                            "rmse.mean" = c(apply(rmse3.sqgee, 1, mean), apply(rmse3.mqgee, 1, mean)),
                            "rmse.sd" = c(apply(rmse3.sqgee, 1, sd), apply(rmse3.mqgee, 1, sd)),
                            "rmse.lower" = c(apply(rmse3.sqgee, 1, quantile, probs = 0.025), apply(rmse3.mqgee, 1, quantile, probs = 0.025)),
                            "rmse.upper" = c(apply(rmse3.sqgee, 1, quantile, probs = 0.975), apply(rmse3.mqgee, 1, quantile, probs = 0.975)))
  bias1.sqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ws)), 
                                  "weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias1.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias1.mqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ws)), 
                                  "weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias1.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias2.sqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ws)), 
                                  "weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias2.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias2.mqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ws)), 
                                  "weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias2.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias3.sqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ws)), 
                                  "weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias3.sqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  bias3.mqgee <- cbind.data.frame("Z1" = rep(Z1.grid, length(ws)), 
                                  "weight" = rep(unlist(lapply(ws, function(vec){paste0("(", paste(round(vec, 2), collapse = ","), ")")})), each = length(Z1.grid)), 
                                  "bias.mean" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, mean)})), 
                                  "bias.lower" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, quantile, probs = 0.025)})), 
                                  "bias.upper" = unlist(lapply(bias3.mqgee, function(df){apply(df, 1, quantile, probs = 0.975)})))
  
  return(list(pci = pci, 
              rmse1 = rmse1, rmse2 = rmse2, rmse3 = rmse3,
              bias1.sqgee = bias1.sqgee, bias1.mqgee = bias1.mqgee,
              bias2.sqgee = bias2.sqgee, bias2.mqgee = bias2.mqgee,
              bias3.sqgee = bias3.sqgee, bias3.mqgee = bias3.mqgee))
}

set.seed(2021)

# correlated outcomes
sim.pos <- simulation(I = 1000, 
                      alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                      sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = 1, lambda3 = 1, misspec = FALSE)
sim.ind <- simulation(I = 1000, 
                      alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                      sigma.v = 3.2, sigma.e = 2.7, lambda1 = 0, lambda2 = 0, lambda3 = 0, misspec = FALSE)
sim.neg23 <- simulation(I = 1000, 
                        alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                        sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = -1, lambda3 = -1, misspec = FALSE)

# model misspecification
sim.pos.mis <- simulation(I = 1000, 
                          alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                          sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = 1, lambda3 = 1, misspec = TRUE)
sim.ind.mis <- simulation(I = 1000, 
                          alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                          sigma.v = 3.2, sigma.e = 2.7, lambda1 = 0, lambda2 = 0, lambda3 = 0, misspec = TRUE)
sim.neg23.mis <- simulation(I = 1000, 
                            alpha1 = -2.0, alpha2 = -1.7, alpha3 = -1.6, beta2 = 0, beta3 = 0, gamma2 = -1.2, gamma3 = -0.8, 
                            sigma.v = 3.2, sigma.e = 2.7, lambda1 = 1, lambda2 = -1, lambda3 = -1, misspec = TRUE)

pci <- cbind.data.frame(rbind.data.frame(sim.pos$pci, sim.ind$pci, sim.neg23$pci,
                                         sim.pos.mis$pci, sim.ind.mis$pci, sim.neg23.mis$pci),
                        "cor" = rep(c(rep("positive", 64), rep("independent", 64), rep("negative", 64)), 2),
                        "mis" = rep(c("(a)", "(b)"), each = 3*64))
pci$stage <- factor(pci$stage, levels = c(1, 2), labels = c("Stage 1", "Stage 2"))
pci$cor <- factor(pci$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))
ggplot(pci, aes(x = weight, y = pci, shape = method, color = method)) +
  geom_point() +
  labs(title = "", x = "Weight", y = "PCI") +
  ylim(0, 1) + 
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
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rmse <- cbind.data.frame(rbind.data.frame(sim.pos$rmse1, sim.pos$rmse2, sim.pos$rmse3,
                                          sim.ind$rmse1, sim.ind$rmse2, sim.ind$rmse3, 
                                          sim.neg23$rmse1, sim.neg23$rmse2, sim.neg23$rmse3, 
                                          sim.pos.mis$rmse1, sim.pos.mis$rmse2, sim.pos.mis$rmse3, 
                                          sim.ind.mis$rmse1, sim.ind.mis$rmse2, sim.ind.mis$rmse3, 
                                          sim.neg23.mis$rmse1, sim.neg23.mis$rmse2, sim.neg23.mis$rmse3), 
                         "time" = rep(c(rep(1, 16), rep(2, 16), rep(3, 16)), 6), 
                         "cor" = rep(c(rep("positive", 48), rep("independent", 48), rep("negative", 48)), 2),
                         "mis" = rep(c("(a)", "(b)"), each = 16*3*3))
rmse$time <- factor(rmse$time, levels = c(1, 2, 3), labels = c("Time 1", "Time 2", "Time 3"))
rmse$cor <- factor(rmse$cor, levels = c("positive", "independent", "negative"), labels = c("Positive", "Independent", "Negative"))
ggplot(rmse, aes(x = weight, y = rmse.mean, shape = method, color = method)) +
  geom_errorbar(aes(ymin = rmse.lower, ymax = rmse.upper), width = 0.5)+
  geom_point() +
  labs(title = "", x = "Weight", y = "RMSE") +
  ylim(0, 7.5) + 
  scale_shape_discrete(name = "Method",
                       breaks = c("sqgee", "mqgee"),
                       labels = c("Q-learning with GEE", "Modified Q-learning with GEE")) +
  scale_color_discrete(name = "Method",
                       breaks = c("sqgee", "mqgee"),
                       labels = c("Q-learning with GEE", "Modified Q-learning with GEE")) +
  facet_grid(cor ~ mis + time) + 
  theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

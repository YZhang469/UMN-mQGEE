### data analysis using Q-learning with GEE

## load data

mbridge <- read.csv("simulated_data.csv")

## load packages

library(geepack)
library(ggplot2)
library(dplyr)
library(forcats)
library(knitr)
library(ggpubr)

# function for variable selection

gee_stepper <- function(fit) {
  preds     <- attributes(terms(formula(fit)))$term.labels
  preds_out <- preds
  preds_in  <- character(0)
  fit0 <- stats::update(fit, formula = . ~ 1, data = fit$data, corstr = fit$corstr)
  while(TRUE) {
    fits <-
      c(list(fit0),
        lapply(preds_out, function(p) {stats::update(fit0, formula = stats::as.formula(paste(". ~ . +", p)), data = fit$data, corstr = fit$corstr) }),
        lapply(preds_in,  function(p) {stats::update(fit0, formula = stats::as.formula(paste(". ~ . -", p)), data = fit$data, corstr = fit$corstr) })
      )
    minqic <- which.min(sapply(fits, MESS::QIC)[1, ])
    fit0 <- fits[[minqic]]
    if (minqic <= 1+length(preds_out)){
      preds_in <- c(preds_in, preds_out[minqic-1])
      preds_out <- dplyr::setdiff(preds, preds_in)
    }
    if (minqic > 1+length(preds_out)){
      preds_out <- c(preds_out, preds_in[minqic-1-length(preds_out)])
      preds_in <- dplyr::setdiff(preds, preds_out)
    }
    if (minqic == 1) {
      break
    }
  }
  # force the time-dependent effects of A1 in stage 1 model
  if (!is.element("A2", preds)){
    if (!is.element("factor(time):A1", preds_in)){
      fit0 <- update(fit0, formula = stats::as.formula(paste(". ~ . +", "factor(time):A1")), data = fit$data, corstr = fit$corstr)
    }
  }
  # force the time-dependent effects of A2 and A1*A2 in stage 2 model
  if (is.element("A2", preds)){
    if (!is.element("A1:A2", preds_in)){
      fit0 <- update(fit0, formula = stats::as.formula(paste(". ~ . +", "A1:A2")), data = fit$data, corstr = fit$corstr)
    }
    if (!is.element("factor(time):A2", preds_in)){
      fit0 <- update(fit0, formula = stats::as.formula(paste(". ~ . +", "factor(time):A2")), data = fit$data, corstr = fit$corstr)
    }
    if (!is.element("factor(time):A1:A2", preds_in)){
      fit0 <- update(fit0, formula = stats::as.formula(paste(". ~ . +", "factor(time):A1:A2")), data = fit$data, corstr = fit$corstr)
    }
  }
  fit0
}

# note that M-bridge is a SMART with embedded tailoring, so the data used to estimate stage 2 Q-function should be restricted to heavy-drinkers: A2 in {-1,1} only if R = 0

analyzeData <- function(dat, 
                        method = "MQGEE", # method can be “SQGEE” or “MQGEE”
                        corstr = c("unstructured", "exchangeable"), 
                        w = c(1/3, 1/3, 1/3)){
  
  # construct stage 2 data (including covariates at stage 1)
  n <- nrow(dat)
  dat.s2 <- cbind.data.frame("ID" = rep(dat$subject, 2), 
                             do.call(rbind, replicate(2, dat[ , c("nonwhite", "female", "yes_greek", "num_days0", "avg_drinks0", "binge0", 
                                                                  "A1", "sm_binge", "HD", "A2")], simplify = FALSE)), 
                             "time" = rep(c(2, 3), each = n), "Y" = c(dat$binge1, dat$binge2))
  dat.s2 <- dat.s2[order(dat.s2$ID), ] # the original data set under actual treatment received at stage 1 and 2
  
  # stage 2 optimal treatment and outcome
  rand.sample <- dat.s2$HD == 1
  mod.s2 <- geeglm(Y ~ (factor(time)*(nonwhite + female + yes_greek + num_days0 + avg_drinks0 + binge0 + A1 + sm_binge)) * A2, 
                   data = dat.s2[rand.sample, ], id = ID, corstr = corstr[2])
  mod.s2 <- geeglm(gee_stepper(mod.s2)$formula, data = dat.s2[rand.sample, ], id = ID, corstr = corstr[2])
  # optimal treatment at stage 2: minimize an weighted average of estimated binge1 and binge2 under A2 = (-1,1)
  dat.s2.opt <- dat.s2 # create a data set to store estimated optimal treatment at stage 2
  # A2 = -1
  dat.s2a <- dat.s2[rand.sample, ]
  dat.s2a$A2 <- -1
  dat.s2a$Y <- predict(mod.s2, newdata = dat.s2a)
  # A2 = 1
  dat.s2b <- dat.s2[rand.sample, ]
  dat.s2b$A2 <- 1
  dat.s2b$Y <- predict(mod.s2, newdata = dat.s2b)
  # identify stage 2 optimal rule
  Ya <- apply(matrix(dat.s2a$Y, 2), 2, weighted.mean, w = w[2:3])
  Yb <- apply(matrix(dat.s2b$Y, 2), 2, weighted.mean, w = w[2:3])
  dat.s2.opt$A2[rand.sample] <- rep(ifelse(Ya < Yb, -1, 1), each = 2)
  # compare the accumulative outcome of Y2 and Y3 (this algorithm assumes equal weights)
  if (method == "SQGEE"){
    dat.s2.opt$Y[rand.sample] <- predict(mod.s2, newdata = dat.s2.opt[rand.sample, ])
  }
  else if (method == "MQGEE"){
    int2 <- (dat.s2b$Y - dat.s2a$Y)/2
    dat.s2.opt$Y[rand.sample] <- dat.s2$Y[rand.sample] + 2 * ifelse(dat.s2.opt$A2[rand.sample] == dat.s2$A2[rand.sample], 0, 1) * dat.s2.opt$A2[rand.sample] * int2
  }
  
  # construct stage 1 data
  dat.s2.opt <- dat.s2.opt[order(dat.s2.opt$time), ]
  dat.s1 <- cbind.data.frame("ID" = rep(dat.s2.opt$ID[1:n], 3), 
                             do.call(rbind, replicate(3, dat.s2.opt[1:n, c("nonwhite", "female", "yes_greek", "num_days0", "avg_drinks0", "binge0", 
                                                                           "A1", "HD")], simplify = FALSE)), 
                             "A2" = rep(dat.s2.opt$A2[1:n], 3),
                             "time" = rep(c(1, 2, 3), each = n), 
                             "Y" = c(dat$sm_binge, dat.s2.opt$Y))
  dat.s1 <- dat.s1[order(dat.s1$ID), ] # the data set under actual treatment received at stage 1 and optimal treatment received at stage 2
  
  # stage 1 optimal treatment and outcome
  mod.s1 <- geeglm(Y ~ (factor(time)*(nonwhite + female + yes_greek + num_days0 + avg_drinks0 + binge0)) * A1, 
                   data = dat.s1, id = ID, corstr = corstr[1])
  mod.s1 <- geeglm(gee_stepper(fit = mod.s1)$formula, data = dat.s1, id = ID, corstr = corstr[1])
  # optimal treatment at stage 1: optimize the average of Y1, Y2 and Y3 under A1 = (-1,1)
  dat.s1.opt <- dat.s1 # create a data set to store estimated optimal treatment at stage 1
  # A2 = -1
  dat.s1a <- dat.s1
  dat.s1a$A1 <- -1
  dat.s1a$Y <- predict(mod.s1, newdata = dat.s1a)
  # A2 = 1
  dat.s1b <- dat.s1
  dat.s1b$A1 <- 1
  dat.s1b$Y <- predict(mod.s1, newdata = dat.s1b)
  Ya <- apply(matrix(dat.s1a$Y, 3), 2, weighted.mean, w = w[1:3])
  Yb <- apply(matrix(dat.s1b$Y, 3), 2, weighted.mean, w = w[1:3])
  dat.s1.opt$A1 <- rep(ifelse(Ya < Yb, -1, 1), each = 3)
  # compare the average outcome at time 1, 2 and 3 (this algorithm assumes equal weights)
  dat.s1.opt$Y <- predict(mod.s1, newdata = dat.s1.opt)
  if (method == "SQGEE"){
    dat.s1.opt$Y <- predict(mod.s1, newdata = dat.s1.opt)
  }
  else if (method == "MQGEE"){
    int1 <- (dat.s1b$Y - dat.s1a$Y)/2
    dat.s1.opt$Y <- dat.s1$Y + 2 * ifelse(dat.s1.opt$A1 == dat.s1$A1, 0, 1) * dat.s1.opt$A1 * int1
  }
  
  dat.opt <- cbind.data.frame(dat.s1.opt[ , c("ID", "nonwhite", "female", "yes_greek", "num_days0", "avg_drinks0", "binge0", 
                                              "A1", "HD", "A2", "time", "Y")])
  # convert dat.opt from long format to wide format
  dat.opt <- reshape(dat.opt, v.names = "Y", idvar = "ID", timevar = "time", direction = "wide", sep = "")
  colnames(dat.opt)[which(colnames(dat.opt) == "binge0")] <- "Y0"
  
  return(list("mod.s1" = mod.s1, "mod.s2" = mod.s2, "dat.opt" = dat.opt))
}

## estimate the (optimal or specific) DTRs as a by-product of Q-learning

estTraj <- function(dat, a1, a2, corstr){
  n <- nrow(dat)
  dat.s2 <- cbind.data.frame("ID" = rep(dat$subject, 2), 
                             do.call(rbind, replicate(2, dat[ , c("nonwhite", "female", "yes_greek", "num_days0", "avg_drinks0", "binge0", 
                                                                  "A1", "sm_binge", "HD", "A2")], simplify = FALSE)), 
                             "time" = rep(c(2, 3), each = n), "Y" = c(dat$binge1, dat$binge2))
  dat.s2 <- dat.s2[order(dat.s2$ID), ]
  rand.sample <- dat.s2$HD == 1
  mod.s2 <- geeglm(Y ~ (factor(time)*(nonwhite + female + yes_greek + num_days0 + avg_drinks0 + binge0 + A1 + sm_binge)) * A2,
                   data = dat.s2[rand.sample, ], id = ID, corstr = corstr[2])
  mod.s2 <- geeglm(gee_stepper(mod.s2)$formula,
                   data = dat.s2[rand.sample, ], id = ID, corstr = corstr[2])
  dat.s2.est <- dat.s2
  dat.s2.est$A2[rand.sample] <- a2
  dat.s2.est$Y[rand.sample] <- predict(mod.s2, newdata = dat.s2.est[rand.sample,])
  
  # stage 1 estimation assuming stage 1 treatment a1
  dat.s2.est <- dat.s2.est[order(dat.s2.est$time), ]
  dat.s1 <- cbind.data.frame("ID" = rep(dat.s2.est$ID[1:n], 3), 
                             do.call(rbind, replicate(3, dat.s2.est[1:n, c("nonwhite", "female", "yes_greek", "num_days0", "avg_drinks0", "binge0", 
                                                                           "A1", "HD")], simplify = FALSE)), 
                             "A2" = rep(dat.s2.est$A2[1:n], 3),
                             "time" = rep(c(1, 2, 3), each = n), 
                             "Y" = c(dat$sm_binge, dat.s2.est$Y))
  dat.s1 <- dat.s1[order(dat.s1$ID), ]
  mod.s1 <- geeglm(Y ~ (factor(time)*(nonwhite + female + yes_greek + num_days0 + avg_drinks0 + binge0)) * A1, 
                   data = dat.s1, id = ID, corstr = corstr[1])
  mod.s1 <- geeglm(gee_stepper(mod.s1)$formula,
                   data = dat.s1, id = ID, corstr = corstr[1])
  dat.s1.est <- dat.s1
  dat.s1.est$A1 <- a1
  dat.s1.est$Y <- predict(mod.s1, newdata = dat.s1.est)
  
  dat.est <- cbind.data.frame(dat.s1.est[ , c("ID", "nonwhite", "female", "yes_greek", "num_days0", "avg_drinks0", "binge0", 
                                              "A1", "HD", "A2", "time", "Y")])
  # convert dat.opt from long format to wide format
  dat.est <- reshape(dat.est, v.names = "Y", idvar = "ID", timevar = "time", direction = "wide", sep = "")
  colnames(dat.est)[which(colnames(dat.est) == "binge0")] <- "Y0"
  
  return(list("mod.s1" = mod.s1, "mod.s2" = mod.s2, "dat.est" = dat.est))
}

## summary and plots of analysis results

# plot the trajectory of DTRs

sumTraj <- function(dat, dat.opt, corstr = c("unstructured", "exchangeable")){
  dat.est1 <- estTraj(dat, a1 = 1, a2 = -1, corstr)$dat.est
  dat.est2 <- estTraj(dat, a1 = 1, a2 = 1, corstr)$dat.est
  dat.est3 <- estTraj(dat, a1 = -1, a2 = -1, corstr)$dat.est
  dat.est4 <- estTraj(dat, a1 = -1, a2 = 1, corstr)$dat.est
  y.opt <- c(mean(dat.opt$Y0, na.rm = T), mean(dat.opt$Y1, na.rm = T), mean(dat.opt$Y2, na.rm = T), mean(dat.opt$Y3, na.rm = T))
  y.est1 <- c(mean(dat.est1$Y0, na.rm = T), mean(dat.est1$Y1, na.rm = T), mean(dat.est1$Y2, na.rm = T), mean(dat.est1$Y3, na.rm = T))
  y.est2 <- c(mean(dat.est2$Y0, na.rm = T), mean(dat.est2$Y1, na.rm = T), mean(dat.est2$Y2, na.rm = T), mean(dat.est2$Y3, na.rm = T))
  y.est3 <- c(mean(dat.est3$Y0, na.rm = T), mean(dat.est3$Y1, na.rm = T), mean(dat.est3$Y2, na.rm = T), mean(dat.est3$Y3, na.rm = T))
  y.est4 <- c(mean(dat.est4$Y0, na.rm = T), mean(dat.est4$Y1, na.rm = T), mean(dat.est4$Y2, na.rm = T), mean(dat.est4$Y3, na.rm = T))
  pred <- cbind.data.frame("time" = rep(c("Baseline", "Self-monitoring", "Follow-up 1", "Follow-up 2"), 5),
                           "y" = c(y.opt, y.est1, y.est2, y.est3, y.est4), 
                           "DTR" = c(rep("Optimal DTR", 4), rep("Early/Email", 4), rep("Early/Coach", 4), rep("Late/Email", 4), rep("Late/Coach", 4)))
  return(pred)
}

output_equal_mqgee <- analyzeData(dat = mbridge, method = "MQGEE")
pred_equal_mqgee <- sumTraj(dat = mbridge, dat.opt = output_equal_mqgee$dat.opt)

pred_equal_mqgee %>%
  mutate(DTR = fct_relevel(DTR, "Optimal DTR", "Early/Email", "Early/Coach", "Late/Email", "Late/Coach")) %>%
  mutate(time = fct_relevel(factor(time), "Baseline", "Self-monitoring", "Follow-up 1", "Follow-up 2")) %>%
  ggplot(aes(x = time, y = y, shape = DTR, linetype = DTR)) + 
  geom_point() + geom_line(aes(group = DTR)) + 
  labs(title = "", x = "", y = "Average frequency of binge drinking") + 
  ylim(0,2.1) + 
  theme(legend.title = element_blank())
dev.off()

# boxplot of estimated CATEs (stage k, time j)

calcCATE <- function(dat, output){
  dat1.temp <- cbind.data.frame("ID" = rep(dat$subject, 3), 
                                do.call(rbind, replicate(3, dat[ , c("nonwhite", "female", "yes_greek", 
                                                                     "num_days0", "avg_drinks0", "binge0", 
                                                                     "A1")], simplify = FALSE)), 
                                "time" = rep(c(1, 2, 3), each = nrow(dat)))
  dat1.temp <- dat1.temp[order(dat1.temp$ID), ]
  dat1a.temp <- dat1.temp
  dat1a.temp$A1 <- -1
  dat1b.temp <- dat1.temp
  dat1b.temp$A1 <- 1
  delta1.hat <- matrix(predict(output$mod.s1, newdata = dat1b.temp), ncol = 3, byrow = TRUE) - matrix(predict(output$mod.s1, newdata = dat1a.temp), ncol = 3, byrow = TRUE)
  dat2.temp <- cbind.data.frame("ID" = rep(dat$subject, 2), 
                                do.call(rbind, replicate(2, dat[ , c("nonwhite", "female", "yes_greek", "num_days0", "avg_drinks0", "binge0", 
                                                                     "A1", "sm_binge", "HD", "A2")], simplify = FALSE)), 
                                "time" = rep(c(2, 3), each = nrow(dat)))
  dat2.temp <- dat2.temp[order(dat2.temp$ID), ]
  dat2a.temp <- dat2.temp
  dat2a.temp$A2[dat2a.temp$HD == 1] <- -1
  dat2b.temp <- dat2.temp
  dat2b.temp$A2[dat2a.temp$HD == 1] <- 1
  delta2.hat <- matrix(predict(output$mod.s2, newdata = dat2b.temp), ncol = 2, byrow = TRUE) - matrix(predict(output$mod.s2, newdata = dat2a.temp), ncol = 2, byrow = TRUE)
  CATE <- cbind.data.frame(delta1.hat, delta2.hat)
  colnames(CATE) <- c("s1t1", "s1t2", "s1t3", "s2t2", "s2t3")
  return(CATE)
}

CATE <- calcCATE(dat = mbridge, output = output_equal_mqgee)
CATE$s2t2[CATE$s2t2 == 0] <- NA
CATE$s2t3[CATE$s2t3 == 0] <- NA

library(reshape2)
CATE1.temp <- melt(CATE[, 1:3])
CATE2.temp <- melt(CATE[, 4:5])
p1 <- ggplot(data = CATE1.temp, aes(x = variable, y = value)) + 
  geom_boxplot() + stat_summary(fun = "mean", shape = 15) +
  labs(title= "Stage 1 intervention", x = "", y = "Heterogeneous causal effect") + 
  ylim(-2, 2) + 
  scale_x_discrete(labels = c("Self-monitoring", "Follow-up 1", "Follow-up 2"))
p2 <- ggplot(data = CATE2.temp, aes(x = variable, y = value)) + 
  geom_boxplot() + stat_summary(fun = "mean", shape = 15) +
  labs(title= "Stage 2 intervention", x = "", y = "Heterogeneous causal effect") + 
  ylim(-4, 4) + 
  scale_x_discrete(labels = c("Follow-up 1", "Follow-up 2"))
ggarrange(p1, p2, ncol = 2, nrow = 1)

########################################
########## Additional Results ##########
########################################

# outputs under unequal weights with different methods

w0 <- c(1/3, 1/3, 1/3)
w1 <- c(0, 1, 0)
w2 <- c(0, 0.5, 0.5)
w3 <- c(0, 0, 1)
w4 <- c(0.25, 0.5, 0.25)
w5 <- c(0.25, 0.25, 0.5)
w6 <- c(0.5, 0.5, 0)
w7 <- c(0.5, 0.25, 0.25)
ws <- list(w0, w1, w2, w3, w4, w5, w6, w7)

perc_ws_method <- data.frame()
traj_ws_method <- data.frame()
for (i in 1:length(ws)){
  output.sqgee <- analyzeData(dat = mbridge, method = "SQGEE", w = ws[[i]])
  output.mqgee <- analyzeData(dat = mbridge, method = "MQGEE", w = ws[[i]])
  perc.sqgee <- c(table(factor(output.sqgee$dat.opt$A1,levels = c(-1,1)))/nrow(mbridge), table(factor(output.sqgee$dat.opt$A2,levels = c(-1,1)))/sum(mbridge$HD))
  perc.mqgee <- c(table(factor(output.mqgee$dat.opt$A1,levels = c(-1,1)))/nrow(mbridge), table(factor(output.mqgee$dat.opt$A2,levels = c(-1,1)))/sum(mbridge$HD))
  perc_ws_method <- rbind.data.frame(perc_ws_method, 
                                     c(paste0("(", paste(round(ws[[i]], 2), collapse = ","), ")"), "SQGEE", round(perc.sqgee, 3)), 
                                     c(paste0("(", paste(round(ws[[i]], 2), collapse = ","), ")"), "MQGEE", round(perc.mqgee, 3)))
  y.opt.sqgee <- c(mean(output.sqgee$dat.opt$Y0, na.rm = T), mean(output.sqgee$dat.opt$Y1, na.rm = T), 
                   mean(output.sqgee$dat.opt$Y2, na.rm = T), mean(output.sqgee$dat.opt$Y3, na.rm = T))
  y.opt.mqgee <- c(mean(output.mqgee$dat.opt$Y0, na.rm = T), mean(output.mqgee$dat.opt$Y1, na.rm = T), 
                   mean(output.mqgee$dat.opt$Y2, na.rm = T), mean(output.mqgee$dat.opt$Y3, na.rm = T))
  traj_ws_method <- rbind.data.frame(traj_ws_method, 
                                     c(paste0("(", paste(round(ws[[i]], 2), collapse = ","), ")"), "SQGEE", y.opt.sqgee), 
                                     c(paste0("(", paste(round(ws[[i]], 2), collapse = ","), ")"), "MQGEE", y.opt.mqgee))
}
colnames(perc_ws_method) <- c("w", "method", "A1 = -1", "A1 = 1", "A2 = -1", "A2 = 1")
colnames(traj_ws_method) <- c("w", "method", "Y0", "Y1", "Y2", "Y3")

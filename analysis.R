### data analysis using Q-learning with GEE

## load data

dat <- read.table("http://supp.apa.org/psycarticles/supplemental/met0000219/5.TutorialDatasetWideFormat.txt", header = TRUE)

## load packages

library(geepack)
library(ggplot2)
library(dplyr)
library(forcats)

## use Q-learning algorithm with GEE to analyze a simulated data set generated from a simple version of ENGAGE study
# note that ENGAGE is a tailored SMART, so the model for stage 2 Q-function should be stratified by responders and non-responders: A2 in {-1,1} only if R = 0

analyzeData <- function(dat, 
                        method = "MQGEE", # method can be “SQGEE” or “MQGEE”
                        corstr = c("exchangeable", "unstructured")){  # if method includes GEE, then need to specify correlation structure
  
  # construct stage 2 data (including covariates at stage 1)
  n <- nrow(dat)
  dat.s2 <- cbind.data.frame(do.call(rbind, replicate(2, dat[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1", "Y1", "R", "A2")], simplify = FALSE)), 
                             "time" = rep(c(2, 3), each = n), "Y" = c(dat$Y2, dat$Y3))
  dat.s2 <- dat.s2[order(dat.s2$SubjectID), ] # the original data set under actual treatment received at stage 1 and 2
  
  # stage 2 optimal treatment and outcome
  rand.sample <- dat.s2$R==0 # patients who were randomized at stage 2
  mod.s2 <- geeglm(Y ~ factor(time) + CenteredAge + Male + Education + HDDays + A1 + Y1 + factor(time)*(CenteredAge + Male + Education + HDDays + A1 + Y1) + 
                     A2 + A2*(factor(time) + CenteredAge + Male + Education + HDDays + A1 + Y1 + factor(time)*(CenteredAge + Male + Education + HDDays + A1 + Y1)),
                   data = dat.s2[rand.sample, ], id = SubjectID, corstr = corstr[1])
  # optimal treatment at stage 2: optimize the accumulative outcome of Y2 and Y3 under A2 = (-1,1)
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
  Ya <- colMeans(matrix(dat.s2a$Y, 2))
  Yb <- colMeans(matrix(dat.s2b$Y, 2))
  dat.s2.opt$A2[rand.sample] <- rep(ifelse(Ya > Yb, -1, 1), each = 2)
  # compare the accumulative outcome of Y2 and Y3 (this algorithm assumes equal weights)
  if (method == "SQGEE"){
    dat.s2.opt$Y[rand.sample] <- predict(mod.s2, newdata = dat.s2.opt[rand.sample, ])
  }
  else if (method == "MQGEE"){
    int2 <- as.matrix(cbind.data.frame(dat.s2$A2, dat.s2$A2*ifelse(dat.s2$time==3, 1, 0), 
                                       dat.s2$A2*dat.s2$CenteredAge, dat.s2$A2*dat.s2$Male, dat.s2$A2*dat.s2$Education, 
                                       dat.s2$A2*dat.s2$HDDays, dat.s2$A2*dat.s2$A1, dat.s2$A2*dat.s2$Y1, 
                                       dat.s2$A2*dat.s2$CenteredAge*ifelse(dat.s2$time==3, 1, 0), dat.s2$A2*dat.s2$Male*ifelse(dat.s2$time==3, 1, 0), dat.s2$A2*dat.s2$Education*ifelse(dat.s2$time==3, 1, 0), 
                                       dat.s2$A2*dat.s2$HDDays*ifelse(dat.s2$time==3, 1, 0), dat.s2$A2*dat.s2$A1*ifelse(dat.s2$time==3, 1, 0), dat.s2$A2*dat.s2$Y1*ifelse(dat.s2$time==3, 1, 0)))[rand.sample, ] %*% coef(mod.s2)[c(9, 16:28)]
    dat.s2.opt$Y[rand.sample] <- dat.s2$Y[rand.sample] + 2 * ifelse(dat.s2.opt$A2[rand.sample] == dat.s2$A2[rand.sample], 0, 1) * abs(int2)
  }
  
  # construct stage 1 data
  dat.s2.opt <- dat.s2.opt[order(dat.s2.opt$time), ]
  dat.s1 <- cbind.data.frame(do.call(rbind, replicate(3, dat[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1", "Y1", "R")], simplify = FALSE)), 
                             "A2" = rep(dat.s2.opt$A2[1:n], 3),
                             "time" = rep(c(1, 2, 3), each = n), 
                             "Y" = c(dat$Y1, dat.s2.opt$Y))
  dat.s1 <- dat.s1[order(dat.s1$SubjectID), ] # the data set under actual treatment received at stage 1 and optimal treatment received at stage 2
  # stage 1 optimal treatment and outcome
  # stage 1 model
  # use Caucasian only to avoid rank deficiency
  mod.s1 <- geeglm(Y ~ factor(time) + CenteredAge + Male + Education + HDDays + factor(time)*(CenteredAge + Male + Education + HDDays) + 
                     A1 + A1 * (factor(time) + CenteredAge + Male + Education + HDDays + factor(time)*(CenteredAge + Male + Education + HDDays)), 
                   data = dat.s1, id = SubjectID, corstr = corstr[2])
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
  Ya <- colMeans(matrix(dat.s1a$Y, 3))
  Yb <- colMeans(matrix(dat.s1b$Y, 3))
  dat.s1.opt$A1 <- rep(ifelse(Ya > Yb, -1, 1), each = 3)
  # compare the average outcome at time 1, 2 and 3 (this algorithm assumes equal weights)
  dat.s1.opt$Y <- predict(mod.s1, newdata = dat.s1.opt)
  if (method == "SQGEE"){
    dat.s1.opt$Y <- predict(mod.s1, newdata = dat.s1.opt)
  }
  else if (method == "MQGEE"){
    int1 <- as.matrix(cbind.data.frame(dat.s1$A1, dat.s1$A1*ifelse(dat.s1$time==2, 1, 0), dat.s1$A1*ifelse(dat.s1$time==3, 1, 0), 
                                       dat.s1$A1*dat.s1$CenteredAge, dat.s1$A1*dat.s1$Male, dat.s1$A1*dat.s1$Education, dat.s1$A1*dat.s1$HDDays, 
                                       dat.s1$A1*dat.s1$CenteredAge*ifelse(dat.s1$time==2, 1, 0), dat.s1$A1*dat.s1$CenteredAge*ifelse(dat.s1$time==3, 1, 0), 
                                       dat.s1$A1*dat.s1$Male*ifelse(dat.s1$time==2, 1, 0), dat.s1$A1*dat.s1$Male*ifelse(dat.s1$time==3, 1, 0), 
                                       dat.s1$A1*dat.s1$Education*ifelse(dat.s1$time==2, 1, 0), dat.s1$A1*dat.s1$Education*ifelse(dat.s1$time==3, 1, 0), 
                                       dat.s1$A1*dat.s1$HDDays*ifelse(dat.s1$time==2, 1, 0), dat.s1$A1*dat.s1$HDDays*ifelse(dat.s1$time==3, 1, 0))) %*% coef(mod.s1)[c(8, 17:30)]
    dat.s1.opt$Y <- dat.s1$Y + 2 * ifelse(dat.s1.opt$A1 == dat.s1$A1, 0, 1) * abs(int1)
  }
  
  dat.opt <- cbind.data.frame(dat.s1.opt[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1", "R", "A2", "time", "Y")])
  # convert dat.opt from long format to wide format
  dat.opt <- reshape(dat.opt, v.names = "Y", idvar = "SubjectID", timevar = "time", direction = "wide", sep = "")
  
  return(list("mod.s1" = mod.s1, "mod.s2" = mod.s2, "dat.opt" = dat.opt))
}

output.mqgee <- analyzeData(dat, method = "MQGEE")
output.sqgee <- analyzeData(dat, method = "SQGEE")

## estimate the (optimal or specific) DTRs as a by-product of Q-learning

estTraj <- function(dat, a1, a2, corstr){
  n <- nrow(dat)
  dat.s2 <- cbind.data.frame(do.call(rbind, replicate(2, dat[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1", "Y1", "R", "A2")], simplify = FALSE)), 
                             "time" = rep(c(2, 3), each = n), "Y" = c(dat$Y2, dat$Y3))
  dat.s2 <- dat.s2[order(dat.s2$SubjectID), ]
  rand.sample <- dat.s2$R==0
  mod.s2 <- geeglm(Y ~ factor(time) + CenteredAge + Male + Education + HDDays + A1 + Y1 + factor(time)*(CenteredAge + Male + Education + HDDays + A1 + Y1) + 
                     A2 + A2*(factor(time) + CenteredAge + Male + Education + HDDays + A1 + Y1 + factor(time)*(CenteredAge + Male + Education + HDDays + A1 + Y1)),
                   data = dat.s2[rand.sample, ], id = SubjectID, corstr = corstr[1])
  dat.s2.est <- dat.s2
  dat.s2.est$A2[rand.sample] <- a2
  dat.s2.est$Y[rand.sample] <- predict(mod.s2, newdata = dat.s2.est[rand.sample,])
  
  # stage 1 estimation assuming stage 1 treatment a1
  dat.s2.est <- dat.s2.est[order(dat.s2.est$time), ]
  dat.s1 <- cbind.data.frame(do.call(rbind, replicate(3, dat[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1", "Y1", "R")], simplify = FALSE)), 
                             "A2" = rep(dat.s2.est$A2[1:n], 3),
                             "time" = rep(c(1, 2, 3), each = n), 
                             "Y" = c(dat$Y1, dat.s2.est$Y))
  dat.s1 <- dat.s1[order(dat.s1$SubjectID), ]
  mod.s1 <- geeglm(Y ~ factor(time) + CenteredAge + Male + Education + HDDays + factor(time)*(CenteredAge + Male + Education + HDDays) + 
                     A1 + A1 * (factor(time) + CenteredAge + Male + Education + HDDays + factor(time)*(CenteredAge + Male + Education + HDDays)), 
                   data = dat.s1, id = SubjectID, corstr = corstr[2])
  dat.s1.est <- dat.s1
  dat.s1.est$A1 <- a1
  dat.s1.est$Y <- predict(mod.s1, newdata = dat.s1.est)
  dat.est <- cbind.data.frame(dat.s1.est[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1", "R", "A2", "time", "Y")])
  # convert dat.opt from long format to wide format
  dat.est <- reshape(dat.est, v.names = "Y", idvar = "SubjectID", timevar = "time", direction = "wide", sep = "")
  return(list("mod.s1" = mod.s1, "mod.s2" = mod.s2, "dat.est" = dat.est))
}

## summary and plots of analysis results

# plot the trajectory of DTRs

sumTraj <- function(dat.opt, corstr = c("exchangeable", "unstructured")){
  dat.est1 <- estTraj(dat, a1 = -1, a2 = 1, corstr)$dat.est
  dat.est2 <- estTraj(dat, a1 = 1, a2 = 1, corstr)$dat.est
  dat.est3 <- estTraj(dat, a1 = -1, a2 = -1, corstr)$dat.est
  dat.est4 <- estTraj(dat, a1 = 1, a2 = -1, corstr)$dat.est
  y.opt <- c(mean(dat.opt$Y1), mean(dat.opt$Y2), mean(dat.opt$Y3))
  y.est1 <- c(mean(dat.est1$Y1), mean(dat.est1$Y2), mean(dat.est1$Y3))
  y.est2 <- c(mean(dat.est2$Y1), mean(dat.est2$Y2), mean(dat.est2$Y3))
  y.est3 <- c(mean(dat.est3$Y1), mean(dat.est3$Y2), mean(dat.est3$Y3))
  y.est4 <- c(mean(dat.est4$Y1), mean(dat.est4$Y2), mean(dat.est4$Y3))
  pred <- cbind.data.frame("time" = rep(c(2, 3, 6), 5),
                           "y" = c(y.opt, y.est1, y.est2, y.est3, y.est4), 
                           "DTR" = c(rep("Optimal DTR", 3), rep("(c) Later choice", 3), rep("(a) Choice throughout", 3), rep("(d) No-choice", 3), rep("(b) Initial-choice", 3)))
  return(pred)
}

pred.mqgee <- sumTraj(output.mqgee$dat.opt)
pred.sqgee <- sumTraj(output.sqgee$dat.opt)

png(filename="OptAndEstENGAGE_mqgee.png", width = 15, height = 10, units = "cm", res = 200)
pred.mqgee %>%
  mutate(DTR = fct_relevel(DTR, "Optimal DTR", "(a) Choice throughout", "(b) Initial-choice", "(c) Later choice", "(d) No-choice")) %>%
ggplot(aes(x = time, y = y, shape = DTR, linetype = DTR)) + 
  geom_point() + geom_line(aes(group = DTR)) + 
  labs(title= "", x = "Time (months)", y = "Treatment readiness") + 
  scale_x_continuous(limits = c(1, 6), breaks = seq(1, 6, by = 1)) +
  theme_bw()
dev.off()

# percentage of patients assigned to a specific decision rule

table(output.mqgee$dat.opt$A1)/nrow(dat)
# -1     1 
# 0.844 0.156 
table(output.mqgee$dat.opt$A2)/nrow(dat)
# -1     0     1 
# 0.068 0.700 0.232 

# boxplot of estimated CATEs (stage k, time j)

calcCATE <- function(output){
  dat1.temp <- cbind.data.frame(do.call(rbind, replicate(3, dat[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1")], simplify = FALSE)), 
                                "time" = rep(c(1, 2, 3), each = nrow(dat)))
  dat1.temp <- dat1.temp[order(dat1.temp$SubjectID), ]
  dat1a.temp <- dat1.temp
  dat1a.temp$A1 <- -1
  dat1b.temp <- dat1.temp
  dat1b.temp$A1 <- 1
  delta1.hat <- matrix(predict(output$mod.s1, newdata = dat1a.temp), ncol = 3, byrow = TRUE) - matrix(predict(output$mod.s1, newdata = dat1b.temp), ncol = 3, byrow = TRUE)
  dat2.temp <- cbind.data.frame(do.call(rbind, replicate(2, dat[ , c("SubjectID", "CenteredAge", "Male", "Education", "HDDays", "A1", "Y1", "R", "A2")], simplify = FALSE)), 
                                "time" = rep(c(2, 3), each = nrow(dat)))
  dat2.temp <- dat2.temp[order(dat2.temp$SubjectID), ]
  dat2a.temp <- dat2.temp
  dat2a.temp$A2 <- -1
  dat2b.temp <- dat2.temp
  dat2b.temp$A2 <- 1
  delta2.hat <- matrix(predict(output$mod.s2, newdata = dat2a.temp), ncol = 2, byrow = TRUE) - matrix(predict(output$mod.s2, newdata = dat2b.temp), ncol = 2, byrow = TRUE)
  CATE <- cbind.data.frame(delta1.hat, delta2.hat)
  colnames(CATE) <- c("s1t1", "s1t2", "s1t3", "s2t2", "s2t3")
  return(CATE)
}

CATE <- calcCATE(output.mqgee)

library(reshape2)
CATE.temp <- melt(CATE)
addline_format <- function(x,...){
  gsub('\\s','\n', x)
}
png(filename="CATEDistENGAGE_mqgee.png", width = 12, height = 8, units = "cm", res = 200)
ggplot(data = CATE.temp, aes(x = variable, y = value)) + 
  geom_boxplot() + 
  labs(title= "", x = "", y = "Heterogeneous Causal Effect") + 
  scale_x_discrete(labels = c("stage 1 Q-fn\n2 months", "stage 1 Q-fn\n3 months", "stage 1 Q-fn\n6 months",
                              "stage 2 Q-fn\n3 months", "stage 2 Q-fn\n6 months")) + 
  theme_bw()
dev.off()

## model diagnostics
QIC(output.sqgee$mod.s1)
QIC(output.mqgee$mod.s1)

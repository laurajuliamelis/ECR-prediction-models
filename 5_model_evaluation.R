
# ============================================================================ #
# ===================        EOR PREDICTIVE MODELS           ================= #
# ===================          Model evaluation              ================= #
# ============================================================================ #

rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary packages
library(openxlsx)
library(skimr)
library(knitr)
library(mice)
library(dplyr)
library(gtsummary)
library(VIM)
library(ggplot2)
library(gridExtra)
library(glmnet)
library(caret)
library(pROC)
library(Epi)
library(OptimalCutpoints)
library(RColorBrewer)
library(caret)
library(predtools)
library(rms)
library(CalibrationCurves)
library(boot)
library(papeR)
library(gbm)
library(naniar)
library(stringr)
library(egg)
library(ggsignif)
library(ggsci)
library(fastshap)
library(shapviz)
library(reshape2)
library(viridis)
library(patchwork)
library(ggpubr)

# Load myfunctions
source("myfunctions.R")

# --------------------------------- #
# 1. Area under the ROC Curve -----
# --------------------------------- #

### LOGISTIC REGRESSION ----
rms.procs <- list(roc.clin, roc.gen, roc.clingen) 
names(rms.procs) <- c(paste0("Clinical (", round(auc.clin[12,1], 2), "-",  round(auc.clin[12,5],2),")"),
                      paste0("Genetic (", round(auc.gen[12,1], 2), "-",  round(auc.gen[12,5],2), ")"),  
                      paste0("Clin+Gen (", round(auc.clingen[12,1], 2), "-",  round(auc.clingen[12,5],2), ")"))
p1 <- ggroc(rms.procs,linewidth = 0.8) + 
  theme_article() + 
  theme(legend.position = c(0.8, 0.13), 
        legend.key.height = unit(0.4,"cm"),
        legend.title = element_text(face="bold", size=10)) + 
  labs(colour="Model (AUC-AUCc)", x="Specificity", y="Sensitivity",
       title = "Logistic Regression", tag="A") +
  scale_color_brewer(palette="Paired") +
  geom_abline(intercept = 1, slope = 1, linetype=2) 


### NAIVE BAYES ----
naive.procs <- list(naive.roc.clin, naive.roc.gen, naive.roc.clingen) 
names(naive.procs) <-c(paste0("Clinical (", round(naive.ci.clin, 2)[2], "-",  round(naive.clin.auc.c,2),")"), 
                       paste0("Genetic (", round(naive.ci.gen, 2)[2], "-",  round(naive.gen.auc.c,2), ")"), 
                       paste0("Clin+Gen  (", round(naive.ci.clingen, 2)[2], "-",  round(naive.clingen.auc.c,2), ")")) 

p2 <- ggroc(naive.procs, linewidth = 0.8) + 
  theme_article() + 
  theme(legend.position = c(0.8, 0.13), 
        legend.key.height = unit(0.4,"cm"),
        legend.title = element_text(face="bold", size=10)) + 
  labs(colour="Model (AUC-AUCc)", x="Specificity", y="Sensitivity",
       title = "Naive Bayes Classifier", tag="B") +
  scale_color_brewer(palette="Paired") +
  geom_abline(intercept = 1, slope = 1, linetype=2)


### GBM ----
gbm.procs <- list(gbm.roc.clin, gbm.roc.gen, gbm.roc.clingen) 
names(gbm.procs) <-c(paste0("Clinical (", round(gbm.ci.clin, 2)[2], "-",  round(gbm.clin.auc.c,2),")"), 
                     paste0("Genetic (", round(gbm.ci.gen, 2)[2], "-",  round(gbm.gen.auc.c,2), ")"),
                     paste0("Clin+Gen  (", round(gbm.ci.clingen, 2)[2], "-",  round(gbm.clingen.auc.c,2), ")")) 

p3 <- ggroc(gbm.procs, linewidth = 0.8) + 
  theme_article() + 
  theme(legend.position = c(0.8, 0.13), 
        legend.key.height = unit(0.4,"cm"),
        legend.title = element_text(face="bold", size=10)) + 
  labs(colour="Model (AUC-AUCc)", x="Specificity", y="Sensitivity",
       title = "Gradient Boosting Machines", tag="C") +
  scale_color_brewer(palette="Paired") +
  geom_abline(intercept = 1, slope = 1, linetype=2)

### SVM ----
svm.procs <- list(svm.roc.clin, svm.roc.gen, svm.roc.clingen) 
names(svm.procs) <-c(paste0("Clinical (", round(svm.ci.clin, 2)[2], "-",  round(svm.clin.auc.c,2),")"), 
                     paste0("Genetic (", round(svm.ci.gen, 2)[2], "-",  round(svm.gen.auc.c,2), ")"), 
                     paste0("Clin+Gen  (", round(svm.ci.clingen, 2)[2], "-",  round(svm.clingen.auc.c,2), ")")) 

p4 <- ggroc(svm.procs, linewidth = 0.8) + 
  theme_article() + 
  theme(legend.position = c(0.8, 0.13), 
        legend.key.height = unit(0.4,"cm"),
        legend.title = element_text(face="bold", size=10)) + 
  labs(colour="Model (AUC-AUCc)", x="Specificity", y="Sensitivity",
       title = "Support Vector Machines", tag="D") +
  scale_color_brewer(palette="Paired") +
  geom_abline(intercept = 1, slope = 1, linetype=2)


grid.arrange(p1,p2,p3,p4, ncol=2)



# ------------------------------- #
# 2. Test and 95%CI for AUC -----
# ------------------------------- #

# We will compute the 95\% confidence interval for the AUC (the uncorrected for optimism) with 2000 stratified bootstrap replicates (stratified means that each bootstrap resample contains exactly the same number of cases and controls than the original sample).
# We will also use the bootstrap method to compare each pair of ROC curves. Particularly, we'll test for a difference in the ROC curves and therefore the null hypothesis establishes that the difference is 0 (i.e. AUC1=AUC2, the AUCof one doesn't fall inside the CI95\% of the other) so if the pvalue of the test is lower than 0.05 we have sufficient statistical evidences to reject that hypothesis. 


## 2.1. Results for models ----
data.auc <- data.frame("Method"= rep(c("Logistic", "Naive", "GBM", "SVM"), each=3),
                       "Model" = rep(c("Clinical", "Genetic", "Clin+Gen"), times=4),
                       "AUC"= c(round(logisitc.ci.clin[[2]], 3),
                                round(logisitc.ci.gen[[2]], 3),
                                round(logisitc.ci.clingen[[2]], 3),
                                round(naive.ci.clin[[2]], 3),
                                round(naive.ci.gen[[2]], 3),
                                round(naive.ci.clingen[[2]], 3),
                                round(gbm.ci.clin[[2]], 3),
                                round(gbm.ci.gen[[2]], 3),
                                round(gbm.ci.clingen[[2]], 3),
                                round(svm.ci.clin[[2]], 3),
                                round(svm.ci.gen[[2]], 3),
                                round(svm.ci.clingen[[2]], 3)),
                       "low.AUC"= c(round(logisitc.ci.clin[[1]], 3),
                                    round(logisitc.ci.gen[[1]], 3),
                                    round(logisitc.ci.clingen[[1]], 3),
                                    round(naive.ci.clin[[1]], 3),
                                    round(naive.ci.gen[[1]], 3),
                                    round(naive.ci.clingen[[1]], 3),
                                    round(gbm.ci.clin[[1]], 3),
                                    round(gbm.ci.gen[[1]], 3),
                                    round(gbm.ci.clingen[[1]], 3),
                                    round(svm.ci.clin[[1]], 3),
                                    round(svm.ci.gen[[1]], 3),
                                    round(svm.ci.clingen[[1]], 3)),
                       "upp.AUC"= c(round(logisitc.ci.clin[[3]], 3),
                                    round(logisitc.ci.gen[[3]], 3),
                                    round(logisitc.ci.clingen[[3]], 3),
                                    round(naive.ci.clin[[3]], 3),
                                    round(naive.ci.gen[[3]], 3),
                                    round(naive.ci.clingen[[3]], 3),
                                    round(gbm.ci.clin[[3]], 3),
                                    round(gbm.ci.gen[[3]], 3),
                                    round(gbm.ci.clingen[[3]], 3),
                                    round(svm.ci.clin[[3]], 3),
                                    round(svm.ci.gen[[3]], 3),
                                    round(svm.ci.clingen[[3]], 3)))
data.auc$Method <- factor(data.auc$Method, levels=c("Logistic", "Naive", "GBM", "SVM"))
data.auc$Model <- factor(data.auc$Model, levels=c("Clinical", "Genetic", "Clin+Gen"))

set.seed(123)
p.log.clin.gen <- roc.test(roc.clin, roc.gen, method="b")
p.log.clin.clingen <- roc.test(roc.clin, roc.clingen, method="b")
p.log.gen.clingen <- roc.test(roc.gen, roc.clingen, method="b")

p.naive.clin.gen <- roc.test(naive.roc.clin, naive.roc.gen, method="b")
p.naive.clin.clingen <- roc.test(naive.roc.clin, naive.roc.clingen, method="b")
p.naive.gen.clingen <- roc.test(naive.roc.gen, naive.roc.clingen, method="b")

p.gbm.clin.gen <- roc.test(gbm.roc.clin, gbm.roc.gen, method="b")
p.gbm.clin.clingen <- roc.test(gbm.roc.clin, gbm.roc.clingen, method="b")
p.gbm.gen.clingen <- roc.test(gbm.roc.gen, gbm.roc.clingen, method="b")

p.svm.clin.gen <- roc.test(svm.roc.clin, svm.roc.gen, method="b")
p.svm.clin.clingen <- roc.test(svm.roc.clin, svm.roc.clingen, method="b")
p.svm.gen.clingen <- roc.test(svm.roc.gen, svm.roc.clingen, method="b")


tests.AUC <- data.frame("Method"= rep(c("Logistic", "Naive", "GBM", "SVM"), each=3),
                        "Model1" = rep(c("Clinical", "Clinical", "Genetic"), times=4),
                        "Model2" = rep(c("Genetic", "Clin+Gen", "Clin+Gen"), times=4),
                        "Difference"=c(round(logisitc.ci.clin[[2]]-logisitc.ci.gen[[2]], 3),
                                       round(logisitc.ci.clin[[2]]-logisitc.ci.clingen[[2]], 3),
                                       round(logisitc.ci.gen[[2]]-logisitc.ci.clingen[[2]], 3),
                                       round(naive.ci.clin[[2]]-naive.ci.gen[[2]], 3),
                                       round(naive.ci.clin[[2]]-naive.ci.clingen[[2]], 3),
                                       round(naive.ci.gen[[2]]-naive.ci.clingen[[2]], 3),
                                       round(gbm.ci.clin[[2]]-gbm.ci.gen[[2]], 3),
                                       round(gbm.ci.clin[[2]]-gbm.ci.clingen[[2]], 3),
                                       round(gbm.ci.gen[[2]]-gbm.ci.clingen[[2]], 3),
                                       round(svm.ci.clin[[2]]-svm.ci.gen[[2]], 3),
                                       round(svm.ci.clin[[2]]-svm.ci.clingen[[2]], 3),
                                       round(svm.ci.gen[[2]]-svm.ci.clingen[[2]], 3)),
                        "z"=c(round(p.log.clin.gen$statistic, 2),
                              round(p.log.clin.clingen$statistic, 2),
                              round(p.log.gen.clingen$statistic, 2),
                              round(p.naive.clin.gen$statistic, 2),
                              round(p.naive.clin.clingen$statistic, 2),
                              round(p.naive.gen.clingen$statistic, 2),
                              round(p.gbm.clin.gen$statistic, 2),
                              round(p.gbm.clin.clingen$statistic, 2),
                              round(p.gbm.gen.clingen$statistic, 2),
                              round(p.svm.clin.gen$statistic, 2),
                              round(p.svm.clin.clingen$statistic, 2),
                              round(p.svm.gen.clingen$statistic, 2)),
                        "p"=c(round(p.log.clin.gen$p.value, 2),
                              round(p.log.clin.clingen$p.value, 2),
                              round(p.log.gen.clingen$p.value, 2),
                              round(p.naive.clin.gen$p.value, 2),
                              round(p.naive.clin.clingen$p.value, 2),
                              round(p.naive.gen.clingen$p.value, 2),
                              round(p.gbm.clin.gen$p.value, 2),
                              round(p.gbm.clin.clingen$p.value, 2),
                              round(p.gbm.gen.clingen$p.value, 2),
                              round(p.svm.clin.gen$p.value, 2),
                              round(p.svm.clin.clingen$p.value, 2),
                              round(p.svm.gen.clingen$p.value, 2)))

annotations <- data.frame(
  Method = factor(tests.AUC$Method, levels = c("Logistic", "Naive", "GBM", "SVM")),
  y_position = rep(c(1.02, 1.065, 1.02), times = 4),
  xmin = rep(c(1, 1, 2), times = 4), 
  xmax = rep(c(1.97, 3, 3), times = 4),
  annotations = c(paste0("p=", tests.AUC$p[1]), paste0("p=", tests.AUC$p[2]), paste0("p=", tests.AUC$p[3]),
                  paste0("p=", tests.AUC$p[4]), paste0("p=", tests.AUC$p[5]), paste0("p=", tests.AUC$p[6]), 
                  rep("p<0.001",3), paste0("p=", tests.AUC$p[10]), paste0("p=", tests.AUC$p[11]), "p<0.001")
)

ggplot(data.auc) +
  geom_bar(aes(x=Model, y= AUC, fill=Model), stat="identity",alpha=0.7, width=0.6)+
  geom_errorbar( aes(x=Model, ymin=low.AUC, ymax=upp.AUC), width=0.1, colour="grey30", size=0.8) +
  scale_fill_brewer(palette = "Paired")+
  facet_wrap(~Method, nrow=1) + theme_article() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1.065), breaks = seq(0, 1, by = 0.25)) + 
  theme(strip.text.x = element_text(size = 8,face="bold")) +
  geom_signif(
    data = annotations,
    aes(xmin = xmin, xmax = xmax, annotations = annotations, 
        y_position = y_position,group=c(rep(1:3,4))),
    manual = TRUE,
    tip_length = 0.02,
    textsize = 2.5
  )


# for models
data.auc.all <- rbind(data.auc, data.auc.ensemble)

## LOGISTIC
data.aux <- subset(data.auc.all, Method == "Logistic")
comparisons <- data.frame(t(data.frame(combn(data.aux$Model, m=2, simplify = F))),row.names = NULL)
names(comparisons) <- c("Model1", "Model2")
comparisons$difference <- comparisons$pval <- comparisons$statistic <- rep(NA, nrow(comparisons))

set.seed(123)
for(i in 1:nrow(comparisons)){
  roc1 <- ifelse(comparisons$Model1[i] == "Clinical", "roc.clin", 
                 ifelse(comparisons$Model1[i] == "Genetic", "roc.gen", "roc.clingen"))
  roc2 <- ifelse(comparisons$Model2[i] == "Clinical", "roc.clin", 
                 ifelse(comparisons$Model2[i] == "Genetic", "roc.gen", "roc.clingen"))
  aux <- roc.test(get(roc1), get(roc2), method="b")
  comparisons$statistic[i] <- round(aux$statistic,2)
  comparisons$pval[i] <- round(aux$p.value,2)
  
  id1 <- which(data.aux$Model == comparisons$Model1[i]) 
  id2 <- which(data.aux$Model == comparisons$Model2[i]) 
  comparisons$difference[i] <- data.aux[id1,3] - data.aux[id2,3]
}

comparisons.log <- comparisons

## NAIVE BAYES
data.aux <- subset(data.auc.all, Method == "Naive")
comparisons <- data.frame(t(data.frame(combn(data.aux$Model, m=2, simplify = F))),row.names = NULL)
names(comparisons) <- c("Model1", "Model2")
comparisons$difference <- comparisons$pval <- comparisons$statistic <- rep(NA, nrow(comparisons))

set.seed(123)
for(i in 1:nrow(comparisons)){
  roc1 <- ifelse(comparisons$Model1[i] == "Clinical", "naive.roc.clin", 
                 ifelse(comparisons$Model1[i] == "Genetic", "naive.roc.gen", "naive.roc.clingen"))
  roc2 <- ifelse(comparisons$Model2[i] == "Clinical", "naive.roc.clin", 
                 ifelse(comparisons$Model2[i] == "Genetic", "naive.roc.gen", "naive.roc.clingen"))
  aux <- roc.test(get(roc1), get(roc2), method="b")
  comparisons$statistic[i] <- round(aux$statistic,2)
  comparisons$pval[i] <- round(aux$p.value,2)
  
  id1 <- which(data.aux$Model == comparisons$Model1[i]) 
  id2 <- which(data.aux$Model == comparisons$Model2[i]) 
  comparisons$difference[i] <- data.aux[id1,3] - data.aux[id2,3]
}

comparisons.naive <- comparisons


## GBM
data.aux <- subset(data.auc.all, Method == "GBM")
comparisons <- data.frame(t(data.frame(combn(data.aux$Model, m=2, simplify = F))),row.names = NULL)
names(comparisons) <- c("Model1", "Model2")
comparisons$difference <- comparisons$pval <- comparisons$statistic <- rep(NA, nrow(comparisons))

set.seed(123)
for(i in 1:nrow(comparisons)){
  roc1 <- ifelse(comparisons$Model1[i] == "Clinical", "gbm.roc.clin", 
                 ifelse(comparisons$Model1[i] == "Genetic", "gbm.roc.gen", "gbm.roc.clingen"))
  roc2 <- ifelse(comparisons$Model2[i] == "Clinical", "gbm.roc.clin", 
                 ifelse(comparisons$Model2[i] == "Genetic", "gbm.roc.gen", "gbm.roc.clingen"))
  aux <- roc.test(get(roc1), get(roc2), method="b")
  comparisons$statistic[i] <- round(aux$statistic,2)
  comparisons$pval[i] <- round(aux$p.value,2)
  
  id1 <- which(data.aux$Model == comparisons$Model1[i]) 
  id2 <- which(data.aux$Model == comparisons$Model2[i]) 
  comparisons$difference[i] <- data.aux[id1,3] - data.aux[id2,3]
}

comparisons.gbm <- comparisons

## SVM
data.aux <- subset(data.auc.all, Method == "SVM")
comparisons <- data.frame(t(data.frame(combn(data.aux$Model, m=2, simplify = F))),row.names = NULL)
names(comparisons) <- c("Model1", "Model2")
comparisons$difference <- comparisons$pval <- comparisons$statistic <- rep(NA, nrow(comparisons))

set.seed(123)
for(i in 1:nrow(comparisons)){
  roc1 <- ifelse(comparisons$Model1[i] == "Clinical", "svm.roc.clin", 
                 ifelse(comparisons$Model1[i] == "Genetic", "svm.roc.gen", "svm.roc.clingen"))
  roc2 <- ifelse(comparisons$Model2[i] == "Clinical", "svm.roc.clin", 
                 ifelse(comparisons$Model2[i] == "Genetic", "svm.roc.gen", "svm.roc.clingen"))
  aux <- roc.test(get(roc1), get(roc2), method="b")
  comparisons$statistic[i] <- round(aux$statistic,2)
  comparisons$pval[i] <- round(aux$p.value,2)
  
  id1 <- which(data.aux$Model == comparisons$Model1[i]) 
  id2 <- which(data.aux$Model == comparisons$Model2[i]) 
  comparisons$difference[i] <- data.aux[id1,3] - data.aux[id2,3]
}

comparisons.svm <- comparisons


comparisons <- rbind(comparisons.log, comparisons.naive, comparisons.svm, comparisons.gbm)
comparisons$Method <- rep(c("Logistic", "Naive Bayes", "SVM", "GBM", ), each=3)
comparisons <- comparisons[,c(6, 1:4)]
kable(comparisons, "pipe", 
      caption = "AUC (95%CI) and comparison of predictive performances for ALL models with each method.")


## 2.2. Results for methods ----
set.seed(123)
p.clin.log.naive <- roc.test(roc.clin, naive.roc.clin, method="b")
p.clin.log.gbm <- roc.test(roc.clin, gbm.roc.clin, method="b")
p.clin.log.svm <- roc.test(roc.clin, svm.roc.clin, method="b")
p.clin.naive.gbm <- roc.test(naive.roc.clin, gbm.roc.clin, method="b")
p.clin.naive.svm <- roc.test(naive.roc.clin, svm.roc.clin, method="b")
p.clin.gbm.svm <- roc.test(gbm.roc.clin, svm.roc.clin, method="b")

p.gen.log.naive <- roc.test(roc.gen, naive.roc.gen, method="b")
p.gen.log.gbm <- roc.test(roc.gen, gbm.roc.gen, method="b")
p.gen.log.svm <- roc.test(roc.gen, svm.roc.gen, method="b")
p.gen.naive.gbm <- roc.test(naive.roc.gen, gbm.roc.gen, method="b")
p.gen.naive.svm <- roc.test(naive.roc.gen, svm.roc.gen, method="b")
p.gen.gbm.svm <- roc.test(gbm.roc.gen, svm.roc.gen, method="b")

p.clingen.log.naive <- roc.test(roc.clingen, naive.roc.clingen, method="b")
p.clingen.log.gbm <- roc.test(roc.clingen, gbm.roc.clingen, method="b")
p.clingen.log.svm <- roc.test(roc.clingen, svm.roc.clingen, method="b")
p.clingen.naive.gbm <- roc.test(naive.roc.clingen, gbm.roc.clingen, method="b")
p.clingen.naive.svm <- roc.test(naive.roc.clingen, svm.roc.clingen, method="b")
p.clingen.gbm.svm <- roc.test(gbm.roc.clingen, svm.roc.clingen, method="b")

tests.AUC2 <- data.frame("Model"= rep(c("Clinical", "Genetic", "Clin+Gen"), each=6),
                         "Method1" = rep(c("Logistic", "Logistic", "Logistic", "Naive", "Naive", "GBM"), 3),
                         "Method2" = rep(c("Naive", "GBM", "SVM", "GBM", "SVM", "SVM"), 3),
                         "Difference"=c(round(logisitc.ci.clin[[2]]-naive.ci.clin[[2]], 3),
                                        round(logisitc.ci.clin[[2]]-gbm.ci.clin[[2]], 3),
                                        round(logisitc.ci.clin[[2]]-svm.ci.clin[[2]], 3),
                                        round(naive.ci.clin[[2]]-gbm.ci.clin[[2]], 3),
                                        round(naive.ci.clin[[2]]-svm.ci.clin[[2]], 3),
                                        round(gbm.ci.clin[[2]]-svm.ci.clin[[2]], 3),
                                        round(logisitc.ci.gen[[2]]-naive.ci.gen[[2]], 3),
                                        round(logisitc.ci.gen[[2]]-gbm.ci.gen[[2]], 3),
                                        round(logisitc.ci.gen[[2]]-svm.ci.gen[[2]], 3),
                                        round(naive.ci.gen[[2]]-gbm.ci.gen[[2]], 3),
                                        round(naive.ci.gen[[2]]-svm.ci.gen[[2]], 3),
                                        round(gbm.ci.gen[[2]]-svm.ci.gen[[2]], 3),
                                        round(logisitc.ci.clingen[[2]]-naive.ci.clingen[[2]], 3),
                                        round(logisitc.ci.clingen[[2]]-gbm.ci.clingen[[2]], 3),
                                        round(logisitc.ci.clingen[[2]]-svm.ci.clingen[[2]], 3),
                                        round(naive.ci.clingen[[2]]-gbm.ci.clingen[[2]], 3),
                                        round(naive.ci.clingen[[2]]-svm.ci.clingen[[2]], 3),
                                        round(gbm.ci.clingen[[2]]-svm.ci.clingen[[2]], 3)),
                         "z"=c(round(p.clin.log.naive$statistic, 2),
                               round(p.clin.log.gbm$statistic, 2),
                               round(p.clin.log.svm$statistic, 2),
                               round(p.clin.naive.gbm$statistic, 2),
                               round(p.clin.naive.svm$statistic, 2),
                               round(p.clin.gbm.svm$statistic, 2),
                               round(p.gen.log.naive$statistic, 2),
                               round(p.gen.log.gbm$statistic, 2),
                               round(p.gen.log.svm$statistic, 2),
                               round(p.gen.naive.gbm$statistic, 2),
                               round(p.gen.naive.svm$statistic, 2),
                               round(p.gen.gbm.svm$statistic, 2),
                               round(p.clingen.log.naive$statistic, 2),
                               round(p.clingen.log.gbm$statistic, 2),
                               round(p.clingen.log.svm$statistic, 2),
                               round(p.clingen.naive.gbm$statistic, 2),
                               round(p.clingen.naive.svm$statistic, 2),
                               round(p.clingen.gbm.svm$statistic, 2)),
                         "p"=c(round(p.clin.log.naive$p.value, 2),
                               round(p.clin.log.gbm$p.value, 2),
                               round(p.clin.log.svm$p.value, 2),
                               round(p.clin.naive.gbm$p.value, 2),
                               round(p.clin.naive.svm$p.value, 2),
                               round(p.clin.gbm.svm$p.value, 2),
                               round(p.gen.log.naive$p.value, 2),
                               round(p.gen.log.gbm$p.value, 2),
                               round(p.gen.log.svm$p.value, 2),
                               round(p.gen.naive.gbm$p.value, 2),
                               round(p.gen.naive.svm$p.value, 2),
                               round(p.gen.gbm.svm$p.value, 2),
                               round(p.clingen.log.naive$p.value, 2),
                               round(p.clingen.log.gbm$p.value, 2),
                               round(p.clingen.log.svm$p.value, 2),
                               round(p.clingen.naive.gbm$p.value, 2),
                               round(p.clingen.naive.svm$p.value, 2),
                               round(p.clingen.gbm.svm$p.value, 2)))


annotations <- data.frame(
  Model = factor(tests.AUC2$Model, levels = c("Clinical", "Genetic", "Clin+Gen")),
  y_position = rep(c(1, 1.04, 1.13, 1, 1.09, 1), times = 3),
  xmin = rep(c(1, 1, 1, 2.03, 2, 3.03), times = 3),
  xmax = rep(c(1.97, 3, 4, 2.97, 4, 4), times = 3),
  annotations = c(paste0("p=", tests.AUC2$p[1:7]), rep("p<0.001", times=5),
                  paste0("p=", tests.AUC2$p[13]), "p<0.001",
                  paste0("p=", tests.AUC2$p[15]),"p<0.001",
                  paste0("p=", tests.AUC2$p[17]),"p<0.001")
)


ggplot(data.auc) +
  geom_bar(aes(x=Method, y= AUC, fill=Method), stat="identity",alpha=0.7, width=0.7)+
  geom_errorbar( aes(x=Method, ymin=low.AUC, ymax=upp.AUC), width=0.1, colour="grey30", size=0.8) +
  scale_fill_brewer(palette = "Set2")+
  facet_wrap(~Model, nrow=1) + theme_article() + 
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1.13), breaks = seq(0, 1, by = 0.25)) + 
  theme(strip.text.x = element_text(size = 8,face="bold")) +
  geom_signif(
    data = annotations,
    aes(xmin = xmin, xmax = xmax, annotations = annotations, 
        y_position = y_position,group=rep(1:6,3)),
    manual = TRUE,
    tip_length = 0.02,
    textsize = 2.5
  )

data.auc.all <- rbind(data.auc, data.auc.ensemble)

## CLINICAL
clin.aux <- subset(data.auc.all, Model == "Clinical")
comparisons <- data.frame(t(data.frame(combn(clin.aux$Method, m=2, simplify = F))),row.names = NULL)
names(comparisons) <- c("Method1", "Method2")
comparisons$difference <- comparisons$pval <- comparisons$statistic <- rep(NA, nrow(comparisons))

set.seed(123)
for(i in 1:nrow(comparisons)){
  roc1 <- ifelse(comparisons$Method1[i] == "Logistic", "roc.clin", 
                 ifelse(comparisons$Method1[i] == "Naive", "naive.roc.clin", 
                        ifelse(comparisons$Method1[i] == "GBM", "gbm.roc.clin","svm.roc.clin")))
  roc2 <- ifelse(comparisons$Method2[i] == "Logistic", "roc.clin", 
                 ifelse(comparisons$Method2[i] == "Naive", "naive.roc.clin", 
                        ifelse(comparisons$Method2[i] == "GBM", "gbm.roc.clin", "svm.roc.clin")))
  aux <- roc.test(get(roc1), get(roc2), method="b")
  comparisons$statistic[i] <- round(aux$statistic,2)
  comparisons$pval[i] <- round(aux$p.value,2)
  
  id1 <- which(clin.aux$Method == comparisons$Method1[i]) 
  id2 <- which(clin.aux$Method == comparisons$Method2[i]) 
  comparisons$difference[i] <- clin.aux[id1,3] - clin.aux[id2,3]
}

comparisons.clin <- comparisons

## GENETIC
gen.aux <- subset(data.auc.all, Model == "Genetic")
comparisons <- data.frame(t(data.frame(combn(gen.aux$Method, m=2, simplify = F))),row.names = NULL)
names(comparisons) <- c("Method1", "Method2")
comparisons$difference <- comparisons$pval <- comparisons$statistic <- rep(NA, nrow(comparisons))

set.seed(123)
for(i in 1:nrow(comparisons)){
  roc1 <- ifelse(comparisons$Method1[i] == "Logistic", "roc.gen", 
                 ifelse(comparisons$Method1[i] == "Naive", "naive.roc.gen", 
                        ifelse(comparisons$Method1[i] == "GBM", "gbm.roc.gen", "svm.roc.gen")))
  roc2 <- ifelse(comparisons$Method2[i] == "Logistic", "roc.gen", 
                 ifelse(comparisons$Method2[i] == "Naive", "naive.roc.gen", 
                        ifelse(comparisons$Method2[i] == "GBM", "gbm.roc.gen","svm.roc.gen")))
  aux <- roc.test(get(roc1), get(roc2), method="b")
  comparisons$statistic[i] <- round(aux$statistic,2)
  comparisons$pval[i] <- round(aux$p.value,2)
  
  id1 <- which(gen.aux$Method == comparisons$Method1[i]) 
  id2 <- which(gen.aux$Method == comparisons$Method2[i]) 
  comparisons$difference[i] <- gen.aux[id1,3] - gen.aux[id2,3]
}

comparisons.gen <- comparisons


## CLIN + GEN
clingen.aux <- subset(data.auc.all, Model == "Clin+Gen")
comparisons <- data.frame(t(data.frame(combn(clingen.aux$Method, m=2, simplify = F))),row.names = NULL)
names(comparisons) <- c("Method1", "Method2")
comparisons$difference <- comparisons$pval <- comparisons$statistic <- rep(NA, nrow(comparisons))

set.seed(123)
for(i in 1:nrow(comparisons)){
  roc1 <- ifelse(comparisons$Method1[i] == "Logistic", "roc.clingen", 
                 ifelse(comparisons$Method1[i] == "Naive", "naive.roc.clingen", 
                        ifelse(comparisons$Method1[i] == "GBM", "gbm.roc.clingen", "svm.roc.clingen")))
  roc2 <- ifelse(comparisons$Method2[i] == "Logistic", "roc.clingen", 
                 ifelse(comparisons$Method2[i] == "Naive", "naive.roc.clingen", 
                        ifelse(comparisons$Method2[i] == "GBM", "gbm.roc.clingen", "svm.roc.clingen")))
  aux <- roc.test(get(roc1), get(roc2), method="b")
  comparisons$statistic[i] <- round(aux$statistic,2)
  comparisons$pval[i] <- round(aux$p.value,2)
  
  id1 <- which(clingen.aux$Method == comparisons$Method1[i]) 
  id2 <- which(clingen.aux$Method == comparisons$Method2[i]) 
  comparisons$difference[i] <- clingen.aux[id1,3] - clingen.aux[id2,3]
}

comparisons.clingen <- comparisons <- comparisons

comparisons <- rbind(comparisons.clin, comparisons.gen, comparisons.clingen)
comparisons$Model <- rep(c("Clinial", "Genetic", "Clin+Gen"), each=21)
comparisons <- comparisons[,c(5, 1:4)]
kable(comparisons, "pipe", 
      caption = "AUC (95%CI) and comparison of predictive performances for ALL methods with each model.")

# ------------------------------------------------ #
# 3. Table of predictive performance metrics -----
# ------------------------------------------------ #

## 3.1. Clinical Models ----
results <- as.data.frame(table.clin)
names(results)[1] <- "Method"
kable(results, "pipe", 
      caption = "Predictive Performance by Model and Dataset.")

## 3.2. Genetic Models ----
results <- as.data.frame(table.gen)
names(results)[1] <- "Method"
kable(results, "pipe", 
      caption = "Predictive Performance by Model and Dataset.")

## 3.3. Clin+Gen Models ----
results <- as.data.frame(table.clingen)
names(results)[1] <- "Method"
kable(results, "pipe", 
      caption = "Predictive Performance by Model and Dataset.")


# --------------------------- #
# 4. Confusion Matrices -----
# --------------------------- #

## 4.1. Clinical Models ----
#### LOGISTIC 
cm.rms.clin <- as.data.frame(metrics.clin$table)
cm.rms.clin <- cm.rms.clin %>%
  mutate(goodbad = ifelse(cm.rms.clin$Prediction == cm.rms.clin$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clin),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
rms1 <- ggplot(cm.rms.clin, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="Logistic")


#### NAIVE BAYES 
cm.naive.clin <- as.data.frame(naive.metrics.clin$table)
cm.naive.clin <- cm.naive.clin %>%
  mutate(goodbad = ifelse(cm.naive.clin$Prediction == cm.naive.clin$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clin),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
naive1 <- ggplot(cm.naive.clin, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="Naive Bayes")


#### GBM 
cm.gbm.clin <- as.data.frame(gbm.metrics.clin$table)
cm.gbm.clin <- cm.gbm.clin %>%
  mutate(goodbad = ifelse(cm.gbm.clin$Prediction == cm.gbm.clin$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clin),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
gbm1 <- ggplot(cm.gbm.clin, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="GBM")


#### SVM 
cm.svm.clin <- as.data.frame(svm.metrics.clin$table)
cm.svm.clin <- cm.svm.clin %>%
  mutate(goodbad = ifelse(cm.svm.clin$Prediction == cm.svm.clin$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clin),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
svm1 <- ggplot(cm.svm.clin, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="Predicted", title="SVM")
grid.arrange(rms1,naive1, gbm1, svm1, ncol=3)


## 4.2. Genetic Models ----
#### LOGISTIC 
cm.rms.gen <- as.data.frame(metrics.gen$table)
cm.rms.gen <- cm.rms.gen %>%
  mutate(goodbad = ifelse(cm.rms.gen$Prediction == cm.rms.gen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.gen),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
rms1 <- ggplot(cm.rms.gen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="Logistic")


#### NAIVE BAYES 
cm.naive.gen <- as.data.frame(naive.metrics.gen$table)
cm.naive.gen <- cm.naive.gen %>%
  mutate(goodbad = ifelse(cm.naive.gen$Prediction == cm.naive.gen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.gen),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
naive1 <- ggplot(cm.naive.gen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="Naive Bayes")


#### GBM 
cm.gbm.gen <- as.data.frame(gbm.metrics.gen$table)
cm.gbm.gen <- cm.gbm.gen %>%
  mutate(goodbad = ifelse(cm.gbm.gen$Prediction == cm.gbm.gen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.gen),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
gbm1 <- ggplot(cm.gbm.gen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="GBM")


#### SVM 
cm.svm.gen <- as.data.frame(svm.metrics.gen$table)
cm.svm.gen <- cm.svm.gen %>%
  mutate(goodbad = ifelse(cm.svm.gen$Prediction == cm.svm.gen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.gen),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
svm1 <- ggplot(cm.svm.gen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="Predicted", title="SVM")

grid.arrange(rms1,naive1, gbm1, svm1, ncol=3)


## 4.3. Clin+Gen Models ----
#### LOGISTIC 
cm.rms.clingen <- as.data.frame(metrics.clingen$table)
cm.rms.clingen <- cm.rms.clingen %>%
  mutate(goodbad = ifelse(cm.rms.clingen$Prediction == cm.rms.clingen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clingen.imp),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
rms1 <- ggplot(cm.rms.clingen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="Logistic")


#### NAIVE BAYES 
cm.naive.clingen <- as.data.frame(naive.metrics.clingen$table)
cm.naive.clingen <- cm.naive.clingen %>%
  mutate(goodbad = ifelse(cm.naive.clingen$Prediction == cm.naive.clingen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clingen.imp),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
naive1 <- ggplot(cm.naive.clingen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="Naive Bayes")


#### GBM 
cm.gbm.clingen <- as.data.frame(gbm.metrics.clingen$table)
cm.gbm.clingen <- cm.gbm.clingen %>%
  mutate(goodbad = ifelse(cm.gbm.clingen$Prediction == cm.gbm.clingen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clingen.imp),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
gbm1 <- ggplot(cm.gbm.clingen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="", title="GBM")


#### SVM 
cm.svm.clingen <- as.data.frame(svm.metrics.clingen$table)
cm.svm.clingen <- cm.svm.clingen %>%
  mutate(goodbad = ifelse(cm.svm.clingen$Prediction == cm.svm.clingen$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/nrow(data.clingen.imp),
         Freq2 = paste0(Freq, " (", round(prop*100,2), "%)"))
svm1 <- ggplot(cm.svm.clingen, aes(x=Prediction,y=Reference, fill = goodbad, alpha = prop))+
  geom_tile() +
  geom_text(aes(label = Freq2), vjust = .5, alpha = 1)+
  scale_fill_manual(values = c(good = "#73D17B", bad = "#EB8484"))  +
  theme_article() + theme(legend.position = "none")+
  scale_x_discrete(labels=c("Non EOR","EOR")) +
  scale_y_discrete(labels=c("Non EOR","EOR")) +
  labs(y="Actual", x="Predicted", title="SVM")


grid.arrange(rms1,naive1, gbm1, svm1, ncol=3)

# ---------------------------- #
# 5. Variable importance -----
# ---------------------------- #

## 5.1. Clinical model ----
naive.varimp <- varImp(naive.clin)$importance[, "Pos"]
gbm.varimp <- varImp(gbm.clin)$importance
svm.varimp <- varImp(svm.clin)$importance[, "Pos"]

varimp.clin <- round(cbind(naive.varimp, gbm.varimp, svm.varimp),2)
names(varimp.clin) <- c("Naive", "GBM", "SVM")
kable(varimp.clin, "pipe", 
      caption = "Variable importance for the Clinical Model per ML Algorithm.")

# Convert the dataframe to long format
varimp.clin$variable <- rownames(varimp.clin)
varimp.clin$variable[6] <- "PTV2M"
varimp_long <- pivot_longer(varimp.clin, cols = c(Naive, GBM, SVM), names_to = "Model", values_to = "Importance")

# Define a color palette from RColorBrewer
palette <- brewer.pal(n = 3, name = "Paired")

# Plot using ggplot2 with facets
p1 <- ggplot(varimp_long, aes(x = Importance, y = fct_reorder(variable, Importance), fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Clinical Model",
       x = "Importance",
       y = "Variable") +
  facet_wrap(~Model, scales = "free", ncol = 1) +
  scale_fill_manual(values = palette) +  # Use the defined color palette
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1), legend.position = "none")


## 5.2. Genetic model ----
naive.varimp <- varImp(naive.gen)$importance[, "Pos"]
gbm.varimp <- varImp(gbm.gen)$importance
svm.varimp <- varImp(svm.gen)$importance[, "Pos"]

varimp.gen <- round(cbind(naive.varimp, gbm.varimp, svm.varimp),2)
names(varimp.gen) <- c("Naive", "GBM", "SVM")
kable(varimp.gen, "pipe", 
      caption = "Variable importance for the Genetic Model per ML Algorithm.")

# Convert the dataframe to long format
varimp.gen$variable <- rownames(varimp.gen)
varimp_long <- pivot_longer(varimp.gen, cols = c(Naive, GBM, SVM), names_to = "Model", values_to = "Importance")

# Define a color palette from RColorBrewer
palette <- brewer.pal(n = 3, name = "Paired")

# Plot using ggplot2 with facets
p2 <- ggplot(varimp_long, aes(x = Importance, y = fct_reorder(variable, Importance), fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Genetic Model",
       x = "Importance",
       y = "") +
  facet_wrap(~Model, scales = "free", ncol = 1) +
  scale_fill_manual(values = palette) +  # Use the defined color palette
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1), legend.position = "none")



## 5.3. Clin+gen model ----
naive.varimp <- varImp(naive.clingen)$importance[, "Pos"]
gbm.varimp <- varImp(gbm.clingen)$importance
svm.varimp <- varImp(svm.clingen)$importance[, "Pos"]

varimp.clingen <- round(cbind(naive.varimp, gbm.varimp, svm.varimp),2)
names(varimp.clingen) <- c("Naive", "GBM", "SVM")
kable(varimp.clingen, "pipe", 
      caption = "Variable importance for the Clinical+Genetic Model per ML Algorithm.")

varimp.clingen$variable <- rownames(varimp.clingen)
varimp.clingen$variable[6] <- "PTV2M"
varimp_long <- pivot_longer(varimp.clingen, cols = c(Naive, GBM, SVM), names_to = "Model", values_to = "Importance")

# Define a color palette from RColorBrewer
palette <- brewer.pal(n = 3, name = "Paired")

# Plot using ggplot2 with facets
p3 <- ggplot(varimp_long, aes(x = Importance, y = fct_reorder(variable, Importance), fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title="Clin+Gen Model",
       x = "Importance",
       y = "") +
  facet_wrap(~Model, scales = "free", ncol = 1) +
  scale_fill_manual(values = palette) +  # Use the defined color palette
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1), legend.position = "none")


## Final plot ----
grid.arrange(p1,p2,p3, ncol=3)


# -------------------- #
# 6. SHAP Values -----
# -------------------- #

## 6.1. Compute Shanpley values ----
pfun <- function(object, newdata) {  # prediction wrapper
  unname(predict(object, type = "fitted", newdata = newdata))
}
X <- subset(data.clin.imp, select = c("DUP", "DTP", "EEAG_Total_VB", "Reserva_Cognitiva", "Insight", "respuestas_perseverativas_PTV2M"))  # features only

set.seed(2129)  # for reproducibility
ex.clin <- explain(rms.clin, X = X, pred_wrapper = pfun, newdata = X,nsim = 500, adjust = TRUE, shap_only = F)



## 6.2. Make plot ----
shap <- ex.clin$shapley_values
features <- ex.clin$feature_values

### Convert to lon format
shap_long <- melt(shap, varnames = c("observation", "variable"), value.name = "shap_value")
features_long <- melt(features, varnames = c("observation", "variable"), value.name = "feature_value")
shap_long$observation <- rep(1:nrow(features), times = ncol(features))
features_long$observation <- rep(1:nrow(features), times = ncol(features))

### Combine and rename variable (to english)
df <- left_join(shap_long, features_long, by = c("observation", "variable"))
df$variable <- ifelse(df$variable == "EEAG_Total_VB","Functionality",
                      ifelse(df$variable == "Reserva_Cognitiva", "Cognitive Reserve",
                             ifelse(df$variable == "respuestas_perseverativas_PTV2M", "Executive Function",
                                    ifelse(df$variable == "DUP","DUP",
                                           ifelse(df$variable == "DTP","DTP",  "Insight")))))

### Compute mean SHAP for each varible + percentage over total
mean_shap <- df %>%
  group_by(variable) %>%
  dplyr::summarize(mean_abs_shap = mean(abs(shap_value), na.rm = TRUE)) %>%
  mutate(variable_label = paste0(variable, " (", round(mean_abs_shap, 3), ")"))
mean_shap <-mean_shap[order(mean_shap$mean_abs_shap, decreasing = TRUE),]
mean_shap <- mean_shap %>%
  mutate(perc = (mean_abs_shap / sum(mean_abs_shap)) * 100)

df <- left_join(df, mean_shap, by = "variable")

### Normalize feature_value by percentiles (1%–99%)
df <- df %>%
  group_by(variable) %>%
  mutate(
    p1 = quantile(feature_value, 0.01, na.rm = TRUE),
    p99 = quantile(feature_value, 0.99, na.rm = TRUE),
    feature_value_scaled = pmin(pmax(feature_value, p1), p99),
    feature_value_scaled = (feature_value_scaled - p1) / (p99 - p1)
  ) %>%
  ungroup()

### Make sure the order of the variables is the same in both datasets
df$variable_label <- factor(df$variable_label, levels = mean_shap$variable_label)
mean_shap$variable_label <- factor(mean_shap$variable_label, levels = mean_shap$variable_label)


### Final plot
p1 <- ggplot(df, aes(x = shap_value, y =  reorder(variable_label, mean_abs_shap), color = feature_value_scaled)) +
  geom_jitter(height = 0.08, width = 0, alpha = 0.9, size = 1.3) +
  scale_color_viridis(option = "plasma", direction = -1, breaks = c(0, 1), labels = c("Low", "High")) +
  labs(
    x = "SHAP value (impact on model output)",
    y = NULL,
    color = "Patient score"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, linewidth = 1.2)

p2 <- ggplot(mean_shap, aes(x = reorder(variable_label, mean_abs_shap), y = perc)) +
  geom_col(colour= "grey30",fill = "grey80", width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 15) +
  labs(x = NULL, y = "% mean SHAP") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

final_plot <- p1 + p2 + plot_layout(widths = c(4, 1))
windows(width = 20, height = 12)
final_plot

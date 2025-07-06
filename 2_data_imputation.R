
# ============================================================================ #
# ===================        EOR PREDICTIVE MODELS           ================= #
# ===================           Data imputation              ================= #
# ============================================================================ #

rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary packages
library(knitr)
library(mice)
library(dplyr)
library(VIM)
library(naniar)
library(gridExtra)
library(ggplot2)
library(pROC)
library(caret)

# Load myfunctions
source("myfunctions.R")

# -------------------------------- #
# 1. CHECK MISSING MECHANISM -----
# -------------------------------- #

## 1.1. Is it MCAR, MAR or MNAR? ----
data.miss <- data.clin[,3:51]
miss.vars <- c("PAS_adolescencia_temprana", "DUP", "DTP","respuestas_perseverativas_PTV2M")
tab <- as.data.frame(matrix(NA, ncol=length(miss.vars), nrow = ncol(data.miss)))
colnames(tab) <- miss.vars
rownames(tab) <- colnames(data.miss)

set.seed(123)
for(var in miss.vars){
  data.miss[,var] <- ifelse(is.na(data.miss[,var]), 1, 0)
  p_vals <- rep(NA, ncol(data.miss))
  
  # Perform logistic regression 
  for(j in 1:ncol(data.miss)){
    if(names(data.miss[j]) == var){
      p_vals[j] <- NA
    }else{
      s <- summary(glm(as.formula(paste0(var, "~", 
                                         names(data.miss[j]))),
                       data=data.miss, 
                       family= binomial))
      p_vals[j] <- s$coefficients[2,4] 
    }
  }
  tab[, var] <- round(p_vals,2)
  data.miss <- data.clin[,3:51]
}

colnames(tab) <- c("PAS_Temprana", "DUP", "DTP","resp_pers")
tab <- tab[rowSums(tab <=0.05, na.rm = T) >0,]
kable(tab, "pipe", 
      caption="P values of each clinical variable for the logistic regression.")

# If at least one pval < 0.05, it is plausible that data is MAR bc we can predict when is missing based on others. So based on *Table 1*, we can assume that data is MAR and we'll use Multiple Imputation to handle the missing data


## 1.2. Visualize missing data pattern ----
aggr(data.clin[,-c(1,2)], 
     col = mdc(1:2), 
     numbers = TRUE, 
     sortVars = TRUE, 
     labels = names(data.clin),
     cex.axis = .7, 
     combined=F,
     only.miss =F, 
     ylab = c("Proportion of missingness", "Missingness Pattern"))

data.clin[,-1] %>% 
  gg_miss_var(show_pct = T, 
              facet = Resposta)


# -------------------- #
# 2. IMPUTE DATA -----
# -------------------- #

# We will generate Multivariate Imputations by Chained Equations (MICE) with the mice package. We'll define $100\times p$ imputations where $p$ is the proportion of cases with incomplete data among variables with missing values (Buuren 2018; Bodner 2008; von Hippel 2009; Ian R. White, Royston, and Wood 2011). In our case, the proportion of cases with missing value is $0.04368$, so we'll generate 5 imputations. 

imp.mod <- mice(data.clin[,-1], 
                seed  = 3, 
                m = nimpute(data.clin[,-1], method = "avg_gt_0"),
                print=F)
data.clin.imputed <- complete(imp.mod, "long", include = TRUE)
data.clin.imputed <- data.clin.imputed %>% 
  mutate(.imp = as.factor(.imp))

### Examine imputed values to validate the imputation
p1 <- bwplot(imp.mod, DUP ~ .imp)
p2 <- bwplot(imp.mod, DTP ~ .imp)
p3 <- bwplot(imp.mod, Reserva_Cognitiva ~ .imp)
p4 <- bwplot(imp.mod, respuestas_perseverativas_PTV2M ~ .imp)

grid.arrange(p1,p2,p3,p4, nrow=2)


# ----------------------------- #
# 3. SENSITIVITY ANALYSIS -----
# ----------------------------- #

# Define the training control for ML algorithms 
fitControl <- trainControl(
  method = 'cv',                   # k-fold cross validation
  number = 10,                      # number of folds
  savePredictions = 'final',       # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary,  # results summary function
  search="grid"
) 

## 3.1. Clinical Model ----

### 3.1.1. Logistic Regression ----
glm.clin <- with(imp.mod, glm(Resposta ~ +DUP+DTP+EEAG_Total_VB+Reserva_Cognitiva+
                                Insight+respuestas_perseverativas_PTV2M, family = "binomial"))
## Predictions
PREDLIST <- lapply(glm.clin$analyses, predict,se.fit=T, type = "response")
pred.clin <- as.data.frame(matrix(NA, nrow= nrow(data.clin),
                                  ncol= length(glm.clin$analyses)))

sens.clin.logit <- data.frame("Imp"=1:length(glm.clin$analyses),
                              "AUC"= rep(NA, length(glm.clin$analyses)),
                              "Threshold"= rep(NA, length(glm.clin$analyses)),
                              "Accuracy"= rep(NA, length(glm.clin$analyses)),
                              "Sensitivity"= rep(NA, length(glm.clin$analyses)),
                              "Specificity"= rep(NA, length(glm.clin$analyses)),
                              "Precision"= rep(NA, length(glm.clin$analyses)),
                              "Recall"= rep(NA, length(glm.clin$analyses)))

for(i in 1:length(glm.clin$analyses)) {
  pred.clin[,i] <- as.numeric(PREDLIST[[i]]$fit)
  roc.clin <- roc(data.clin$Resposta, pred.clin[,i], auc=TRUE, plot=F)
  ci.clin <- ci.auc(data.clin$Resposta, pred.clin[,i]) 
  cutoff.clin <- pROC::coords(roc.clin, x="best", best.method="youden")[[1]][1]
  PredResponse <- as.factor(ifelse(pred.clin[,i] > cutoff.clin, 1,0))
  metrics.clin <- confusionMatrix(reference=data.clin$Resposta, PredResponse, 
                                  positive="1")
  
  sens.clin.logit$AUC[i] <- ci.clin[2]
  sens.clin.logit$Threshold[i] <- cutoff.clin
  sens.clin.logit$Accuracy[i] <- metrics.clin$overall[1]
  sens.clin.logit$Sensitivity[i] <- metrics.clin$byClass[[1]]
  sens.clin.logit$Specificity[i] <- metrics.clin$byClass[[2]]
  sens.clin.logit$Precision[i] <- metrics.clin$byClass[[5]]
  sens.clin.logit$Recall[i] <- metrics.clin$byClass[[6]]
}

### 3.1.2. Naive Bayes Classifier ----
dat.clin.imp <- complete(imp.mod, "all")
formula <- as.formula(paste("Resposta2 ~", paste(selected.clin, collapse = " + ")))
sens.clin.naive <- data.frame("Imp"=1:length(glm.clin$analyses),
                              "AUC"= rep(NA, length(glm.clin$analyses)),
                              "Threshold"= rep(NA, length(glm.clin$analyses)),
                              "Accuracy"= rep(NA, length(glm.clin$analyses)),
                              "Sensitivity"= rep(NA, length(glm.clin$analyses)),
                              "Specificity"= rep(NA, length(glm.clin$analyses)),
                              "Precision"= rep(NA, length(glm.clin$analyses)),
                              "Recall"= rep(NA, length(glm.clin$analyses)))

for(i in 1:length(dat.clin.imp)){
  dat.clin.imp[[i]]$Resposta2 <- factor(ifelse(dat.clin.imp[[i]]$Resposta==1,
                                               "Pos","Neg"))
  
  set.seed(123)
  naive.clin = train(formula, data=dat.clin.imp[[i]], method="naive_bayes",
                     trControl = fitControl, metric="ROC")
  naive.pred.clin =predict(naive.clin, type = "prob")
  naive.roc.clin <- roc(data.clin.imp$Resposta, naive.pred.clin[,2], auc=TRUE, plot=F)
  naive.ci.clin <- ci.auc(data.clin.imp$Resposta, naive.pred.clin[,2]) 
  naive.cutoff.clin <- pROC::coords(naive.roc.clin, x="best", best.method="youden")[[1]]
  naive.resp.clin <- as.factor(ifelse(naive.pred.clin[,2] > naive.cutoff.clin, 1,0))
  naive.metrics.clin <- confusionMatrix(reference=data.clin.imp$Resposta, 
                                        naive.resp.clin, positive="1")
  
  
  sens.clin.naive$AUC[i] <- naive.ci.clin[2]
  sens.clin.naive$Threshold[i] <- naive.cutoff.clin
  sens.clin.naive$Accuracy[i] <- naive.metrics.clin$overall[1]
  sens.clin.naive$Sensitivity[i] <- naive.metrics.clin$byClass[[1]]
  sens.clin.naive$Specificity[i] <- naive.metrics.clin$byClass[[2]]
  sens.clin.naive$Precision[i] <- naive.metrics.clin$byClass[[5]]
  sens.clin.naive$Recall[i] <- naive.metrics.clin$byClass[[6]]
}


### 3.1.3. Gradient Boosting Machines ----
sens.clin.gbm <- data.frame("Imp"=1:length(glm.clin$analyses),
                            "AUC"= rep(NA, length(glm.clin$analyses)),
                            "Threshold"= rep(NA, length(glm.clin$analyses)),
                            "Accuracy"= rep(NA, length(glm.clin$analyses)),
                            "Sensitivity"= rep(NA, length(glm.clin$analyses)),
                            "Specificity"= rep(NA, length(glm.clin$analyses)),
                            "Precision"= rep(NA, length(glm.clin$analyses)),
                            "Recall"= rep(NA, length(glm.clin$analyses)))

for(i in 1:length(dat.clin.imp)){
  dat.clin.imp[[i]]$Resposta2 <- factor(ifelse(dat.clin.imp[[i]]$Resposta==1,
                                               "Pos","Neg"))
  
  set.seed(123)
  gbm.clin = train(formula, data=dat.clin.imp[[i]], method="gbm",
                   trControl = fitControl, metric="ROC", verbose=F)
  gbm.pred.clin =predict(gbm.clin, type = "prob")
  gbm.roc.clin <- roc(data.clin.imp$Resposta, gbm.pred.clin[,2], auc=TRUE, plot=F)
  gbm.ci.clin <- ci.auc(data.clin.imp$Resposta, gbm.pred.clin[,2]) 
  gbm.cutoff.clin <- pROC::coords(gbm.roc.clin, x="best", best.method="youden")[[1]][1]
  gbm.resp.clin <- as.factor(ifelse(gbm.pred.clin[,2] > gbm.cutoff.clin, 1,0))
  gbm.metrics.clin <- confusionMatrix(reference=data.clin.imp$Resposta, 
                                      gbm.resp.clin, positive="1")
  
  
  sens.clin.gbm$AUC[i] <- gbm.ci.clin[2]
  sens.clin.gbm$Threshold[i] <- gbm.cutoff.clin
  sens.clin.gbm$Accuracy[i] <- gbm.metrics.clin$overall[1]
  sens.clin.gbm$Sensitivity[i] <- gbm.metrics.clin$byClass[[1]]
  sens.clin.gbm$Specificity[i] <- gbm.metrics.clin$byClass[[2]]
  sens.clin.gbm$Precision[i] <- gbm.metrics.clin$byClass[[5]]
  sens.clin.gbm$Recall[i] <- gbm.metrics.clin$byClass[[6]]
}


### 3.1.4. Support Vector Machines ----
sens.clin.svm <- data.frame("Imp"=1:length(glm.clin$analyses),
                            "AUC"= rep(NA, length(glm.clin$analyses)),
                            "Threshold"= rep(NA, length(glm.clin$analyses)),
                            "Accuracy"= rep(NA, length(glm.clin$analyses)),
                            "Sensitivity"= rep(NA, length(glm.clin$analyses)),
                            "Specificity"= rep(NA, length(glm.clin$analyses)),
                            "Precision"= rep(NA, length(glm.clin$analyses)),
                            "Recall"= rep(NA, length(glm.clin$analyses)))

for(i in 1:length(dat.clin.imp)){
  dat.clin.imp[[i]]$Resposta2 <- factor(ifelse(dat.clin.imp[[i]]$Resposta==1,
                                               "Pos","Neg"))
  
  set.seed(123)
  svm.clin = train(formula, data=dat.clin.imp[[i]], method="svmLinear",
                   trControl = fitControl, metric="ROC")
  svm.pred.clin =predict(svm.clin, type = "prob")
  svm.roc.clin <- roc(data.clin.imp$Resposta, svm.pred.clin[,2], auc=TRUE, plot=F)
  svm.ci.clin <- ci.auc(data.clin.imp$Resposta, svm.pred.clin[,2]) 
  svm.cutoff.clin <- pROC::coords(svm.roc.clin, x="best", best.method="youden")[[1]][1]
  svm.resp.clin <- as.factor(ifelse(svm.pred.clin[,2] > svm.cutoff.clin, 1,0))
  svm.metrics.clin <- confusionMatrix(reference=data.clin.imp$Resposta, 
                                      svm.resp.clin, positive="1")
  
  
  sens.clin.svm$AUC[i] <- svm.ci.clin[2]
  sens.clin.svm$Threshold[i] <- svm.cutoff.clin
  sens.clin.svm$Accuracy[i] <- svm.metrics.clin$overall[1]
  sens.clin.svm$Sensitivity[i] <- svm.metrics.clin$byClass[[1]]
  sens.clin.svm$Specificity[i] <- svm.metrics.clin$byClass[[2]]
  sens.clin.svm$Precision[i] <- svm.metrics.clin$byClass[[5]]
  sens.clin.svm$Recall[i] <- svm.metrics.clin$byClass[[6]]
}

## 3.2. Clin+Gen Model ----
data.clin.imputed$ID <- rep(data.clin$ID, times=length(unique(data.clin.imputed$.imp)))
n= 5  # number of imputed datasets

# get all imputed datasets with for() loop
data.clingen.imp1 <- subset(data.clin.imputed, .imp==1)
data.clingen.imp2 <- subset(data.clin.imputed, .imp==2)
data.clingen.imp3 <- subset(data.clin.imputed, .imp==3)
data.clingen.imp4 <- subset(data.clin.imputed, .imp==4)
data.clingen.imp5 <- subset(data.clin.imputed, .imp==5)

# merge those imputed datasets with data.gen by variables ID using a loop
data.clingen.imp1 <- merge(data.clingen.imp1[,], data.gen[,-2], by="ID")
data.clingen.imp2 <- merge(data.clingen.imp2[,], data.gen[,-2], by="ID")
data.clingen.imp3 <- merge(data.clingen.imp3[,], data.gen[,-2], by="ID")
data.clingen.imp4 <- merge(data.clingen.imp4[,], data.gen[,-2], by="ID")
data.clingen.imp5 <- merge(data.clingen.imp5[,], data.gen[,-2], by="ID")

data.all <- list(data.clingen.imp1, data.clingen.imp2, data.clingen.imp3,
                 data.clingen.imp4, data.clingen.imp5)

selected.clingen <- c("DUP","DTP","EEAG_Total_VB", "Reserva_Cognitiva", "Insight",
                      "respuestas_perseverativas_PTV2M", "ASD","CP", "CPD","EA",
                      "IQ", "MIF","IL16")

### 3.2.1. Logistic Regression ----
formula <- as.formula(paste("Resposta ~", paste(selected.clingen, collapse = " + ")))
sens.clingen.logit <- data.frame("Imp"=1:n,
                                 "AUC"= rep(NA, n),
                                 "Threshold"= rep(NA, n),
                                 "Accuracy"= rep(NA, n),
                                 "Sensitivity"= rep(NA, n),
                                 "Specificity"= rep(NA, n),
                                 "Precision"= rep(NA, n),
                                 "Recall"= rep(NA, n))

for(i in 1:n){
  set.seed(123)
  logit.clingen = glm(formula,family = "binomial", data=data.all[[i]])
  logit.pred.clingen =predict(logit.clingen, type = "response")
  logit.roc.clingen <- roc(data.clingen.imp$Resposta, logit.pred.clingen, auc=TRUE, plot=F)
  logit.ci.clingen <- ci.auc(data.clingen.imp$Resposta, logit.pred.clingen) 
  logit.cutoff.clingen <- pROC::coords(logit.roc.clingen, x="best", best.method="youden")[[1]]
  logit.resp.clingen <- as.factor(ifelse(logit.pred.clingen > logit.cutoff.clingen, 1,0))
  logit.metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta, 
                                           logit.resp.clingen, positive="1")
  
  
  sens.clingen.logit$AUC[i] <- logit.ci.clingen[2]
  sens.clingen.logit$Threshold[i] <- logit.cutoff.clingen
  sens.clingen.logit$Accuracy[i] <- logit.metrics.clingen$overall[1]
  sens.clingen.logit$Sensitivity[i] <- logit.metrics.clingen$byClass[[1]]
  sens.clingen.logit$Specificity[i] <- logit.metrics.clingen$byClass[[2]]
  sens.clingen.logit$Precision[i] <- logit.metrics.clingen$byClass[[5]]
  sens.clingen.logit$Recall[i] <- logit.metrics.clingen$byClass[[6]]
}

### 3.2.2. Naive Bayes Classifier ----
formula <- as.formula(paste("Resposta2 ~", paste(selected.clingen, collapse = " + ")))
sens.clingen.naive <- data.frame("Imp"=1:n,
                                 "AUC"= rep(NA, n),
                                 "Threshold"= rep(NA, n),
                                 "Accuracy"= rep(NA, n),
                                 "Sensitivity"= rep(NA, n),
                                 "Specificity"= rep(NA, n),
                                 "Precision"= rep(NA, n),
                                 "Recall"= rep(NA, n))

for(i in 1:n){
  data.all[[i]]$Resposta2 <- factor(ifelse(data.all[[i]]$Resposta==1,
                                           "Pos","Neg"))
  
  set.seed(123)
  naive.clingen = train(formula, data=data.all[[i]], method="naive_bayes",
                        trControl = fitControl, metric="ROC")
  naive.pred.clingen =predict(naive.clingen, type = "prob")
  naive.roc.clingen <- roc(data.clingen.imp$Resposta, naive.pred.clingen[,2], auc=TRUE, plot=F)
  naive.ci.clingen <- ci.auc(data.clingen.imp$Resposta, naive.pred.clingen[,2]) 
  naive.cutoff.clingen <- pROC::coords(naive.roc.clingen, x="best", best.method="youden")[[1]]
  naive.resp.clingen <- as.factor(ifelse(naive.pred.clingen[,2] > naive.cutoff.clingen, 1,0))
  naive.metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta, 
                                           naive.resp.clingen, positive="1")
  
  
  sens.clingen.naive$AUC[i] <- naive.ci.clingen[2]
  sens.clingen.naive$Threshold[i] <- naive.cutoff.clin
  sens.clingen.naive$Accuracy[i] <- naive.metrics.clingen$overall[1]
  sens.clingen.naive$Sensitivity[i] <- naive.metrics.clingen$byClass[[1]]
  sens.clingen.naive$Specificity[i] <- naive.metrics.clingen$byClass[[2]]
  sens.clingen.naive$Precision[i] <- naive.metrics.clingen$byClass[[5]]
  sens.clingen.naive$Recall[i] <- naive.metrics.clingen$byClass[[6]]
}



### 3.2.3. Gradient Boosting Machines ----
sens.clingen.gbm <- data.frame("Imp"=1:n,
                               "AUC"= rep(NA, n),
                               "Threshold"= rep(NA, n),
                               "Accuracy"= rep(NA, n),
                               "Sensitivity"= rep(NA, n),
                               "Specificity"= rep(NA, n),
                               "Precision"= rep(NA, n),
                               "Recall"= rep(NA, n))

for(i in 1:n){
  data.all[[i]]$Resposta2 <- factor(ifelse(data.all[[i]]$Resposta==1,
                                           "Pos","Neg"))
  
  set.seed(123)
  gbm.clingen = train(formula, data=data.all[[i]], method="gbm",
                      trControl = fitControl, metric="ROC", verbose=F)
  gbm.pred.clingen =predict(gbm.clingen, type = "prob")
  gbm.roc.clingen <- roc(data.clingen.imp$Resposta, gbm.pred.clingen[,2], auc=TRUE, plot=F)
  gbm.ci.clingen <- ci.auc(data.clingen.imp$Resposta, gbm.pred.clingen[,2]) 
  gbm.cutoff.clingen <- pROC::coords(gbm.roc.clingen, x="best", best.method="youden")[[1]]
  gbm.resp.clingen <- as.factor(ifelse(gbm.pred.clingen[,2] > gbm.cutoff.clingen, 1,0))
  gbm.metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta, 
                                         gbm.resp.clingen, positive="1")
  
  
  sens.clingen.gbm$AUC[i] <- gbm.ci.clingen[2]
  sens.clingen.gbm$Threshold[i] <- gbm.cutoff.clin
  sens.clingen.gbm$Accuracy[i] <- gbm.metrics.clingen$overall[1]
  sens.clingen.gbm$Sensitivity[i] <- gbm.metrics.clingen$byClass[[1]]
  sens.clingen.gbm$Specificity[i] <- gbm.metrics.clingen$byClass[[2]]
  sens.clingen.gbm$Precision[i] <- gbm.metrics.clingen$byClass[[5]]
  sens.clingen.gbm$Recall[i] <- gbm.metrics.clingen$byClass[[6]]
}



### 3.2.4. Support Vector Machines ----
sens.clingen.svm <- data.frame("Imp"=1:n,
                               "AUC"= rep(NA, n),
                               "Threshold"= rep(NA, n),
                               "Accuracy"= rep(NA, n),
                               "Sensitivity"= rep(NA, n),
                               "Specificity"= rep(NA, n),
                               "Precision"= rep(NA, n),
                               "Recall"= rep(NA, n))

for(i in 1:n){
  data.all[[i]]$Resposta2 <- factor(ifelse(data.all[[i]]$Resposta==1,
                                           "Pos","Neg"))
  
  set.seed(123)
  svm.clingen = train(formula, data=data.all[[i]], method="svmLinear",
                      trControl = fitControl, metric="ROC")
  svm.pred.clingen =predict(svm.clingen, type = "prob")
  svm.roc.clingen <- roc(data.clingen.imp$Resposta, svm.pred.clingen[,2], auc=TRUE, plot=F)
  svm.ci.clingen <- ci.auc(data.clingen.imp$Resposta, svm.pred.clingen[,2]) 
  svm.cutoff.clingen <- pROC::coords(svm.roc.clingen, x="best", best.method="youden")[[1]]
  svm.resp.clingen <- as.factor(ifelse(svm.pred.clingen[,2] > svm.cutoff.clingen, 1,0))
  svm.metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta, 
                                         svm.resp.clingen, positive="1")
  
  
  sens.clingen.svm$AUC[i] <- svm.ci.clingen[2]
  sens.clingen.svm$Threshold[i] <- svm.cutoff.clin
  sens.clingen.svm$Accuracy[i] <- svm.metrics.clingen$overall[1]
  sens.clingen.svm$Sensitivity[i] <- svm.metrics.clingen$byClass[[1]]
  sens.clingen.svm$Specificity[i] <- svm.metrics.clingen$byClass[[2]]
  sens.clingen.svm$Precision[i] <- svm.metrics.clingen$byClass[[5]]
  sens.clingen.svm$Recall[i] <- svm.metrics.clingen$byClass[[6]]
}

## 3.3. Boxplots ----
model = rep(c("Logistic", "Naive", "GBM", "SVM"), each=10)
Data = rep(c(rep("Clinical", 5), rep("Clin+Gen", 5)), 4)
auc <- c(sens.clin.logit$AUC,
         sens.clingen.logit$AUC,
         sens.clin.naive$AUC,
         sens.clingen.naive$AUC,
         sens.clin.gbm$AUC,
         sens.clingen.gbm$AUC,
         sens.clin.svm$AUC,
         sens.clingen.svm$AUC)
data=data.frame(model, Data ,  auc)
data$model <- factor(data$model, levels = c("Logistic", "Naive", "GBM", "SVM"))
#data <- subset(data, Data == "Clinical")
ggplot(data, aes(x = model, y = auc, fill = Data)) +
  geom_boxplot() +
  labs(x = "", y = "AUC") +
  theme_minimal() + ylim(0.65, 1) +
  scale_fill_brewer(palette="BuPu") +
  theme(legend.position = "top")

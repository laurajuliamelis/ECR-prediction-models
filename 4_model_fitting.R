
# ============================================================================ #
# ===================        EOR PREDICTIVE MODELS           ================= #
# ===================            Model fitting               ================= #
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
library(stringr)
library(egg)
library(ggsignif)
library(ggsci)


# Load myfunctions
source("myfunctions.R")

# ------------------------------------- #
# 1. LOGISTIC REGRESSION AFTER MI -----
# ------------------------------------- #

## 1.1. Clinical Model ----
rms.clin <- lrm(Resposta ~ DUP+DTP+EEAG_Total_VB+Reserva_Cognitiva+
                  Insight+respuestas_perseverativas_PTV2M, 
                data= data.clin.imp, x=TRUE, y=TRUE)

### 1.1.1. Predictions ---- 
data.clin.imp$Prediction <- as.numeric(predict(rms.clin, type = "fitted"))

### 1.1.2. Validation Metrics ---- 
#### ROC Curve
roc.clin <- roc(data.clin.imp$Resposta, data.clin.imp$Prediction, auc=TRUE, plot=F)
auc.clin <- validate(rms.clin, B=50)
auc.clin <- CalculateAucFromDxy(auc.clin)
logisitc.ci.clin <- ci.auc(data.clin.imp$Resposta, data.clin.imp$Prediction, method= "b")


#### Confussion matrix 
cutoff.clin <- pROC::coords(roc.clin, x="best", best.method="youden")[[1]][1]
data.clin.imp$PredResponse <- as.factor(ifelse(data.clin.imp$Prediction > cutoff.clin, 1,0))
metrics.clin <- confusionMatrix(reference=data.clin.imp$Resposta,
                                data.clin.imp$PredResponse, positive="1")

#### Calibration plot
cal.clin <- rms::calibrate(rms.clin, B=50)
p1 = val.prob(data.clin.imp$Prediction,as.numeric(as.character(data.clin.imp$Resposta)))
title("Clinical Model")

### OR Table 
res <- round(exp(cbind(coef(rms.clin),confint.default(rms.clin))),3)
res <-  data.frame("OR"= paste0(res[,1], " (", res[,2], "-", res[,3],")"))
x <- round(get_model_stats(rms.clin)$coefs[,4:5],3)
res <- as.data.frame(cbind(res, x))
colnames(res)[1] <- c("OR (95%CI)")
rownames(res) <- names(coef(rms.clin))
kable(res, "pipe", 
      caption = "Results of Logistic Regression in OR and 95%CI")


## 1.2. Genetic Model ----
rms.gen <- lrm(Resposta ~ ASD+CP+CPD+EA+IQ+MIF+IL16, data= data.gen, x=TRUE, y=TRUE)

### 1.2.1. Predictions ---- 
data.gen$Prediction <- as.numeric(predict(rms.gen, type = "fitted"))

### 1.2.2. Validation Metrics ---- 
#### ROC Curve
roc.gen <- roc(data.gen$Resposta, data.gen$Prediction, auc=TRUE, plot=F)
auc.gen <- validate(rms.gen, B=50)
auc.gen <- CalculateAucFromDxy(auc.gen)
logisitc.ci.gen <- ci.auc(data.gen$Resposta, data.gen$Prediction, method= "b")

#### Confussion matrix 
cutoff.gen <- pROC::coords(roc.gen, x="best", best.method="youden")[[1]][1]
data.gen$PredResponse <- as.factor(ifelse(data.gen$Prediction > cutoff.gen, 1,0))
metrics.gen <- confusionMatrix(reference=data.gen$Resposta,
                               data.gen$PredResponse, positive="1")

#### Calibration plot
cal.gen <- rms::calibrate(rms.gen, B=50)
p1 = val.prob(data.gen$Prediction,as.numeric(as.character(data.gen$Resposta)))
title("Genetic Model")

### OR Table 
res <- round(exp(cbind(coef(rms.gen),confint.default(rms.gen))),2)
res <-  data.frame("OR"= paste0(res[,1], " (", res[,2], "-", res[,3],")"))
x <- round(get_model_stats(rms.gen)$coefs[,4:5],3)
res <- as.data.frame(cbind(res, x))
colnames(res)[1] <- c("OR (95%CI)")
rownames(res) <- names(coef(rms.gen))
kable(res, "pipe", 
      caption = "Results of Logistic Regression in OR and 95%CI")


## 1.3. Clin+Gen Model ----
rms.clingen <- lrm(Resposta ~ DUP+DTP+EEAG_Total_VB+Reserva_Cognitiva+
                     Insight+respuestas_perseverativas_PTV2M+ASD+CP+CPD+EA+IQ+MIF+IL16,
                   data= data.clingen.imp, x=TRUE, y=TRUE)


### 1.3.1. Predictions ---- 
data.clingen.imp$Prediction <- as.numeric(predict(rms.clingen, type = "fitted"))

### 1.3.2. Validation Metrics ---- 
#### ROC Curve
roc.clingen <- roc(data.clingen.imp$Resposta, data.clingen.imp$Prediction, auc=TRUE, plot=F)
auc.clingen <- validate(rms.clingen, B=50)
auc.clingen <- CalculateAucFromDxy(auc.clingen)
logisitc.ci.clingen <- ci.auc(data.clingen.imp$Resposta, data.clingen.imp$Prediction, method= "b")

#### Confussion matrix 
cutoff.clingen <- pROC::coords(roc.clingen, x="best", best.method="youden")[[1]][1]
data.clingen.imp$PredResponse <- as.factor(ifelse(data.clingen.imp$Prediction > cutoff.clingen, 1,0))
metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta,
                                   data.clingen.imp$PredResponse, positive="1")

#### Calibration plot
cal.clingen <- rms::calibrate(rms.clingen, B=50)
p1 = val.prob(data.clingen.imp$Prediction,as.numeric(as.character(data.clingen.imp$Resposta)))
title("Clin+Gen Model (LASSO clin + LASSO gen)")

### OR Table 
res <- round(exp(cbind(coef(rms.clingen),confint.default(rms.clingen))),2)
res <-  data.frame("OR"= paste0(res[,1], " (", res[,2], "-", res[,3],")"))
x <- round(get_model_stats(rms.clingen)$coefs[,4:5],3)
res <- as.data.frame(cbind(res, x))
colnames(res)[1] <- c("OR (95%CI)")
rownames(res) <- names(coef(rms.clingen))
kable(res, "pipe", 
      caption = "Results of Logistic Regression in OR and 95%CI")

## 1.4. Calibration plot ----
par(mfrow=c(3,1))
p1 = val.prob(data.clin.imp$Prediction,as.numeric(as.character(data.clin.imp$Resposta)))
title("Clinical Model")
p2 = val.prob(data.gen$Prediction,as.numeric(as.character(data.gen$Resposta)))
title("Genetic Model")
p4 = val.prob(data.clingen.imp$Prediction,as.numeric(as.character(data.clingen.imp$Resposta)))
title("Clin+Gen Model") 


# ----------------------- #
# 2. ML ALGORITHMS  -----
# ----------------------- #

# Define the training control
fitControl <- trainControl(
  method = 'cv',                   # k-fold cross validation
  number = 10,                      # number of folds
  savePredictions = 'final',       # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary,  # results summary function
  search="grid"
) 


## 2.1. Clinical Model ----
formula <- as.formula(paste("Resposta2 ~", paste(selected.clin, collapse = " + ")))
data.clin.imp$Resposta2 <- factor(ifelse(data.clin.imp$Resposta==1,"Pos","Neg")) #VN

### 2.1.1. Naive Bayes ----
set.seed(123)
naive.clin = train(formula, data=data.clin.imp, method="naive_bayes",
                   trControl = fitControl, metric="ROC")

#### Predictions and validation metrics
naive.pred.clin =predict(naive.clin, type = "prob")
naive.roc.clin <- roc(data.clin.imp$Resposta, naive.pred.clin[,2], auc=TRUE, plot=F)
naive.ci.clin <- ci.auc(data.clin.imp$Resposta, naive.pred.clin[,2], method= "b") 
naive.cutoff.clin <- pROC::coords(naive.roc.clin, x="best", best.method="youden")[[1]]
naive.resp.clin <- as.factor(ifelse(naive.pred.clin[,2] > naive.cutoff.clin, 1,0))
naive.metrics.clin <- confusionMatrix(reference=data.clin.imp$Resposta, 
                                      naive.resp.clin, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
naive.boot.clin <- boot(data = data.clin.imp, statistic = naive.auc, 
                        R = 100, parallel = "multicore", ncpus = 8)
naive.boot.clin.table <- as.data.frame(naive.boot.clin$t)
names(naive.boot.clin.table) <- c("training.auc","testing.auc","optimism")

naive.clean.mean.optimism <- mean(naive.boot.clin.table$optimism)
naive.clin.auc.c <- naive.roc.clin$auc - naive.clean.mean.optimism



### 2.1.2. Gradient Boosting Machine ----
set.seed(123)
gbm.clin = train(formula, data=data.clin.imp, method="gbm",
                 trControl = fitControl, metric="ROC", verbose=F)

#### Predictions and validation metrics
gbm.pred.clin =predict(gbm.clin, type = "prob")
gbm.roc.clin <- roc(data.clin.imp$Resposta, gbm.pred.clin[,2], auc=TRUE, plot=F)
gbm.ci.clin <- ci.auc(data.clin.imp$Resposta, gbm.pred.clin[,2], method= "b") 
gbm.cutoff.clin <- pROC::coords(gbm.roc.clin, x="best", best.method="youden")[[1]]
gbm.resp.clin <- as.factor(ifelse(gbm.pred.clin[,2] > gbm.cutoff.clin, 1,0))
gbm.metrics.clin <- confusionMatrix(reference=data.clin.imp$Resposta, 
                                    gbm.resp.clin, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
gbm.boot.clin <- boot(data = data.clin.imp, statistic = gbm.auc, 
                      R = 100, parallel = "multicore", ncpus = 8)
gbm.boot.clin.table <- as.data.frame(gbm.boot.clin$t)
names(gbm.boot.clin.table) <- c("training.auc","testing.auc","optimism")

gbm.clean.mean.optimism <- mean(gbm.boot.clin.table$optimism)
gbm.clin.auc.c <- gbm.roc.clin$auc - gbm.clean.mean.optimism


### 2.1.3. Support Vector Machine ----
set.seed(123)
svm.clin = train(formula, data=data.clin.imp, method="svmLinear",
                 trControl = fitControl, verbose=F, names=F)

#### Predictions and validation metrics
svm.pred.clin =predict(svm.clin, type = "prob")
svm.roc.clin <- roc(data.clin.imp$Resposta, svm.pred.clin[,2], auc=TRUE, plot=F)
svm.ci.clin <- ci.auc(data.clin.imp$Resposta, svm.pred.clin[,2], method= "b") 
svm.cutoff.clin <- pROC::coords(svm.roc.clin, x="best", best.method="youden")[[1]]
svm.resp.clin <- as.factor(ifelse(svm.pred.clin[,2] > svm.cutoff.clin, 1,0))
svm.metrics.clin <- confusionMatrix(reference=data.clin.imp$Resposta, svm.resp.clin, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
svm.boot.clin <- boot(data = data.clin.imp, statistic = svm.auc, 
                      R = 100, parallel = "multicore", ncpus = 8)
svm.boot.clin.table <- as.data.frame(svm.boot.clin$t)
names(svm.boot.clin.table) <- c("training.auc","testing.auc","optimism")

svm.clean.mean.optimism <- mean(svm.boot.clin.table$optimism)
svm.clin.auc.c <- svm.roc.clin$auc - svm.clean.mean.optimism


### 2.1.4. Model comparison. ----
table.clin <- data.frame("Model"= c("Logistic","Naive", "GBM", "SVM"),
                         "Threshold"=c(round(cutoff.clin,2),
                                       round(naive.cutoff.clin,2),
                                       round(gbm.cutoff.clin,2), 
                                       round(svm.cutoff.clin,2)),
                         "AUC"=c(paste0(round(auc.clin[12,1], 2), "-",  round(auc.clin[12,5],2)),
                                 paste0(round(naive.ci.clin,2)[2], "-",  round(naive.clin.auc.c,2)),
                                 paste0(round(gbm.ci.clin,2)[2], "-",  round(gbm.clin.auc.c,2)),
                                 paste0(round(svm.ci.clin,2)[2], "-",  round(svm.clin.auc.c,2))),
                         "Accuracy"=c(round(metrics.clin$overall[1],2),
                                      round(naive.metrics.clin$overall[1],2), 
                                      round(gbm.metrics.clin$overall[1],2),
                                      round(svm.metrics.clin$overall[1],2)),
                         "Sensitivity"=c(round(metrics.clin$byClass[[1]],2),
                                         round(naive.metrics.clin$byClass[[1]],2),
                                         round(gbm.metrics.clin$byClass[[1]],2),
                                         round(svm.metrics.clin$byClass[[1]],2)),
                         "Specificity"=c(round(metrics.clin$byClass[[2]],2), 
                                         round(naive.metrics.clin$byClass[[2]],2), 
                                         round(gbm.metrics.clin$byClass[[2]],2),
                                         round(svm.metrics.clin$byClass[[2]],2)),
                         "Precision"=c(round(metrics.clin$byClass[[5]],2),
                                       round(naive.metrics.clin$byClass[[5]],2),
                                       round(gbm.metrics.clin$byClass[[5]],2),
                                       round(svm.metrics.clin$byClass[[5]],2)),
                         "Recall"=c(round(metrics.clin$byClass[[6]],2),
                                    round(naive.metrics.clin$byClass[[6]],2),
                                    round(gbm.metrics.clin$byClass[[6]],2),
                                    round(svm.metrics.clin$byClass[[6]],2)))

## 2.2. Genetic Model ----
formula <- as.formula(paste("Resposta2 ~", paste(selected.gen, collapse = " + ")))
data.gen$Resposta2 <- factor(ifelse(data.gen$Resposta==1,"Pos","Neg")) #VN

### 2.2.1. Naive Bayes ----
set.seed(123)
naive.gen = train(formula, data=data.gen, method="naive_bayes",  
                  trControl = fitControl, metric="ROC")

#### Predictions and validation metrics
naive.pred.gen =predict(naive.gen, type = "prob")
naive.roc.gen <- roc(data.gen$Resposta, naive.pred.gen[,2], auc=TRUE, plot=F)
naive.ci.gen <- ci.auc(data.gen$Resposta, naive.pred.gen[,2]) 
naive.cutoff.gen <- pROC::coords(naive.roc.gen, x="best", best.method="youden")[[1]]
naive.resp.gen <- as.factor(ifelse(naive.pred.gen[,2] > naive.cutoff.gen, 1,0))
naive.metrics.gen <- confusionMatrix(reference=data.gen$Resposta, naive.resp.gen, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
naive.boot.gen <- boot(data = data.gen, statistic = naive.auc, 
                       R = 100, parallel = "multicore", ncpus = 8)
naive.boot.gen.table <- as.data.frame(naive.boot.gen$t)
names(naive.boot.gen.table) <- c("training.auc","testing.auc","optimism")

naive.clean.mean.optimism <- mean(naive.boot.gen.table$optimism)
naive.gen.auc.c <- naive.roc.gen$auc - naive.clean.mean.optimism


### 2.2.2. Gradient Boosting Machine ----
set.seed(123)
gbm.gen = train(formula, data=data.gen, method="gbm",  
                trControl = fitControl, metric="ROC", verbose=F)

#### Predictions and validation metrics
gbm.pred.gen =predict(gbm.gen, type = "prob")
gbm.roc.gen <- roc(data.gen$Resposta, gbm.pred.gen[,2], auc=TRUE, plot=F)
gbm.ci.gen <- ci.auc(data.gen$Resposta, gbm.pred.gen[,2]) 
gbm.cutoff.gen <- pROC::coords(gbm.roc.gen, x="best", best.method="youden")[[1]]
gbm.resp.gen <- as.factor(ifelse(gbm.pred.gen[,2] > gbm.cutoff.gen, 1,0))
gbm.metrics.gen <- confusionMatrix(reference=data.gen$Resposta, 
                                   gbm.resp.gen, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
gbm.boot.gen <- boot(data = data.gen, statistic = gbm.auc, 
                     R = 100, parallel = "multicore", ncpus = 8)
gbm.boot.gen.table <- as.data.frame(gbm.boot.gen$t)
names(gbm.boot.gen.table) <- c("training.auc","testing.auc","optimism")

gbm.clean.mean.optimism <- mean(gbm.boot.gen.table$optimism)
gbm.gen.auc.c <- gbm.roc.gen$auc - gbm.clean.mean.optimism


### 2.2.3. Support Vector Machine ----
set.seed(123)
svm.gen = train(formula, data=data.gen, method="svmLinear",
                trControl = fitControl)

#### Predictions and validation metrics
svm.pred.gen =predict(svm.gen, type = "prob")
svm.roc.gen <- roc(data.gen$Resposta, svm.pred.gen[,2], auc=TRUE, plot=F)
svm.ci.gen <- ci.auc(data.gen$Resposta, svm.pred.gen[,2]) 
svm.cutoff.gen <- pROC::coords(svm.roc.gen, x="best", best.method="youden")[[1]][1]
svm.resp.gen <- as.factor(ifelse(svm.pred.gen[,2] > svm.cutoff.gen, 1,0))
svm.metrics.gen <- confusionMatrix(reference=data.gen$Resposta, svm.resp.gen, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
svm.boot.gen <- boot(data = data.gen, statistic = svm.auc, 
                     R = 100, parallel = "multicore", ncpus = 8)
svm.boot.gen.table <- as.data.frame(svm.boot.gen$t)
names(svm.boot.gen.table) <- c("training.auc","testing.auc","optimism")

svm.clean.mean.optimism <- mean(svm.boot.gen.table$optimism)
svm.gen.auc.c <- svm.roc.gen$auc - svm.clean.mean.optimism


### 2.2.4. Model comparison. ----
table.gen <- data.frame("Model"= c("Logistic","Naive", "GBM", "SVM"),
                        "Threshold"=c(round(cutoff.gen,2),
                                      round(naive.cutoff.gen,2),
                                      round(gbm.cutoff.gen,2),
                                      round(svm.cutoff.gen,2)),
                        "AUC"=c(paste0(round(auc.gen[12,1], 2), "-",  round(auc.gen[12,5],2)),
                                paste0(round(naive.ci.gen,2)[2], "-",  round(naive.gen.auc.c,2)),
                                paste0(round(gbm.ci.gen,2)[2], "-",  round(gbm.gen.auc.c,2)),
                                paste0(round(svm.ci.gen,2)[2], "-",  round(svm.gen.auc.c,2))),
                        "Accuracy"=c(round(metrics.gen$overall[1],2),
                                     round(naive.metrics.gen$overall[1],2), 
                                     round(gbm.metrics.gen$overall[1],2), 
                                     round(svm.metrics.gen$overall[1],2)),
                        "Sensitivity"=c(round(metrics.gen$byClass[[1]],2),
                                        round(naive.metrics.gen$byClass[[1]],2),
                                        round(gbm.metrics.gen$byClass[[1]],2),
                                        round(svm.metrics.gen$byClass[[1]],2)),
                        "Specificity"=c(round(metrics.gen$byClass[[2]],2), 
                                        round(naive.metrics.gen$byClass[[2]],2), 
                                        round(gbm.metrics.gen$byClass[[2]],2),
                                        round(svm.metrics.gen$byClass[[2]],2)),
                        "Precision"=c(round(metrics.gen$byClass[[5]],2),
                                      round(naive.metrics.gen$byClass[[5]],2),
                                      round(gbm.metrics.gen$byClass[[5]],2),
                                      round(svm.metrics.gen$byClass[[5]],2)),
                        "Recall"=c(round(metrics.gen$byClass[[6]],2),
                                   round(naive.metrics.gen$byClass[[6]],2),
                                   round(gbm.metrics.gen$byClass[[6]],2),
                                   round(svm.metrics.gen$byClass[[6]],2)))


## 2.3. Clinical+Genetic Model ----
formula <- as.formula(paste("Resposta2 ~", paste(selected.clingen, collapse = " + ")))
data.clingen.imp$Resposta2 <- factor(ifelse(data.clingen.imp$Resposta==1,"Pos","Neg")) #VN

### 2.3.1. Naive Bayes ----
set.seed(123)
naive.clingen = train(formula, data=data.clingen.imp, method="naive_bayes",  
                      trControl = fitControl, metric="ROC")

#### Predictions and validation metrics
naive.pred.clingen =predict(naive.clingen, type = "prob")
naive.roc.clingen <- roc(data.clingen.imp$Resposta, naive.pred.clingen[,2], auc=TRUE, plot=F)
naive.ci.clingen <- ci.auc(data.clingen.imp$Resposta, naive.pred.clingen[,2]) 
naive.cutoff.clingen <- pROC::coords(naive.roc.clingen, x="best", best.method="youden")[[1]]
naive.resp.clingen <- as.factor(ifelse(naive.pred.clingen[,2] > naive.cutoff.clingen, 1,0))
naive.metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta, naive.resp.clingen, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
naive.boot.clingen <- boot(data = data.clingen.imp, statistic = naive.auc, 
                           R = 100, parallel = "multicore", ncpus = 8)
naive.boot.clingen.table <- as.data.frame(naive.boot.clingen$t)
names(naive.boot.clingen.table) <- c("training.auc","testing.auc","optimism")

naive.clean.mean.optimism <- mean(naive.boot.clingen.table$optimism)
naive.clingen.auc.c <- naive.roc.clingen$auc - naive.clean.mean.optimism


### 2.3.2. Gradient Boosting Machine ----
set.seed(123)
gbm.clingen = train(formula, data=data.clingen.imp, method="gbm",  
                    trControl = fitControl, metric="ROC", verbose=F)

#### Predictions and validation metrics
gbm.pred.clingen =predict(gbm.clingen, type = "prob")
gbm.roc.clingen <- roc(data.clingen.imp$Resposta, gbm.pred.clingen[,2], auc=TRUE, plot=F)
gbm.ci.clingen <- ci.auc(data.clingen.imp$Resposta, gbm.pred.clingen[,2]) 
gbm.cutoff.clingen <- pROC::coords(gbm.roc.clingen, x="best", best.method="youden")[[1]]
gbm.resp.clingen <- as.factor(ifelse(gbm.pred.clingen[,2] > gbm.cutoff.clingen, 1,0))
gbm.metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta, 
                                       gbm.resp.clingen, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
gbm.boot.clingen <- boot(data = data.clingen.imp, statistic = gbm.auc, 
                         R = 100, parallel = "multicore", ncpus = 8)
gbm.boot.clingen.table <- as.data.frame(gbm.boot.clingen$t)
names(gbm.boot.clingen.table) <- c("training.auc","testing.auc","optimism")

gbm.clean.mean.optimism <- mean(gbm.boot.clingen.table$optimism)
gbm.clingen.auc.c <- gbm.roc.clingen$auc - gbm.clean.mean.optimism


### 2.3.3. Support Vector Machine ----
set.seed(123)
svm.clingen = train(formula, data=data.clingen.imp, method="svmLinear",
                    trControl = fitControl)

#### Predictions and validation metrics
svm.pred.clingen =predict(svm.clingen, type = "prob")
svm.roc.clingen <- roc(data.clingen.imp$Resposta, svm.pred.clingen[,2], auc=TRUE, plot=F)
svm.ci.clingen <- ci.auc(data.clingen.imp$Resposta, svm.pred.clingen[,2]) 
svm.cutoff.clingen <- pROC::coords(svm.roc.clingen, x="best", best.method="youden")[[1]]
svm.resp.clingen <- as.factor(ifelse(svm.pred.clingen[,2] > svm.cutoff.clingen, 1,0))
svm.metrics.clingen <- confusionMatrix(reference=data.clingen.imp$Resposta, svm.resp.clingen, positive="1")

#### Optimism-corrected AUC = Uncorrected AUC - Mean optimism
svm.boot.clingen <- boot(data = data.clingen.imp, statistic = svm.auc, 
                         R = 100, parallel = "multicore", ncpus = 8)
svm.boot.clingen.table <- as.data.frame(svm.boot.clingen$t)
names(svm.boot.clingen.table) <- c("training.auc","testing.auc","optimism")

svm.clean.mean.optimism <- mean(svm.boot.clingen.table$optimism)
svm.clingen.auc.c <- svm.roc.clingen$auc - svm.clean.mean.optimism


### 2.3.4. Model comparison.. ----
table.clingen <- data.frame("Model"= c("Logistic","Naive", "GBM", "SVM"),
                            "Threshold"=c(round(cutoff.clingen,2),
                                          round(naive.cutoff.clingen,2),
                                          round(gbm.cutoff.clingen,2),
                                          round(svm.cutoff.clingen,2)),
                            "AUC"=c(paste0(round(auc.clingen[12,1], 2), "-",  round(auc.clingen[12,5],2)),
                                    paste0(round(naive.ci.clingen,2)[2], "-",  round(naive.clingen.auc.c,2)),
                                    paste0(round(gbm.ci.clingen,2)[2], "-",  round(gbm.clingen.auc.c,2)),
                                    paste0(round(svm.ci.clingen,2)[2], "-",  round(svm.clingen.auc.c,2))),
                            "Accuracy"=c(round(metrics.clingen$overall[1],2),
                                         round(naive.metrics.clingen$overall[1],2), 
                                         round(gbm.metrics.clingen$overall[1],2), 
                                         round(svm.metrics.clingen$overall[1],2)),
                            "Sensitivity"=c(round(metrics.clingen$byClass[[1]],2),
                                            round(naive.metrics.clingen$byClass[[1]],2),
                                            round(gbm.metrics.clingen$byClass[[1]],2),
                                            round(svm.metrics.clingen$byClass[[1]],2)),
                            "Specificity"=c(round(metrics.clingen$byClass[[2]],2), 
                                            round(naive.metrics.clingen$byClass[[2]],2), 
                                            round(gbm.metrics.clingen$byClass[[2]],2),
                                            round(svm.metrics.clingen$byClass[[2]],2)),
                            "Precision"=c(round(metrics.clingen$byClass[[5]],2),
                                          round(naive.metrics.clingen$byClass[[5]],2),
                                          round(gbm.metrics.clingen$byClass[[5]],2),
                                          round(svm.metrics.clingen$byClass[[5]],2)),
                            "Recall"=c(round(metrics.clingen$byClass[[6]],2),
                                       round(naive.metrics.clingen$byClass[[6]],2),
                                       round(gbm.metrics.clingen$byClass[[6]],2), 
                                       round(svm.metrics.clingen$byClass[[6]],2)))

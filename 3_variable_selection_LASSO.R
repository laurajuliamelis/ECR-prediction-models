
# ============================================================================ #
# ===================        EOR PREDICTIVE MODELS           ================= #
# ===================         Variable Selection             ================= #
# ============================================================================ #

rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary packages
library(dplyr)
library(glmnet)


# Logistic Lasso Regression was used for variable selection, with the optimal regularization parameter determined via 10-fold cross validation on loss funcion AUC with package glmnet.

# -------------------------- #
# 1. CLINICAL DATASET  -----
# -------------------------- #
data.clin.imp <- complete(imp.mod) # Extract first imputed dataset for the clinical

set.seed(1234)
cvfit = cv.glmnet(x = model.matrix(~.-1,data.clin.imp[,-1]),
                  y = data.clin.imp$Resposta,
                  family = "binomial",
                  type.measure = "auc")
results <- as.vector(coef(cvfit, s = "lambda.min"))
selected.clin <-  rownames(coef(cvfit, s = "lambda.min"))[which(results !=0)][-1]
selected.clin


# ------------------------- #
# 2. GENETIC DATASET  -----
# ------------------------- #

set.seed(1234)
cvfit = cv.glmnet(x = as.matrix(data.gen[,-c(1,2)]),
                  y = data.gen$Resposta,
                  family = "binomial",
                  type.measure = "auc")
results <- as.vector(coef(cvfit, s = "lambda.min"))
selected.gen <- names(data.gen[,-c(1)])[which(coef(cvfit, s = "lambda.min") !=0)][-1]
selected.gen


# --------------------------------- #
# 2. COMBINED MODELS DATASET  -----
# --------------------------------- #

data.clin.imp$ID <- data.clin$ID
data.clin.imp <- data.clin.imp[,c(ncol(data.clin.imp),1:(ncol(data.clin.imp)-1))]

# merge clinical imputed data with genetic data
data.clingen.imp <- merge(data.clin.imp, data.gen[,-2], by = "ID")

# remove variables ID from imputed datasets
data.clin.imp <- data.clin.imp[,-1]
data.clingen.imp <- data.clingen.imp[,-1]

selected.clingen <- c(selected.clin, selected.gen)
selected.clingen

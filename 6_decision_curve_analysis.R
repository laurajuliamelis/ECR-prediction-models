
# ============================================================================ #
# ===================        EOR PREDICTIVE MODELS           ================= #
# ===================       Decision Curve Analysis          ================= #
# ============================================================================ #

rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary packages
library(dcurves)
library(gtsummary)
library(dplyr)
library(tidyr)

# "By hand"
prevalence=0.7242991 #(EOR/N=155/214)
threshold=0.5        # probability threshold
net_ben = prevalence-(1-prevalence)*(threshold/(1-threshold)) # Net benefit
net_ben # 0.45


# ------------------------------ #
# 1. DCA OF PREDICTING EOR -----
# ------------------------------ #

# 1.1. Clinical Model ----
rms.clin <- lrm(Resposta ~ DUP+DTP+EEAG_Total_VB+Reserva_Cognitiva+
                  Insight+respuestas_perseverativas_PTV2M, 
                data= data.clin.imp, x=TRUE, y=TRUE)
data.clin.imp$Prediction <- as.numeric(predict(rms.clin, type = "fitted"))

p1 <- dca(Resposta ~ Prediction, 
          data= data.clin.imp, 
          thresholds = seq(0, 1, by = 0.01),
          label = list(Prediction = "Clinical Model")) %>%
  plot(smooth = TRUE)+
  ggtitle("Logistic Model (Clinical)")

# 1.2. Genetic Model ----
rms.gen <- lrm(Resposta ~ ASD+CP+CPD+EA+IQ+MIF+IL16, data= data.gen, x=TRUE, y=TRUE)
data.gen$Prediction <- as.numeric(predict(rms.gen, type = "fitted"))

p2 <- dca(Resposta ~ Prediction, 
          data= data.gen, 
          thresholds = seq(0, 1, by = 0.01),
          label = list(Prediction = "Genetic Model")) %>%
  plot(smooth = TRUE)+
  ggtitle("Logistic Model (Genetic)")

# 1.3. Clin+Gen Model ----

rms.clingen <- lrm(Resposta ~ DUP+DTP+EEAG_Total_VB+Reserva_Cognitiva+
                     Insight+respuestas_perseverativas_PTV2M+ASD+CP+CPD+EA+IQ+MIF+IL16,
                   data= data.clingen.imp, x=TRUE, y=TRUE)
data.clingen.imp$Prediction <- as.numeric(predict(rms.clingen, type = "fitted"))
p3 <- dca(Resposta ~ Prediction, 
          data= data.clingen.imp, 
          thresholds = seq(0, 1, by = 0.01),
          label = list(Prediction = "Clin+Gen Model")) %>%
  plot(smooth = TRUE)+
  ggtitle("Logistic Model (Clin+Gen)")


# Save plos in pdf ----
pdf("Plots_EOR.pdf", paper="a4", width = 8, height = 25)
grid.arrange(p1,p2,p3, nrow=3)
dev.off()



# ---------------------------------- #
# 2. DCA OF PREDICTING NON-EOR -----
# ---------------------------------- #

# 2.1. Clinical Model ----
data.clin.imp$PredictionINV <- 1- data.clin.imp$Prediction
data.clin.imp$RespostaINV <- ifelse(data.clin.imp$Resposta == 1, 0, 1)
p4 <- dca(RespostaINV ~ PredictionINV, 
          data= data.clin.imp, 
          thresholds = seq(0, 1, by = 0.01),
          label = list(Prediction = "Clinical Model")) %>%
  plot(smooth = TRUE)+
  ggtitle("Logistic Model (Clinical)")

# roc.clinINV <- roc(data.clin.imp$RespostaINV, data.clin.imp$PredictionINV, 
#                    auc=TRUE, plot=F)
# cutoff.clinINV <- pROC::coords(roc.clinINV, x="best", best.method="youden")[[1]][1]
# data.clin.imp$PredResponseINV <- as.factor(ifelse(data.clin.imp$PredictionINV > cutoff.clinINV, 1,0))
# metrics.clin <- confusionMatrix(reference=as.factor(data.clin.imp$RespostaINV),
#                                 data.clin.imp$PredResponseINV, positive="1")

# 2.2. Genetic Model ----
data.gen$PredictionINV <- 1- data.gen$Prediction
data.gen$RespostaINV <- ifelse(data.gen$Resposta == 1, 0, 1)
p5 <- dca(RespostaINV ~ PredictionINV, 
          data= data.gen, 
          thresholds = seq(0, 1, by = 0.01),
          label = list(Prediction = "Genetic Model")) %>%
  plot(smooth = TRUE)+
  ggtitle("Logistic Model (Genetic)")

# 2.3. Clin+Gen Model ----
data.clingen.imp$PredictionINV <- 1- data.clingen.imp$Prediction
data.clingen.imp$RespostaINV <- ifelse(data.clingen.imp$Resposta == 1, 0, 1)
p6 <- dca(RespostaINV ~ PredictionINV, 
          data= data.clingen.imp, 
          thresholds = seq(0, 1, by = 0.01),
          label = list(Prediction = "Clin+Gen Model")) %>%
  plot(smooth = TRUE)+
  ggtitle("Logistic Model (Clin+Gen)")+
  coord_cartesian(ylim = c(-1, 1))


# Save plos in pdf ----
pdf("Plots_NONEOR.pdf", paper="a4", width = 8, height = 25)
grid.arrange(p4,p5,p6, nrow=3)
dev.off()

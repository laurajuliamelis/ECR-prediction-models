
# ============================================================================ #
# ===================        EOR PREDICTIVE MODELS           ================= #
# ===================         Data preprocessing             ================= #
# ============================================================================ #

rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary packages
library(openxlsx)
library(knitr)
library(papeR)

# ---------------------------- #
# 1. READ AND CLEAN DATA -----
# ---------------------------- #

data <- read.xlsx("data/dades_PNS_EOR_clinical_v2.xlsx") # 335 patients
data$PNS <- NULL
colnames(data)[1:2] <- c("ID","Resposta")

## 1.1. Preprocess data ----
# Recodify Response variable to "1" if response and "0" if no response (instead of 1/2).
data$Resposta <- ifelse(data$Resposta == 1, 1, 0)

#Code as NAs categories that represent missings
data[data== "No evaluado"] <- NA
data[data== "Unknown"] <- NA
data[data== "Desconocido"] <- NA
data[data== "No valorable"] <- NA
data[data== ""] <- NA

# Convert to factor the character variables
data[,c(2:5,7:14,23,25:37)] <- lapply(data[,c(2:5,7:14,23,25:37)] , as.factor)

# Remove variables with more than 30% of missings [skim(data)]` or (colSums(is.na(data))/236)*100
data$Weight_At_birth <- NULL
data$TOTAL_FAS_correcta2años <- NULL

# Remove globuls blancs and lewis to mantain a larger sample size
data[, c(7:11,61:66)] <- NULL

# Remove FAST_Total, CGI_Gravedas, tabaco and alucinogenos
data$Tabaco_VB <- NULL
data$Alucinogenos_VB <- NULL
data$FAST_Total_VB <- NULL
data$CGI_Gravedad_VB <- NULL

# ## 1.2. Normalize genetic variables and EEAG/Insight (DON'T RUN BC THEY'RE ALREADY NORMALIZED)
# selected.gen <- c("ASD", "CP", "CPD", "EA", "IQ", "MIF", "IL16")
# scale_params <- data.frame(Variable = selected.gen,
#                            Media = sapply(data[selected.gen], mean, na.rm = TRUE),
#                            SD = sapply(data[selected.gen], sd, na.rm = TRUE))
# 
# data[, 52:ncol(data)] <- scale(data[, 52:ncol(data)])


## 1.3. Create datasets ----
var.clin <- names(data)[c(3:51)] # algunes són categòriques
var.gen <- names(data)[52:ncol(data)] # totes numèriques OK

# Remove patients without "Resposta"
data <- data[!is.na(data$Resposta),] # 236 patients

# Remove patients with NAs in clinical vars PAS TOTAL =NA
data.clin <- data[, c("ID", "Resposta", var.clin)]
data.clin <- data.clin[!is.na(data.clin$PAS_Total),] # no tenen proves PAS
data.clin <- data.clin[rowSums(is.na(data.clin[, 41:51])) != 11,] # no tenen test cognició
data.clin$Etnia <- droplevels(data.clin$Etnia)

# Remove patients without genetic information
data.gen <- data[complete.cases(data[, var.gen]), c("ID", "Resposta", var.gen)]


# ------------------------------- #
# 2. GET DESCRIPTIVE TABLES -----
# ------------------------------- #

## 2.1. Descriptive statistics: clin data ----
kable(papeR::summarize(data.clin, type ="factor")[,-3], "pipe", 
      caption = "Resum numèric de les variables categòriques.")
kable(papeR::summarize(data.clin)[, -c(7,8,9)], "pipe", 
      caption = "Resum numèric de les variables numèriques.")

## 2.2. Descriptive statistics: gen data ----
kable(papeR::summarize(data.gen, type ="factor")[,-3], "pipe", 
      caption = "Resum numèric de les variables categòriques.")
kable(papeR::summarize(data.gen)[1:20, -c(7,8,9)], "pipe", 
      caption = "Resum numèric de les variables numèriques.")


## 2.3. Descriptive statistics: clin+gen data ----
data.both <- merge(data.clin, data.gen[,-2], by="ID")
kable(papeR::summarize(data.both, type ="factor")[,-3], "pipe", 
      caption = "Resum numèric de les variables categòriques.")
kable(papeR::summarize(data.both)[1:15, -c(7,8,9)], "pipe", 
      caption = "Resum numèric de les variables numèriques.")



# ----------------------------------------------------- #
# 3. GET CHARACTERIZATION TABLE BY RESPONSE GROUP -----
# ----------------------------------------------------- #

variables_cat <- names(dataEOR)[sapply(dataEOR, is.factor)][-1]
variables_num <- names(dataEOR)[sapply(dataEOR, is.numeric)]
covariates <- c(variables_cat, variables_num)

data1 <- dataEOR[dataEOR$EOR == 1,]
data2 <- dataEOR[dataEOR$EOR == 0,]


## Inicialitzar vectors
c <- vector()
n <- vector()
p <- vector()
c1 <- vector()
c2 <- vector()
miss <- vector()

## Iterar per a cada variable 
# Linia 124: paste0("$\\chi^2$=",round(chisq.test(table(dataEOR$EOR,dataEOR[,i]))$statistic[[1]],2), "", paste0(", p=",round(chisq.test(table(dataEOR$EOR,dataEOR[,i]))$p.value,3)))
for (i in covariates){
  
  if(i %in% variables_cat){
    n <- c(n, paste0("**", i ,":**"), paste0(" ",names(table(dataEOR[,i]))))
    c <- c(c, "", paste0(table(dataEOR[,i]), " (", round((table(dataEOR[,i])/sum(table(dataEOR[,i])))*100,2), "%)"))
    c1 <- c(c1, "", paste0(table(data1[,i]), " (", round((table(data1[,i])/sum(table(data1[,i])))*100,2), "%)"))
    c2 <- c(c2, "", paste0(table(data2[,i]), " (", round((table(data2[,i])/sum(table(data2[,i])))*100,2), "%)"))
    p <- c(p, "",  paste0("Fisher exact test",
                          paste0(", p=",round(fisher.test(table(dataEOR$EOR,dataEOR[,i]))$p.value,3))), 
           rep("",length(table(dataEOR[,i]))-1))
    miss <- c(miss, "", 
              paste0(sum(is.na(dataEOR[,i]))," (",round((sum(is.na(dataEOR[,i]))/nrow(dataEOR))*100,2), "%)"), 
              rep("",length(table(dataEOR[,i]))-1))
    
  }else if (i %in% variables_num){
    n <- c(n, paste0("**",i,"**, mean (SD)"))
    c <- c(c, paste0(round(mean(dataEOR[,i], na.rm=T),2), " (", round(sd(dataEOR[,i], na.rm=T),2),")"))
    c1 <- c(c1, paste0(round(mean(data1[,i], na.rm=T),2), " (", round(sd(data1[,i], na.rm=T),2),")"))
    c2 <- c(c2,paste0(round(mean(data2[,i], na.rm=T),2), " (", round(sd(data2[,i], na.rm=T),2),")"))
    p <- c(p, paste0("W=",wilcox.test(dataEOR[,i]~dataEOR$EOR)$statistic[[1]], "",
                     ifelse(wilcox.test(dataEOR[,i]~dataEOR$EOR)$p.value < 0.001, ", p<0.001", 
                            paste0(", p=",round(wilcox.test(dataEOR[,i]~dataEOR$EOR)$p.value,3)))))
    miss <- c(miss, paste0(sum(is.na(dataEOR[,i]))," (",round((sum(is.na(dataEOR[,i]))/nrow(dataEOR))*100,2),"%)"))
  }
}
df2 <-data.frame("v" = n, "c" = c, "c1" = c1, "c2" = c2, "p" = p, "miss" = miss)
rownames(df2) <- NULL
colnames(df2) <- c("**Variables**",
                   paste0("**All (n=",nrow(dataEOR),")**"),
                   paste0("**EOR = 1 (n=",nrow(data1),")**"),
                   paste0("**EOR = 0 (n=",nrow(data2),")**"),
                   "**Statistic, p-value**",
                   "**Missings**")

write.xlsx(df2, "EOR_characterization.xlsx")
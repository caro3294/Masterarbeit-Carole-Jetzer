# Gentamicin nlmixr k10 and V
# autor: Carole Jetzer
# next run Gentamicin_3_saemix.R

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)

## ----libraries, cache=F, message=F, warning=F, include=F-----------------

#Installiere alle n?tigen Packages falls noch nicht installiert und lade sie dann
if(!("saemix" %in% rownames(installed.packages()))){
  install.packages("saemix")
}
library(saemix)
if(!("readxl" %in% rownames(installed.packages()))){
  install.packages("readxl")
}
library(readxl)
if(!("nlme" %in% rownames(installed.packages()))){
  install.packages("nlme")
}
library(nlme)
if(!("lattice" %in% rownames(installed.packages()))){
  install.packages("lattice")
}
library(lattice)
if(!("chron" %in% rownames(installed.packages()))){
  install.packages("chron")
}
library(chron)
if(!("deSolve" %in% rownames(installed.packages()))){
  install.packages("deSolve")
}
if(!("ggrepel" %in% rownames(installed.packages()))){
  install.packages("ggrepel")
}
library(ggrepel)
library(deSolve)
library(gridExtra)
library(RxODE)
library(MASS)
library(data.table)
library(ggplot2)
library(nlmixr)
library(xpose)
library(xpose.nlmixr)
library(shinyMixR)
library(chron)


rm(list=ls()) #Entferne alle gespeicherten Elemente f?r reproduzierbare Scripts
par.default <- par(no.readonly = TRUE) #Speicher die standart Einstellungen f?r die Plots
setwd("C:/Users/carol/OneDrive/Studium/ETH/Master/Masterarbeit/R") #Setze den Working-directory

mixr_dat <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\Gentamicin_dataframe.csv", sep=",")

##Covariaten
dfCOV <- data.frame(ID=mixr_dat$ID,
                    lAGE=mixr_dat$lAGE,
                    lKOF=mixr_dat$lKOF,
                    lHEIGHT=mixr_dat$lHEIGHT,
                    lWT=mixr_dat$lWT,
                    lCREA=mixr_dat$lCREA,
                    lUREA=mixr_dat$lUREA,
                    lGA=mixr_dat$lGA,
                    SEX=mixr_dat$SEX
                    
)

dfCOV$lCREA[dfCOV$lCREA==0] <- NA
dfCOV$lUREA[dfCOV$lUREA==0] <- NA
dfCOV$lGA[dfCOV$lGA==0] <- NA
dfCOV$SEX <- as.numeric(dfCOV$SEX) #male= 2, female=1

#Berechnung der differenz vom mean
dfCOV$dAGE <- dfCOV$lAGE - mean(dfCOV$lAGE,na.rm=T)
dfCOV$dKOF <- dfCOV$lKOF - mean(dfCOV$lKOF,na.rm=T)
dfCOV$dHEIGHT <- dfCOV$lHEIGHT - mean(dfCOV$lHEIGHT,na.rm=T)
dfCOV$dWT <- dfCOV$lWT- mean(dfCOV$lWT, na.rm=T)
dfCOV$dCREA <- dfCOV$lCREA - mean(dfCOV$lCREA, na.rm=T)
dfCOV$dUREA <- dfCOV$lUREA- mean(dfCOV$lUREA, na.rm=T)
dfCOV$dGA <- dfCOV$lGA- mean(dfCOV$lGA, na.rm=T)

mixr_dat$dAGE <- dfCOV$dAGE 
mixr_dat$dKOF <- dfCOV$dKOF 
mixr_dat$dHEIGHT <- dfCOV$dHEIGHT 
mixr_dat$dWT <- dfCOV$dWT 
mixr_dat$dCREA <- dfCOV$dCREA 
mixr_dat$dUREA <- dfCOV$dUREA 
mixr_dat$dGA <- dfCOV$dGA

dfCOV <- unique(dfCOV)
mixr_dat[is.na(mixr_dat)] <- 0 # no NA in the mixr.dat! 


###FIT
mixr.fit <- function(){ 
  ini({ #model parameter --> for comparibility reasons expressed on logarithmic scale
    lK10 <- 0.5 
    lV <- -1
    eta.lK10 ~ 0.01
    eta.lv ~ 0.01
    add.err <- 0.1
  })
  model({ #model-block: specifies the model 
    K10 <- exp(lK10 + eta.lK10) #fixed effect + random effect
    v <- exp(lV + eta.lv)
    d/dt(centr) = - K10 * centr
    cp = centr/v
    cp ~ add(add.err)
  })
}

fit_ohneC <- nlmixr(mixr.fit, mixr_dat, est="saem")
plot(augPred(fit_ohneC))
print(fit_ohneC)

#V corr with different Covariates
mixr.fit_KOF <- function(){ 
  ini({ #model parameter --> for comparibility reasons expressed on logarithmic scale
    lK10 <- 0.5 
    lV <- -1
    eta.lK10 ~ 0.01
    eta.lv ~ 0.01
    effKOF <- 1
    add.err <- 0.1
  })
  model({ #model-block: specifies the model 
    K10 <- exp(lK10 + eta.lK10) #fixed effect + random effect
    v <- exp(lV + eta.lv + effKOF*dKOF)
    d/dt(centr) = - K10 * centr
    cp = centr/v
    cp ~ add(add.err)
  })
}
fit_KOF <- nlmixr(mixr.fit_KOF, mixr_dat, est="saem")

#V corr with KOF and k10 corr with differnt Covariates
mixr.fit_KOF_AGE <- function(){ 
  ini({ #model parameter --> for comparibility reasons expressed on logarithmic scale
    lK10 <- 0.5 
    lV <- -1
    eta.lK10 ~ 0.01
    eta.lv ~ 0.01
    effKOF <- 1
    effAGE <- 1
    add.err <- 0.1
  })
  model({ #model-block: specifies the model 
    K10 <- exp(lK10 + eta.lK10 + effAGE*dAGE) #fixed effect + random effect
    v <- exp(lV + eta.lv + effKOF*dKOF)
    d/dt(centr) = - K10 * centr
    cp = centr/v
    cp ~ add(add.err)
  })
}
fit_KOF_AGE <- nlmixr(mixr.fit_KOF_AGE, mixr_dat, est="saem")

#AIC = Akaike information criterion --> AIC estimates the quality of each model, relative to each of the other models.
#CORR V : AIC(fit) / BIC(fit) ----------------------------------------------------------------------------------------------------------------

AIC(fit_KOF)
#KOF, AGE, WT, HEIGHT, GA, UREA, CREA
#AIC0  = 1719.794
#AIC V: c(2279.408, 1453.891, 2292.326, 1836.387, 1774.652, 2140.489, 1955.264)

AIC(fit_KOF_AGE)
#KOF, AGE, WT, HEIGHT, GA, UREA, CREA
#AIC V: c(2380.834, 1640.218, 2757.317, 2369.696, 2097.925, 2611.796, 2149.341)



#AIC corr V fit: KOF= 2201.327, AGE= 1465.858, WT= 2235.868, Height= 1803.44, GA= 1734.147, UREA= 2307.455, CREA= 1920.447
#AIC0  = 1719.794
#BIC corr V fir: KOF= 2221.323, AGE= 1485.854, WT= 2255.864, Height= 1823.437, GA= 1754.144, UREA= 2327.452, CREA= 1940.443 
#BIC0  = 1736.457

#CORR k10 and V with KOF --------------------------------------------------------------------------------------------------------
#AIC: KOF= 2250.981, AGE= 12250.981, WT= 2664.829, Height= 2288.639, GA= 2081.238, UREA= 2461.987, CREA= 2084.67
#AICV = 1719.794

#AIC(fit_KOF_AGE)
#calc.2LL = function(fit_KOF, nnodes.gq=8, nsd.gq=4) 
#-----------------------------------------------------------------------------------------------------------------------

##CORR V ---------------------------------------------------------------------------------------------------
###AIC
#AIC_corr_V <- c(2201.327, 1465.858, 2235.868, 1803.44, 1734.147, 2307.455, 1920.447) -1400
AIC_corr_V <- c(2279.408, 1453.891, 2292.326, 1836.387, 1774.652, 2140.489, 1955.264) -1400
AIC0 <- rep(1719.794,7) -1400
AIC_V_DIFF <- AIC0 - AIC_corr_V
AIC_V_DIFF[AIC_V_DIFF<=0] <- 0
AIC_V <- matrix(c(AIC_corr_V, AIC_V_DIFF),byrow=TRUE,nrow=2)
COV <- c("BSA","AGE","WT","HEIGHT", "GA", "UREA", "CREA")

par(mfrow=c(1,1))
barplot(AIC_V, names.arg=COV, main = "AIC after correction of V with different covariates - nlmixr fit", 
        col = c("lightblue", "grey"), ylab="AIC", ylim=c(0,1800), yaxt="n")
legend("topleft",legend=c("AIC after the correction with covariates", "AIC no covariates"),
       pch=c(22,22), pt.bg=c("lightblue","grey"))
axis(2, at=c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500), labels=c(1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900))
abline(h=AIC0, col="grey30",)
#(*) schwach signifikant (p<0.1) * signifikant (p<0.05) ** hoch signifikant (p<0.01) *** höchst signifikant (p<0.001) 

#CORR k10 and V with KOF-------------------------------------------------------------------------------------------------
###AIC
#AIC_corr_V <- c(2250.981,12250.981,2664.829,2288.639,2081.238,2461.987,2084.67) - 1400
AIC_corr_V <- c(2380.834, 1640.218, 2757.317, 2369.696, 2097.925, 2611.796, 2149.341) - 1400
AIC0 <- rep(1719.794,7) - 1400
AICV <- rep(2201.327,7) -1400
AIC_V_DIFF <- AICV-AIC_corr_V
AIC_V_DIFF[AIC_V_DIFF<=0] <- 0
AIC_V <- matrix(c(AIC_corr_V, AIC_V_DIFF),byrow=TRUE,nrow=2)
COV <- c("BSA","AGE","WT","HEIGHT", "GA", "UREA", "CREA")

barplot(AIC_V, names.arg=COV, main = "AIC - V corr with BSA and K10 with different covariate - nlmixr fit", 
        col = c("lightblue4", "lightblue"), ylab="AIC",ylim=c(0,1800), yaxt="n")
legend("topleft",legend=c("V corr with BSA, K10 corr with different covariates", "V corr BSA"),
       pch=c(22,22), pt.bg=c("lightblue4", "lightblue"))
axis(2, at=c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500), labels=c(1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900))
abline(h=AICV, col="lightblue4",)
#---------------------------------------------------------------------------------------------------------------------

dfCOV <- unique(dfCOV)

dfCOV$eta.lK10 <- unique(fit_ohneC $eta.lK10)
dfCOV$eta.lv <- unique(fit_ohneC $eta.lv)
dfCOV$eta.lK10_corrKOF <- unique(fit_KOF$eta.lK10)
dfCOV$eta.lv_corrKOF <- unique(fit_KOF$eta.lv)
dfCOV$eta.lK10_corrKOFAGE <- unique(fit_KOF_AGE$eta.lK10)
dfCOV$eta.lv_corrKOFAGE <- unique(fit_KOF_AGE$eta.lv)

results  <- as.numeric(as.character((fit_KOF_AGE$parFixed$Est.)))

nlmixr.CorrETA_k10 <- unique(fit_KOF_AGE$K10)
nlmixr.CorrETA_V <- unique(fit_KOF_AGE$v)

df.nlmixr_results <- data.frame(ID = dfCOV$ID,
                                nlmixr.CorrETA_k10 = nlmixr.CorrETA_k10,
                                nlmixr.CorrETA_V = nlmixr.CorrETA_V
)

write.csv(df.nlmixr_results,"D:/BioMed IT CSV/CSV_Gentamicin\\nlmixr_dataframe_results.csv")

mixr_dat$DV[mixr_dat$DV==0] <- NA

fit_ohneC$objf
fit_KOF$objf
fit_KOF_AGE$objf

AIC(fit_ohneC)
AIC(fit_KOF)
AIC(fit_KOF_AGE)

FITConc <- data.frame(ID=mixr_dat$ID,
                      AMT=mixr_dat$AMT,
                      CONC=as.numeric(mixr_dat$DV))

remove <- c()
for(i in seq(1,length(FITConc$ID))){
  if(all(FITConc[i,2]!=0)){
    remove <- c(remove,i)
  }
}
FITConc <- FITConc[-remove,]

FITPRED <- data.frame(ID=fit_ohneC$ID,
                      PRED=as.numeric(as.character(fit_ohneC$PRED)),
                      IPRED=as.numeric(as.character(fit_ohneC$IPRED)),
                      PRED_CORR=as.numeric(as.character(fit_KOF$PRED)),
                      IPRED_CORR=as.numeric(as.character(fit_KOF$IPRED)),
                      PRED_CORR_2=as.numeric(as.character(fit_KOF_AGE$PRED)),
                      IPRED_CORR_2=as.numeric(as.character(fit_KOF_AGE$IPRED)),
                      CONC=as.numeric(as.character(FITConc$CONC))
)

#SSR Fit1
SSRfit <- sum((predict(fit_ohneC,level=1)-FITPRED$CONC[!is.na(FITPRED$CONC)])^2,na.rm=TRUE)
print(SSRfit)

#SSR V corr with KOF
SSRfit <- sum((predict(fit_KOF,level=1)-FITPRED$CONC[!is.na(FITPRED$CONC)])^2,na.rm=TRUE)
print(SSRfit)

#SSR V corr with KOF and k10 corr with AGE
SSRfit <- sum((predict(fit_KOF_AGE,level=1)-FITPRED$CONC[!is.na(FITPRED$CONC)])^2,na.rm=TRUE)
print(SSRfit)

#Predicted vs. real Conc. Population and Individual
par(mfrow=c(1,1))
lmCONC <- lm(FITPRED$PRED ~ FITPRED$CONC)
plot(FITPRED$CONC, FITPRED$PRED , type="n",
     main ="Predicted vs. real concentrations - nlmixr Fit", las=1,
     xlab=expression(paste("measured conc [mg/L]")), xlim=c(0,30),
     ylab=expression(paste("predicted conc [mg/L]")), ylim=c(0,30)
)
abline(c(0,1))
points(FITPRED$CONC, FITPRED$PRED, pch=21, bg="grey")
points(FITPRED$CONC, FITPRED$IPRED, pch=21, bg="grey40")
abline(lmCONC, col="red") 
legend("topleft",legend=c("Population prediction (PRED)","Individual prediction (IPRED)"),
       pch=21, pt.bg=c("grey","grey40"))
text(5, 25, "adj. R-sqared = 0.8892", cex=1, col="red")

mean(c(0.8892, 0.8878, 0.8714))
sqrt(var(c(0.8892, 0.8878, 0.8714)))

     

#population prediction
par(mfrow=c(1,3))

conc1 <- FITPRED$CONC
conc.pred1 <- FITPRED$PRED
lm1 <- lm(conc1 ~ conc.pred1)
plot(FITPRED$CONC, FITPRED$PRED, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - nlmixr fit")
abline(c(0,1),col="black")
points(FITPRED$CONC, FITPRED$PRED, pch=21, bg="grey")
abline(lm1, col="red")
text(8, 30, "adj. R-sqared = 0.8892", cex=1, col="red")

conc.pred2 <- FITPRED$PRED_CORR
lm2 <- lm(conc1 ~ conc.pred2)
plot(FITPRED$CONC, FITPRED$PRED, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - nlmixr fit
V corrected with BSA")
abline(c(0,1),col="black")
points(FITPRED$CONC, FITPRED$PRED_CORR, pch=21, bg="lightblue")
abline(lm2, col="red")
text(8, 30, "adj. R-sqared = 0.9325", cex=1, col="red")

conc.pred3 <-  FITPRED$PRED_CORR_2
lm3 <- lm(conc1 ~ conc.pred3)
plot(FITPRED$CONC, FITPRED$PRED, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]", type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - nlmixr fit
V corrected with BSA
K10 corr with AGE")
abline(c(0,1),col="black")
points(FITPRED$CONC, FITPRED$PRED_CORR_2, pch=21, bg="lightblue4")
abline(lm3, col="red")
text(8, 30, "adj. R-sqared = 0.9311", cex=1, col="red")


#---------------------------------------------------------------------------------------------------------------------
#PLOT Covariates against V
par(mfrow=c(2,2))
plot(dfCOV$dKOF, dfCOV$eta.lv ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects",  pch=21, bg="grey80",
     main ="V - BSA")
lnv_KOF <- lm(dfCOV$eta.lv ~ dfCOV$dKOF)
summary(lnv_KOF)
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dKOF,na.rm=TRUE),0.95*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnv_KOF)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dKOF,na.rm=TRUE),0.78*max(dfCOV$eta.lv,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnv_KOF$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dKOF,na.rm=TRUE),0.61*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("slope = ",round(coef(lnv_KOF)[2],5),sep=""),adj=0)

plot(dfCOV$dAGE, dfCOV$eta.lv ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects",  pch=21, bg="grey80",
     main ="V - AGE")
lnv_AGE <- lm(dfCOV$eta.lv ~ dfCOV$dAGE)
summary(lnv_AGE)
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dAGE,na.rm=TRUE),0.95*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnv_AGE)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dAGE,na.rm=TRUE),0.78*max(dfCOV$eta.lv,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnv_AGE$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dAGE,na.rm=TRUE),0.61*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("slope = ",round(coef(lnv_AGE)[2],5),sep=""),adj=0)

plot(dfCOV$dWT, dfCOV$eta.lv ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  pch=21, bg="grey80",
     main ="V - WT")
lnv_WT <- lm(dfCOV$eta.lv ~ dfCOV$dWT)
summary(lnv_WT)
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT)[1]+coef(lnv_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT)[1]+coef(lnv_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dWT,na.rm=TRUE),0.95*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnv_WT)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dWT,na.rm=TRUE),0.78*max(dfCOV$eta.lv,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnv_WT$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dWT,na.rm=TRUE),0.61*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("slope = ",round(coef(lnv_WT)[2],5),sep=""),adj=0)

plot(dfCOV$dHEIGHT, dfCOV$eta.lv ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  pch=21, bg="grey80",
     main ="V - HEIGHT")
lnv_HEIGHT <- lm(dfCOV$eta.lv ~ dfCOV$dHEIGHT)
summary(lnv_HEIGHT)
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dHEIGHT,na.rm=TRUE),0.95*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnv_HEIGHT)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dHEIGHT,na.rm=TRUE),0.78*max(dfCOV$eta.lv,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnv_HEIGHT$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dHEIGHT,na.rm=TRUE),0.61*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("slope = ",round(coef(lnv_HEIGHT)[2],5),sep=""),adj=0)

plot(dfCOV$dGA, dfCOV$eta.lv ,xlab="ln(GA) - mean(ln(GA))",ylab="ln(V) random effects",  pch=21, bg="grey80",
     main ="V - GA")
lnv_GA <- lm(dfCOV$eta.lv ~ dfCOV$dGA)
summary(lnv_GA)
points(c(min(dfCOV$dGA,na.rm=TRUE),max(dfCOV$dGA,na.rm=TRUE)),
       c(coef(lnv_GA)[1]+coef(lnv_GA)[2]*min(dfCOV$dGA,na.rm=TRUE),
         coef(lnv_GA)[1]+coef(lnv_GA)[2]*max(dfCOV$dGA,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dGA,na.rm=TRUE),0.95*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnv_GA)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dGA,na.rm=TRUE),0.78*max(dfCOV$eta.lv,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnv_GA$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dGA,na.rm=TRUE),0.61*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("slope = ",round(coef(lnv_GA)[2],5),sep=""),adj=0)

plot(dfCOV$dCREA, dfCOV$eta.lv ,xlab="ln(CREA) - mean(ln(CREA))",ylab="ln(V) random effects",  pch=21, bg="grey80",
     main ="V - CREA")
lnv_CREA <- lm(dfCOV$eta.lv ~ dfCOV$dCREA)
summary(lnv_CREA)
points(c(min(dfCOV$dCREA,na.rm=TRUE),max(dfCOV$dCREA,na.rm=TRUE)),
       c(coef(lnv_CREA)[1]+coef(lnv_CREA)[2]*min(dfCOV$dCREA,na.rm=TRUE),
         coef(lnv_CREA)[1]+coef(lnv_CREA)[2]*max(dfCOV$dCREA,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dCREA,na.rm=TRUE),0.95*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnv_CREA)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dCREA,na.rm=TRUE),0.78*max(dfCOV$eta.lv,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnv_CREA$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dCREA,na.rm=TRUE),0.61*max(dfCOV$eta.lv,na.rm=TRUE),
     paste("slope = ",round(coef(lnv_CREA)[2],5),sep=""),adj=0)

plot(dfCOV$dUREA, dfCOV$eta.lv ,xlab="ln(UREA) - mean(ln(UREA))",ylab="ln(V) random effects",  pch=21, bg="grey80",
     main ="V - UREA")
lnv_UREA <- lm(dfCOV$eta.lv ~ dfCOV$dUREA)
summary(lnv_UREA)
points(c(min(dfCOV$dUREA,na.rm=TRUE),max(dfCOV$dUREA,na.rm=TRUE)),
       c(coef(lnv_UREA)[1]+coef(lnv_UREA)[2]*min(dfCOV$dUREA,na.rm=TRUE),
         coef(lnv_UREA)[1]+coef(lnv_UREA)[2]*max(dfCOV$dUREA,na.rm=TRUE)),
       type='l', col="grey60")
#text(min(dfCOV$dUREA,na.rm=TRUE),0.95*max(dfCOV$eta.lv,na.rm=TRUE),
#     paste("adj R-squared= ",round(summary(lnv_UREA)$adj.r.squared,5),sep=""),adj=0)
#text(min(dfCOV$dUREA,na.rm=TRUE),0.78*max(dfCOV$eta.lv,na.rm=TRUE), 
#     paste("SSR = ",round(sum(lnv_UREA$residuals^2),5),sep=""),adj=0)
#text(min(dfCOV$dUREA,na.rm=TRUE),0.61*max(dfCOV$eta.lv,na.rm=TRUE),
#     paste("slope = ",round(coef(lnv_UREA)[2],5),sep=""),adj=0)

#---------------------------------------------------------------------------------------------------------------------
#PLOT Covariates against K10
par(mfrow=c(2,2))
plot(dfCOV$dKOF, dfCOV$eta.lK10 ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - BSA", ylim=c(-0.3,0.3))
lnk10_KOF <- lm(dfCOV$eta.lK10 ~ dfCOV$dKOF)
summary(lnk10_KOF)
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dKOF,na.rm=TRUE),1.4*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_KOF)$adj.r.squared,5),sep=""),adj=0, bg="white")
text(min(dfCOV$dKOF,na.rm=TRUE),1.2*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_KOF$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dKOF,na.rm=TRUE),1.0*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_KOF)[2],5),sep=""),adj=0)

plot(dfCOV$dAGE, dfCOV$eta.lK10 ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - AGE", ylim=c(-0.3,0.3))
lnk10_AGE <- lm(dfCOV$eta.lK10 ~ dfCOV$dAGE)
summary(lnk10_AGE)
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dAGE,na.rm=TRUE),1.4*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_AGE)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dAGE,na.rm=TRUE),1.2*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_AGE$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dAGE,na.rm=TRUE),1.0*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_AGE)[2],5),sep=""),adj=0)

plot(dfCOV$dWT, dfCOV$eta.lK10 ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - WT", ylim=c(-0.3,0.3))
lnk10_WT <- lm(dfCOV$eta.lK10 ~ dfCOV$dWT)
summary(lnk10_WT)
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dWT,na.rm=TRUE),1.4*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_WT)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dWT,na.rm=TRUE),1.2*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_WT$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dWT,na.rm=TRUE),1.0*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_WT)[2],5),sep=""),adj=0)

plot(dfCOV$dHEIGHT, dfCOV$eta.lK10 ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - HEIGHT", ylim=c(-0.3,0.3))
lnk10_HEIGHT <- lm(dfCOV$eta.lK10 ~ dfCOV$dHEIGHT)
summary(lnk10_HEIGHT)
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dHEIGHT,na.rm=TRUE),1.4*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_HEIGHT)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dHEIGHT,na.rm=TRUE),1.20*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_HEIGHT$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dHEIGHT,na.rm=TRUE),1.0*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_HEIGHT)[2],5),sep=""),adj=0)

plot(dfCOV$dGA, dfCOV$eta.lK10 ,xlab="ln(GA) - mean(ln(GA))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - GA", ylim=c(-0.3,0.3))
lnk10_GA <- lm(dfCOV$eta.lK10 ~ dfCOV$dGA)
summary(lnk10_GA)
points(c(min(dfCOV$dGA,na.rm=TRUE),max(dfCOV$dGA,na.rm=TRUE)),
       c(coef(lnk10_GA)[1]+coef(lnk10_GA)[2]*min(dfCOV$dGA,na.rm=TRUE),
         coef(lnk10_GA)[1]+coef(lnk10_GA)[2]*max(dfCOV$dGA,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dGA,na.rm=TRUE),1.4*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_GA)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dGA,na.rm=TRUE),1.20*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_GA$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dGA,na.rm=TRUE),1.00*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_GA)[2],5),sep=""),adj=0)

plot(dfCOV$dCREA, dfCOV$eta.lK10 ,xlab="ln(CREA) - mean(ln(CREA))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - CREA", ylim=c(-0.3,0.3))
lnk10_CREA <- lm(dfCOV$eta.lK10 ~ dfCOV$dCREA)
summary(lnk10_CREA)
points(c(min(dfCOV$dCREA,na.rm=TRUE),max(dfCOV$dCREA,na.rm=TRUE)),
       c(coef(lnk10_CREA)[1]+coef(lnk10_CREA)[2]*min(dfCOV$dCREA,na.rm=TRUE),
         coef(lnk10_CREA)[1]+coef(lnk10_CREA)[2]*max(dfCOV$dCREA,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dCREA,na.rm=TRUE),1.4*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_CREA)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dCREA,na.rm=TRUE),1.2*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_CREA$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dCREA,na.rm=TRUE),1.0*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_CREA)[2],5),sep=""),adj=0)

plot(dfCOV$dUREA, dfCOV$eta.lK10 ,xlab="ln(UREA) - mean(ln(UREA))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - UREA", ylim=c(-0.3,0.3))
lnk10_UREA <- lm(dfCOV$eta.lK10 ~ dfCOV$dUREA)
summary(lnk10_UREA)
points(c(min(dfCOV$dUREA,na.rm=TRUE),max(dfCOV$dUREA,na.rm=TRUE)),
       c(coef(lnk10_UREA)[1]+coef(lnk10_UREA)[2]*min(dfCOV$dUREA,na.rm=TRUE),
         coef(lnk10_UREA)[1]+coef(lnk10_UREA)[2]*max(dfCOV$dUREA,na.rm=TRUE)),
       type='l', col="grey60")
#text(min(dfCOV$dUREA,na.rm=TRUE),0.95*max(dfCOV$eta.lK10,na.rm=TRUE),
#     paste("adj R-squared= ",round(summary(lnk10_UREA)$adj.r.squared,5),sep=""),adj=0)
#text(min(dfCOV$dUREA,na.rm=TRUE),0.78*max(dfCOV$eta.lK10,na.rm=TRUE), 
#     paste("SSR = ",round(sum(lnk10_UREA$residuals^2),5),sep=""),adj=0)
#text(min(dfCOV$dUREA,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
#     paste("slope = ",round(coef(lnk10_UREA)[2],5),sep=""),adj=0)

#-----------------------------------------------------------------------------------------------------------------------

#PLOT Covariates against V, 3 plots 
par(mfrow=c(2,3))
lnv_KOF <- lm(dfCOV$eta.lv ~ dfCOV$dKOF)
lnv_KOF_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dKOF)
lnv_KOF_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dKOF)
plot(dfCOV$dKOF, dfCOV$eta.lv ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="nlmixr V - BSA", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dKOF, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dKOF, dfCOV$eta.lv_corrKOF ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="nlmixr V - BSA - V corr with BSA", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF_corrKOF)[1]+coef(lnv_KOF_corrKOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF_corrKOF)[1]+coef(lnv_KOF_corrKOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dKOF, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dKOF, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="nlmixr V - BSA - V corr with BSA, 
     k10 corr with AGE ", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF_corrKOFAGE)[1]+coef(lnv_KOF_corrKOFAGE)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF_corrKOFAGE)[1]+coef(lnv_KOF_corrKOFAGE)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dKOF, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")

lnv_AGE <- lm(dfCOV$eta.lv ~ dfCOV$dAGE)
lnv_AGE_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dAGE)
lnv_AGE_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dAGE)
plot(dfCOV$dAGE, dfCOV$eta.lv ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="nlmixr V - AGE", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dAGE, dfCOV$eta.lv_corrKOF ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="nlmixr V - AGE - V corr with BSA", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE_corrKOF)[1]+coef(lnv_AGE_corrKOF)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE_corrKOF)[1]+coef(lnv_AGE_corrKOF)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dAGE, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dAGE, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="nlmixr V - AGE - V corr with BSA, 
     k10 corr with AGE ", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE_corrKOFAGE)[1]+coef(lnv_AGE_corrKOFAGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE_corrKOFAGE)[1]+coef(lnv_AGE_corrKOFAGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dAGE, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")

lnv_WT <- lm(dfCOV$eta.lv ~ dfCOV$dWT)
lnv_WT_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dWT)
lnv_WT_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dWT)
plot(dfCOV$dWT, dfCOV$eta.lv ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="nlmixr V - WT")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT)[1]+coef(lnv_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT)[1]+coef(lnv_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dWT, dfCOV$eta.lv_corrKOF ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="nlmixr V - WT - V corr with BSA")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT_corrKOF)[1]+coef(lnv_WT_corrKOF)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT_corrKOF)[1]+coef(lnv_WT_corrKOF)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dWT, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="nlmixr V - WT - V corr with BSA, 
     k10 corr with AGE")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT_corrKOFAGE)[1]+coef(lnv_WT_corrKOFAGE)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT_corrKOFAGE)[1]+coef(lnv_WT_corrKOFAGE)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dWT, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")


lnv_HEIGHT <- lm(dfCOV$eta.lv ~ dfCOV$dHEIGHT)
lnv_HEIGHT_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dHEIGHT)
lnv_HEIGHT_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dHEIGHT)
plot(dfCOV$dHEIGHT, dfCOV$eta.lv ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="nlmixr V - HEIGHT")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dHEIGHT, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOF ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="nlmixr V - HEIGHT - V corr with BSA")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT_corrKOF)[1]+coef(lnv_HEIGHT_corrKOF)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT_corrKOF)[1]+coef(lnv_HEIGHT_corrKOF)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="nlmixr V - HEIGHT - V corr with BSA,
     k10 corr with AGE")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT_corrKOFAGE)[1]+coef(lnv_HEIGHT_corrKOFAGE)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT_corrKOFAGE)[1]+coef(lnv_HEIGHT_corrKOFAGE)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")


#---------------------------------------------------------------------------------------------------------------------
#PLOT Covariates against K10, 3 plots 
par(mfrow=c(2,3))

lnk10_KOF <- lm(dfCOV$eta.lK10 ~ dfCOV$dKOF)
lnk10_KOF_corrKOF <- lm(dfCOV$eta.lK10_corrKOF ~ dfCOV$dKOF)
lnk10_KOF_corrKOFAGE <- lm(dfCOV$eta.lK10_corrKOFAGE ~ dfCOV$dKOF)
plot(dfCOV$dKOF, dfCOV$eta.lK10 ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - BSA")
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dKOF, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dKOF, dfCOV$eta.lK10_corrKOF ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - BSA - V corr with BSA")
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF_corrKOF)[1]+coef(lnk10_KOF_corrKOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF_corrKOF)[1]+coef(lnk10_KOF_corrKOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dKOF, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dKOF, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - BSA - V corr with BSA,
     K10 corr with AGE")
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF_corrKOFAGE)[1]+coef(lnk10_KOF_corrKOFAGE)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF_corrKOFAGE)[1]+coef(lnk10_KOF_corrKOFAGE)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dKOF, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue4")


lnk10_AGE <- lm(dfCOV$eta.lK10 ~ dfCOV$dAGE)
lnk10_AGE_corrKOF <- lm(dfCOV$eta.lK10_corrKOF ~ dfCOV$dAGE)
lnk10_AGE_corrKOFAGE <- lm(dfCOV$eta.lK10_corrKOFAGE ~ dfCOV$dAGE)
plot(dfCOV$dAGE, dfCOV$eta.lK10 ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - AGE")
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dAGE, dfCOV$eta.lK10_corrKOF ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - AGE - V corr with BSA")
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE_corrKOF)[1]+coef(lnk10_AGE_corrKOF)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE_corrKOF)[1]+coef(lnk10_AGE_corrKOF)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue")
plot(dfCOV$dAGE, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - AGE - V corr with BSA,
     K10 corr with AGE")
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE_corrKOFAGE)[1]+coef(lnk10_AGE_corrKOFAGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE_corrKOFAGE)[1]+coef(lnk10_AGE_corrKOFAGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue4")



lnk10_WT <- lm(dfCOV$eta.lK10 ~ dfCOV$dWT)
lnk10_WT_corrKOF <- lm(dfCOV$eta.lK10_corrKOF ~ dfCOV$dWT)
lnk10_WT_corrKOFAGE <- lm(dfCOV$eta.lK10_corrKOFAGE ~ dfCOV$dWT)
plot(dfCOV$dWT, dfCOV$eta.lK10 ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - WT")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dWT, dfCOV$eta.lK10_corrKOF ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - WT - V corr with BSA")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT_corrKOF)[1]+coef(lnk10_WT_corrKOF)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT_corrKOF)[1]+coef(lnk10_WT_corrKOF)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dWT, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - WT - V corr with BSA,
     K10 corr with AGE")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT_corrKOFAGE)[1]+coef(lnk10_WT_corrKOFAGE)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT_corrKOFAGE)[1]+coef(lnk10_WT_corrKOFAGE)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue4")

lnk10_HEIGHT <- lm(dfCOV$eta.lK10 ~ dfCOV$dHEIGHT)
lnk10_HEIGHT_corrKOF <- lm(dfCOV$eta.lK10_corrKOF ~ dfCOV$dHEIGHT)
lnk10_HEIGHT_corrKOFAGE <- lm(dfCOV$eta.lK10_corrKOFAGE ~ dfCOV$dHEIGHT)
plot(dfCOV$dHEIGHT, dfCOV$eta.lK10 ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - HEIGHT")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dHEIGHT, dfCOV$eta.lK10_corrKOF ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - HEIGHT - V corr with BSA")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT_corrKOF)[1]+coef(lnk10_HEIGHT_corrKOF)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT_corrKOF)[1]+coef(lnk10_HEIGHT_corrKOF)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dHEIGHT, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  type="n",
     main ="nlmixr K10 - HEIGHT - V corr with BSA,
     K10 corr with AGE")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT_corrKOFAGE)[1]+coef(lnk10_HEIGHT_corrKOFAGE)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT_corrKOFAGE)[1]+coef(lnk10_HEIGHT_corrKOFAGE)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue4")

#-----------------------------------------------------------------------------------------------------------------------



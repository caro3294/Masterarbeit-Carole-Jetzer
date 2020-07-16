# Gentamicin nlme
# autor: Carole Jetzer
# next run Gentamicin_5_doseprediction.R

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
Corr_data <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\Gentamicin_CORRdata.csv", sep=",")
dfCOV <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\Covariates_dataframe.csv", sep=",")
dfsaemix <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\saemix_dataframe.csv", sep=",")

df_nlme <- dfsaemix

df_nlme$Ti1 <- df_nlme$Ti1*60 # h in min --> hours do not work in nlme
df_nlme$Ti2 <- df_nlme$Ti2*60
df_nlme$Ti3 <- df_nlme$Ti3*60
df_nlme$TIME_P <- df_nlme$TIME_P*60
df_nlme$t2 <- df_nlme$t2*60
df_nlme$t3 <- df_nlme$t3*60

df_nlme$CONC <- df_nlme$CONC*1000 #-> mug/L
df_nlme$d1 <- df_nlme$d1*1000 #mug
df_nlme$d2 <- df_nlme$d2*1000 #mug
df_nlme$d3 <- df_nlme$d3*1000 #mug

df_nlme[is.na(df_nlme)] <- 9 # value 0 or NA is not suitable for nlme fit

Gentainfusions_nlme <- function(d1, d2, d3, Ti1, Ti2, Ti3, TIME_P, t2, t3, INF_NUMB, k10, V){
  k10 <- exp(k10)
  V <- exp(V)
  ypred <- ((d1/(Ti1*k10*V)*(1-exp(-k10*(TIME_P*(TIME_P<=Ti1)+Ti1*(TIME_P>Ti1))))*
            exp(-k10*(TIME_P-Ti1)*(TIME_P>Ti1))*1000)
            +
            (d2/(Ti2*k10*V)*(1-exp(-k10*((TIME_P-t2)*((TIME_P-t2)<=Ti2)+Ti2*((TIME_P-t2)>Ti2))))*
            exp(-k10*((TIME_P-t2)-Ti2)*((TIME_P-t2)>Ti2))*1000)*(INF_NUMB>=2)
            +
            (d3/(Ti3*k10*V)*(1-exp(-k10*((TIME_P-t3)*((TIME_P-t3)<=Ti3)+Ti3*((TIME_P-t3)>Ti3))))*
            exp(-k10*((TIME_P-t3)-Ti3)*((TIME_P-t3)>Ti3))*1000)*(INF_NUMB==3))
  return(ypred)
}

df_nlme_fit <- groupedData(CONC ~ d1 + d2 + d3 + Ti1 + Ti2 + Ti3 + TIME_P + t2 + t3 + INF_NUMB | ID,  data = df_nlme)


#startvalues
nlme.fit1 <- nlme(CONC ~ Gentainfusions_nlme(d1, d2, d3, Ti1, Ti2, Ti3, TIME_P, 
                                             t2, t3, INF_NUMB, k10, V),
                  data = df_nlme_fit,
                  fixed  = list(k10 + V ~ 1),
                  random = k10 + V ~ 1|ID,
                  groups = ~ ID ,
                  start  = list(fixed=c(log(0.05), log(1000))),
                  na.action = na.pass,
                  verbose = TRUE, 
                  control=list(returnObject= TRUE, msMaxIter=250, pnlsMaxIter=7,
                               msMaxIter =500, minScale= 0.001, tolerance=1e-6, 
                               niterEM=25, pnlsTol = 0.001, msTol=1e-7, msVerbose=TRUE))


par(mfrow=c(1,1))
conc <- df_nlme$CONC/1000
conc.pred <- nlme.fit1$fitted[,1]/1000
lm <- lm(conc ~ conc.pred)
plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc [mg/L]",type="n" ,ylab="predicted conc [mg/L]",
     main="Predicted vs. real concentrations - nlme Fit",xlim=c(0,30), ylim=c(0,30))
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, pch=21, bg="grey")
points(df_nlme$CONC/1000, nlme.fit1$fitted[,2]/1000, pch=21, bg="grey40")
abline(lm, col="red")
legend("topleft",legend=c("Population prediction (PRED)","Individual prediction (IPRED"),
       pch=c(21,21),pt.bg =c("grey","grey40"))
text(5, 25, "adj. R-sqared = 0.8878", cex=1, col="red")


##CORR COVARITIES --> try out different Covarities and see which one fits best
nlme.fit_CORR_KOF <- nlme(CONC~Gentainfusions_nlme(d1, d2, d3, Ti1, Ti2, Ti3, 
                                                   TIME_P, t2, t3, INF_NUMB, k10, V), 
                  data=df_nlme_fit,
                  fixed = list(k10 ~ 1, V ~ dKOF),
                  random = k10 + V ~ 1|ID,
                  groups = ~ ID ,
                  start = c(log(0.05),log(1000) ,log(1)),
                  verbose = TRUE, 
                  control=list(returnObject= TRUE, msMaxIter=250, pnlsMaxIter=7,
                               msMaxIter =500, minScale= 0.001, tolerance=1e-6, 
                               niterEM=25, pnlsTol = 0.001, msTol=1e-7, msVerbose=TRUE))


#AIC ohne Covariaten: 3651.568
#CORR of V
#p-value: KOF= <.0001, AGE= 0.0003, WT= <.0001, Height= <.0001, GA= 0.0257, UREA= 0.3763, CREA= 0.0033

###CORR of V--------------------------------------------------------------------------------------------
#AIC 
AIC_corr_V <- c(3614.832, 3640.721, 3615.035, 3626.867, 3648.617, 3652.786, 3644.924) -3600
AIC0 <- rep(3651.568,7) -3600
AIC_V_DIFF <- AIC0-AIC_corr_V
AIC_V_DIFF[AIC_V_DIFF<=0] <- 0
AIC_V <- matrix(c(AIC_corr_V, AIC_V_DIFF),byrow=TRUE,nrow=2)
COV <- c("BSA","AGE","WT","HEIGHT", "GA", "UREA", "CREA")

barplot(AIC_V, names.arg=COV, main = "AIC - correction of V with different covariates - nlme fit", 
        col = c("lightblue","grey"), ylab="AIC", ylim=c(0,65), yaxt="n")
legend("topleft",legend=c("AIC after the correction with covariates", "AIC no covariates"),
       pch=c(22,22), pt.bg=c("lightblue","grey"))
axis(2, at=c(0,5,10,15,20,25,30,35,40,45,50,55), labels=c(3600,3605,3610,3615,3620,3625,3630,3635,3640,3645,3650,3655))
text(0.7,55, "***", cex=2) #*** höchst signifikant (p<0.001) 
text(1.9,55, "***", cex=2)
text(3.1,55, "***", cex=2)
text(4.3,55, "***", cex=2)

text(0.7,13, 3614.832, cex=1)
text(3.1,13, 3615.035, cex=1)

#(*) schwach signifikant (p<0.1) * signifikant (p<0.05) ** hoch signifikant (p<0.01) *** höchst signifikant (p<0.001) 

###CORR of V only vs. CORR V KOF and k10 with different Cov.--------------------------------------------------------------------------------------------

nlme.fit_CORR_KOF_AGE <- nlme(CONC ~ Gentainfusions_nlme(d1, d2, d3, Ti1, Ti2, Ti3, TIME_P, t2, t3, INF_NUMB, k10, V), 
                      data = df_nlme_fit,
                      fixed = list(k10 ~ dAGE, V ~ dKOF),
                      random = k10 + V ~ 1|ID,
                      groups = ~ ID ,
                      start = c(log(0.05),log(1), log(1000),log(1)),
                      verbose = TRUE, 
                      control = list(returnObject= TRUE, msMaxIter=250, pnlsMaxIter=7,
                                   msMaxIter =500, minScale= 0.001, tolerance=1e-6, 
                                   niterEM=25, pnlsTol = 0.001, msTol=1e-7, msVerbose=TRUE))

anova(nlme.fit1, nlme.fit_CORR_KOF, nlme.fit_CORR_KOF_AGE)

#AIC:  KOF= 3612.282, AGE= 3606.640, WT= 3612.160, Height= 3613.964, GA= 3616.695, UREA= 3616.324, CREA= 3610.681
#p-value: KOF= 0.0327, AGE= 0.0014, WT= 0.0305, Height= 0.0899, GA= 0.7019, UREA=  0.472, CREA= 0.0131

AIC_corr_V <- c(3612.282, 3606.640, 3612.160, 3613.964, 3616.695, 3616.324, 3610.681) -3600
AICV <- rep(3614.841,7) -3600
AIC0 <- rep(3651.568,7) -3600
AIC_V_DIFF <- AICV-AIC_corr_V
AIC_V_DIFF[AIC_V_DIFF<=0] <- 0
AIC_V <- matrix(c(AIC_corr_V, AIC_V_DIFF),byrow=TRUE,nrow=2)
COV <- c("BSA","AGE","WT","HEIGHT", "GA", "UREA", "CREA")

barplot(AIC_V, names.arg=COV, main = "AIC - V corr with BSA and K10 with different covariates - nlme fit", 
        col = c("lightblue4","lightblue"), ylab="AIC", ylim=c(0,65), yaxt="n")
legend("topleft",legend=c("V corr with BSA, K10 corr with different covariates", "V corr BSA"),
       pch=c(22,22), pt.bg=c("lightblue4","lightblue"))
axis(2, at=c(0,5,10,15,20,25,30,35,40,45,50,55), labels=c(3600,3605,3610,3615,3620,3625,3630,3635,3640,3645,3650,3655))
text(0.7,17, "*", cex=2)  #* signifikant (p<0.05)
text(1.9,17, "**", cex=2) #** hoch signifikant (p<0.01)
text(3.1,17, "*", cex=2)  #* signifikant (p<0.05)
#----------------------------------------------------------------------------------------------------------------------------------

#population prediction
par(mfrow=c(1,1))
plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc, mg/l",ylab="predicted conc, mg/l",main="Predicted vs. real concentrations corrected with KOF - nlme Fit", type="n", xlim=c(0,30), ylim=c(0,30))
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, pch=21, bg="grey")
points(df_nlme$CONC/1000, nlme.fit_CORR_KOF$fitted[,1]/1000, pch=21, bg="lightblue")
points(df_nlme$CONC/1000, nlme.fit_CORR_KOF_AGE$fitted[,1]/1000, pch=21, bg="lightblue4")
legend("topleft",legend=c("Pop. prediction","Pop. prediction V corrected with BSA","Pop. prediction V corrected with BSA, K10 corr with AGE"),
       pch=c(21,21),pt.bg =c("grey","lightblue", "lightblue4"))

#population prediction
par(mfrow=c(1,3))

plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30))
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, pch=21, bg="grey")
text(12,30, "Population prediction - nlme fit", cex=1)

plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30))
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit_CORR_KOF$fitted[,1]/1000, pch=21, bg="lightblue")
text(12,30, "Population prediction - nlme fit
V corrected with BSA", cex=1)

plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]", type="n", xlim=c(0,30), ylim=c(0,30))
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit_CORR_KOF_AGE$fitted[,1]/1000, pch=21, bg="lightblue4")
text(12,30, "Population prediction - nlme fit
V corrected with BSA
K10 corr with AGE", cex=1)

#population prediction
par(mfrow=c(1,3))

conc1 <- df_nlme$CONC/1000
conc.pred1 <- nlme.fit1$fitted[,1]/1000
lm1 <- lm(conc1 ~ conc.pred1)
plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - nlme fit")
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, pch=21, bg="grey")
abline(lm1, col="red")
text(8, 30, "adj. R-sqared = 0.8878", cex=1, col="red")

conc.pred2 <- nlme.fit_CORR_KOF$fitted[,1]/1000
lm2 <- lm(conc1 ~ conc.pred2)
plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - nlme fit
V corrected with BSA")
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit_CORR_KOF$fitted[,1]/1000, pch=21, bg="lightblue")
abline(lm2, col="red")
text(8, 30, "adj. R-sqared = 0.9282", cex=1, col="red")

conc.pred3 <-  nlme.fit_CORR_KOF_AGE$fitted[,1]/1000
lm3 <- lm(conc1 ~ conc.pred3)
plot(df_nlme$CONC/1000, nlme.fit1$fitted[,1]/1000, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]", type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - nlme fit
V corrected with BSA
K10 corr with AGE")
abline(c(0,1),col="black")
points(df_nlme$CONC/1000, nlme.fit_CORR_KOF_AGE$fitted[,1]/1000, pch=21, bg="lightblue4")
abline(lm3, col="red")
text(8, 30, "adj. R-sqared = 0.9276", cex=1, col="red")

#random effects verteilung
nlme.fit_randoms <- nlme.fit1$coefficients$random$ID[order(as.numeric(as.character(names(nlme.fit1$coefficients$random$ID[,1]))),nlme.fit1$coefficients$random$ID[,1]),]
nlme.fit_CORR_KOF_randoms <- nlme.fit_CORR_KOF$coefficients$random$ID[order(as.numeric(as.character(names(nlme.fit_CORR_KOF$coefficients$random$ID[,1]))),nlme.fit_CORR_KOF$coefficients$random$ID[,1]),]
nlme.fit_CORR_KOF_AGE_randoms <- nlme.fit_CORR_KOF_AGE$coefficients$random$ID[order(as.numeric(as.character(names(nlme.fit_CORR_KOF_AGE$coefficients$random$ID[,1]))),nlme.fit_CORR_KOF_AGE$coefficients$random$ID[,1]),]

nlme.Corr_K10_popMean <- as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["k10.(Intercept)"]))
nlme.Corr_K10_COVeff <- as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["k10.dAGE"]))
nlme.Corr_V_popMean <- as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["V.(Intercept)"]))
nlme.Corr_V_COVeff <- as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["V.dKOF"]))

dfCOV$eta.lK10 <- as.numeric(as.character(nlme.fit_randoms[,1]))
dfCOV$eta.lv <- as.numeric(as.character(nlme.fit_randoms[,2]))
dfCOV$eta.lK10_corrKOF <- as.numeric(as.character(nlme.fit_CORR_KOF_randoms[,1]))
dfCOV$eta.lv_corrKOF <- as.numeric(as.character(nlme.fit_CORR_KOF_randoms[,2]))
dfCOV$eta.lK10_corrKOFAGE <- as.numeric(as.character(nlme.fit_CORR_KOF_AGE_randoms[,1]))
dfCOV$eta.lv_corrKOFAGE <- as.numeric(as.character(nlme.fit_CORR_KOF_AGE_randoms[,2]))

nlme.Corr_k10 <- exp(nlme.Corr_K10_popMean + nlme.Corr_K10_COVeff*dfCOV$dAGE)
nlme.CorrETA_k10 <- exp(nlme.Corr_K10_popMean + dfCOV$eta.lK10_corrKOFAGE + nlme.Corr_K10_COVeff*dfCOV$dAGE)
nlme.Corr_V <- exp(nlme.Corr_V_popMean + nlme.Corr_V_COVeff*dfCOV$dAGE)
nlme.CorrETA_V <- exp(nlme.Corr_V_popMean + dfCOV$eta.lv_corrKOFAGE  + nlme.Corr_V_COVeff*dfCOV$dKOF)


df.nlme_results <- data.frame(ID = dfCOV$ID,
                              nlme.Corr_k10 = nlme.Corr_k10,
                              nlme.CorrETA_k10 = nlme.CorrETA_k10,
                              nlme.Corr_V = nlme.Corr_V,
                              nlme.CorrETA_V = nlme.CorrETA_V,
                              nlme.Corr_K10_popMean = as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["k10.(Intercept)"])),
                              nlme.Corr_K10_COVeff = as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["k10.dAGE"])),
                              nlme.Corr_V_popMean = as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["V.(Intercept)"])),
                              nlme.Corr_V_COVeff = as.numeric(as.character(nlme.fit_CORR_KOF_AGE$coefficients$fixed["V.dKOF"]))
                              
)


write.csv(df.nlme_results,"D:/BioMed IT CSV/CSV_Gentamicin\\nlme_dataframe_results.csv")

#SSR Fit1
SSRfit <- sum((predict(nlme.fit1,level=1)-df_nlme_fit$CONC[!is.na(df_nlme_fit$CONC)])^2,na.rm=TRUE)
print(SSRfit)
SSR_1 <- sum((nlme.fit1$residuals)^2)
print(SSR_1)
SSR_2 = sum((nlme.fit1$fitted - df_nlme_fit$CONC[!is.na(df_nlme_fit$CONC)])^2, na.rm=TRUE)
print(SSR_2)

#SSR V corr with KOF
SSRfit <- sum((predict(nlme.fit_CORR_KOF,level=1)-df_nlme_fit$CONC[!is.na(df_nlme_fit$CONC)])^2,na.rm=TRUE)
print(SSRfit)

#SSR V corr with KOF and k10 corr with AGE
SSRfit <- sum((predict(nlme.fit_CORR_KOF_AGE,level=1)-df_nlme_fit$CONC[!is.na(df_nlme_fit$CONC)])^2,na.rm=TRUE)
print(SSRfit)


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
     main ="K10 - BSA")
lnk10_KOF <- lm(dfCOV$eta.lK10 ~ dfCOV$dKOF)
summary(lnk10_KOF)
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dKOF,na.rm=TRUE),0.95*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_KOF)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dKOF,na.rm=TRUE),0.78*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_KOF$residuals^2),12),sep=""),adj=0)
text(min(dfCOV$dKOF,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_KOF)[2],10),sep=""),adj=0)

plot(dfCOV$dAGE, dfCOV$eta.lK10 ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - AGE")
lnk10_AGE <- lm(dfCOV$eta.lK10 ~ dfCOV$dAGE)
summary(lnk10_AGE)
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dAGE,na.rm=TRUE),0.95*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_AGE)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dAGE,na.rm=TRUE),0.78*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_AGE$residuals^2),12),sep=""),adj=0)
text(min(dfCOV$dAGE,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_AGE)[2],10),sep=""),adj=0)

plot(dfCOV$dWT, dfCOV$eta.lK10 ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - WT")
lnk10_WT <- lm(dfCOV$eta.lK10 ~ dfCOV$dWT)
summary(lnk10_WT)
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dWT,na.rm=TRUE),0.95*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_WT)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dWT,na.rm=TRUE),0.78*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_WT$residuals^2),12),sep=""),adj=0)
text(min(dfCOV$dWT,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_WT)[2],10),sep=""),adj=0)

plot(dfCOV$dHEIGHT, dfCOV$eta.lK10 ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - HEIGHT")
lnk10_HEIGHT <- lm(dfCOV$eta.lK10 ~ dfCOV$dHEIGHT)
summary(lnk10_HEIGHT)
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dHEIGHT,na.rm=TRUE),0.95*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_HEIGHT)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dHEIGHT,na.rm=TRUE),0.78*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_HEIGHT$residuals^2),12),sep=""),adj=0)
text(min(dfCOV$dHEIGHT,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_HEIGHT)[2],10),sep=""),adj=0)

plot(dfCOV$dGA, dfCOV$eta.lK10 ,xlab="ln(GA) - mean(ln(GA))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - GA")
lnk10_GA <- lm(dfCOV$eta.lK10 ~ dfCOV$dGA)
summary(lnk10_GA)
points(c(min(dfCOV$dGA,na.rm=TRUE),max(dfCOV$dGA,na.rm=TRUE)),
       c(coef(lnk10_GA)[1]+coef(lnk10_GA)[2]*min(dfCOV$dGA,na.rm=TRUE),
         coef(lnk10_GA)[1]+coef(lnk10_GA)[2]*max(dfCOV$dGA,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dGA,na.rm=TRUE),0.95*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_GA)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dGA,na.rm=TRUE),0.78*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_GA$residuals^2),12),sep=""),adj=0)
text(min(dfCOV$dGA,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_GA)[2],10),sep=""),adj=0)

plot(dfCOV$dCREA, dfCOV$eta.lK10 ,xlab="ln(CREA) - mean(ln(CREA))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - CREA")
lnk10_CREA <- lm(dfCOV$eta.lK10 ~ dfCOV$dCREA)
summary(lnk10_CREA)
points(c(min(dfCOV$dCREA,na.rm=TRUE),max(dfCOV$dCREA,na.rm=TRUE)),
       c(coef(lnk10_CREA)[1]+coef(lnk10_CREA)[2]*min(dfCOV$dCREA,na.rm=TRUE),
         coef(lnk10_CREA)[1]+coef(lnk10_CREA)[2]*max(dfCOV$dCREA,na.rm=TRUE)),
       type='l', col="grey60")
text(min(dfCOV$dCREA,na.rm=TRUE),0.95*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("adj R-squared= ",round(summary(lnk10_CREA)$adj.r.squared,5),sep=""),adj=0)
text(min(dfCOV$dCREA,na.rm=TRUE),0.78*max(dfCOV$eta.lK10,na.rm=TRUE), 
     paste("SSR = ",round(sum(lnk10_CREA$residuals^2),12),sep=""),adj=0)
text(min(dfCOV$dCREA,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_CREA)[2],10),sep=""),adj=0)

plot(dfCOV$dUREA, dfCOV$eta.lK10 ,xlab="ln(UREA) - mean(ln(UREA))",ylab="ln(k10) random effects",  pch=21, bg="grey80",
     main ="K10 - UREA")
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
plot(dfCOV$dKOF, dfCOV$eta.lv ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="nlme V - BSA", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dKOF, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dKOF, dfCOV$eta.lv_corrKOF ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="nlme V - BSA - V corr with BSA", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF_corrKOF)[1]+coef(lnv_KOF_corrKOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF_corrKOF)[1]+coef(lnv_KOF_corrKOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dKOF, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dKOF, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="nlme V - BSA - V corr with BSA, 
     k10 corr with AGE ", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF_corrKOFAGE)[1]+coef(lnv_KOF_corrKOFAGE)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF_corrKOFAGE)[1]+coef(lnv_KOF_corrKOFAGE)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dKOF, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")

lnv_AGE <- lm(dfCOV$eta.lv ~ dfCOV$dAGE)
lnv_AGE_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dAGE)
lnv_AGE_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dAGE)
plot(dfCOV$dAGE, dfCOV$eta.lv ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="nlme V - AGE", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dAGE, dfCOV$eta.lv_corrKOF ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="nlme V - AGE - V corr with BSA", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE_corrKOF)[1]+coef(lnv_AGE_corrKOF)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE_corrKOF)[1]+coef(lnv_AGE_corrKOF)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dAGE, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dAGE, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="nlme V - AGE - V corr with BSA, 
     k10 corr with AGE ", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE_corrKOFAGE)[1]+coef(lnv_AGE_corrKOFAGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE_corrKOFAGE)[1]+coef(lnv_AGE_corrKOFAGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dAGE, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")

lnv_WT <- lm(dfCOV$eta.lv ~ dfCOV$dWT)
lnv_WT_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dWT)
lnv_WT_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dWT)
plot(dfCOV$dWT, dfCOV$eta.lv ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="nlme V - WT")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT)[1]+coef(lnv_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT)[1]+coef(lnv_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dWT, dfCOV$eta.lv_corrKOF ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="nlme V - WT - V corr with BSA")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT_corrKOF)[1]+coef(lnv_WT_corrKOF)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT_corrKOF)[1]+coef(lnv_WT_corrKOF)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dWT, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="nlme V - WT - V corr with BSA, 
     k10 corr with AGE")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT_corrKOFAGE)[1]+coef(lnv_WT_corrKOFAGE)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT_corrKOFAGE)[1]+coef(lnv_WT_corrKOFAGE)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dWT, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")


lnv_HEIGHT <- lm(dfCOV$eta.lv ~ dfCOV$dHEIGHT)
lnv_HEIGHT_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dHEIGHT)
lnv_HEIGHT_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dHEIGHT)
plot(dfCOV$dHEIGHT, dfCOV$eta.lv ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="nlme V - HEIGHT")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dHEIGHT, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOF ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="nlme V - HEIGHT - V corr with BSA")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT_corrKOF)[1]+coef(lnv_HEIGHT_corrKOF)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT_corrKOF)[1]+coef(lnv_HEIGHT_corrKOF)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="nlme V - HEIGHT - V corr with BSA,
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
     main ="nlme K10 - BSA")
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dKOF, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dKOF, dfCOV$eta.lK10_corrKOF ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - BSA - V corr with BSA")
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF_corrKOF)[1]+coef(lnk10_KOF_corrKOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF_corrKOF)[1]+coef(lnk10_KOF_corrKOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dKOF, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dKOF, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - BSA - V corr with BSA,
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
     main ="nlme K10 - AGE")
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dAGE, dfCOV$eta.lK10_corrKOF ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - AGE - V corr with BSA")
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE_corrKOF)[1]+coef(lnk10_AGE_corrKOF)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE_corrKOF)[1]+coef(lnk10_AGE_corrKOF)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue")
plot(dfCOV$dAGE, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - AGE - V corr with BSA,
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
     main ="nlme K10 - WT")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dWT, dfCOV$eta.lK10_corrKOF ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - WT - V corr with BSA")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT_corrKOF)[1]+coef(lnk10_WT_corrKOF)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT_corrKOF)[1]+coef(lnk10_WT_corrKOF)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dWT, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - WT - V corr with BSA,
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
     main ="nlme K10 - HEIGHT")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dHEIGHT, dfCOV$eta.lK10_corrKOF ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - HEIGHT - V corr with BSA")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT_corrKOF)[1]+coef(lnk10_HEIGHT_corrKOF)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT_corrKOF)[1]+coef(lnk10_HEIGHT_corrKOF)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dHEIGHT, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  type="n",
     main ="nlme K10 - HEIGHT - V corr with BSA,
     K10 corr with AGE")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT_corrKOFAGE)[1]+coef(lnk10_HEIGHT_corrKOFAGE)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT_corrKOFAGE)[1]+coef(lnk10_HEIGHT_corrKOFAGE)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue4")

#-----------------------------------------------------------------------------------------------------------------------
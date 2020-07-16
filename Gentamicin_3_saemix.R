# Gentamicin saemix
# autor: Carole Jetzer
# next run Gentamicin_4_nlme.R

## ----setup, include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)

## ----libraries, cache=F, message=F, warning=F, include=F-----------------------

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

Corr_data$ID <- Corr_data$ï..Studien.ID
Corr_data$CONC1 <- as.numeric(as.character(Corr_data$Gentamicin.Spiegel.nach.30.min))
Corr_data$CONC2 <- as.numeric(as.character(Corr_data$Gentamicin.Spiegel.nach.4h))
Corr_data$CONC3 <- as.numeric(as.character(Corr_data$Gentamicin.Spiegel.nach.24h))

#for saemix a data.frame is need with the length of the number of predictions (measured concentrations)
#3 rows --> 3 measurements. in each row define all the infusion parameters: time, duration, dose of each infusion. 

df_saemix <- data.frame(ID=rep(as.numeric(as.character(Corr_data$ID)),3), 
                        TIME_P=c(Corr_data$TimeGenta1*24, Corr_data$TimeGenta2*24 ,Corr_data$TimeGenta3*24), #h
                        CONC=c(Corr_data$CONC1, Corr_data$CONC2, Corr_data$CONC3),  #mug/L
                        t2=rep(as.numeric(as.character(Corr_data$TimeDose2*24)),3), #Time of the infusion #2-3 ; t1=0 
                        t3=rep(as.numeric(as.character(Corr_data$TimeDose3*24)),3),
                        Ti1=rep(as.numeric(as.character(Corr_data$APPLIKATIONSDAUER/60)),3), #duration of infusion #1-#3
                        Ti2=rep(as.numeric(as.character(Corr_data$APPLIKATIONSDAUER.2/60)),3),
                        Ti3=rep(as.numeric(as.character(Corr_data$APPLIKATIONSDAUER.3/60)),3),
                        d1=rep(as.numeric(as.character(Corr_data$DOSIERUNG.1)),3), #dose of infusion #1-#3 #mu
                        d2=rep(as.numeric(as.character(Corr_data$DOSIERUNG.2)),3),
                        d3=rep(as.numeric(as.character(Corr_data$DOSIERUNG.3)),3),
                        INF_NUMB=1,
                        dKOF=rep(as.numeric(as.character(dfCOV$dKOF)),3),
                        dAGE=rep(as.numeric(as.character(dfCOV$dAGE)),3),
                        dHEIGHT=rep(as.numeric(as.character(dfCOV$dHEIGHT)),3),
                        dWT=rep(as.numeric(as.character(dfCOV$dWT)),3),
                        dCREA=rep(as.numeric(as.character(dfCOV$dCREA)),3),
                        dUREA=rep(as.numeric(as.character(dfCOV$dUREA)),3),
                        dGA=rep(as.numeric(as.character(dfCOV$dGA)),3)
                        )

df_saemix = df_saemix[order(df_saemix$ID, df_saemix$TIME_P),] #order the time points 
df_saemix$INF_NUMB[df_saemix$TIME_P>=df_saemix$t2] <- 2
df_saemix$INF_NUMB[df_saemix$TIME_P>=df_saemix$t3] <- 3

#remove time point of no measurement  
remove <- c()
for(i in seq(1,length(df_saemix$ID))){
  if(all(is.na(df_saemix[i,2]))){
    remove <- c(remove,i)
  }
}
df_saemix <- df_saemix[-remove,]

write.csv(df_saemix,"D:/BioMed IT CSV/CSV_Gentamicin\\saemix_dataframe.csv")

#set all NA in the data.frame = 9 --> could not be NA or 0 that the fuction can work. 
df_saemix[is.na(df_saemix)] <- 9

#saemix model:
#psi = matrix with the number columns equal to the number of parameters in the model
#id = vector of indices matching subject index
#xdep = matrix with as many columns as predictors
"calculate the conc of the first dosing. If the dosing is infusion number two 
add the calc_conc(t) from the first dosing; if it is infusion number three
add calc_conc(t) from first and second dosing." 

tracker1 <- NULL
tracker2 <- NULL

# conc~dose/(Ti*k10*V)*(1-exp(-k10*(t*(t<=Ti)+Ti*(t>Ti))))*exp(-k10*(t-Ti)*(t>Ti))*1000

Gentainfusions <- function(psi, id, xdep){ 
  tracker1 <<- psi
  t <- xdep[,1]      # time point of the measurement
  t2 <- xdep[,2]     # time of the infusions #2-3 ; t1 = 0
  t3 <- xdep[,3]
  d1 <- xdep[,4]     # dosing of the infusions #1-#3
  d2 <- xdep[,5]
  d3 <- xdep[,6]
  Ti1 <- xdep[,7]    # duration of the infusions #1-#3
  Ti2 <- xdep[,8]
  Ti3 <- xdep[,9]
  NUMB <- xdep[,10]  #number of infusion at time point of measurement
  k10 <- exp(psi[id,1])
  V <- exp(psi[id,2])
  conc_calc <- ((d1/(Ti1*k10*V)*(1-exp(-k10*(t*(t<=Ti1)+Ti1*(t>Ti1))))*exp(-k10*(t-Ti1)*(t>Ti1))*1000)
                +
                (d2/(Ti2*k10*V)*(1-exp(-k10*((t-t2)*((t-t2)<=Ti2)+Ti2*((t-t2)>Ti2))))*
                exp(-k10*((t-t2)-Ti2)*((t-t2)>Ti2))*1000)*(NUMB>=2)
                +
                (d3/(Ti3*k10*V)*(1-exp(-k10*((t-t3)*((t-t3)<=Ti3)+Ti3*((t-t3)>Ti3))))*
                exp(-k10*((t-t3)-Ti3)*((t-t3)>Ti3))*1000)*(NUMB==3)
                )
  tracker2 <<- conc_calc
  return(conc_calc)
}

#saemix data:
#name.data = which data.frame
#name.predictors = Name (or number) of the column(s) containing the predictors (the algorithm requires at least one predictor x).
#name.group = Name (or number) of the column containing the subject id.
#name.response: Name (or number) of the column containing the response variable y modeled as a function of predictor(s) x.
#name.covariates: Name (or number) of the column(s) containing the covariates, if present (otherwise missing).
#name.X = "Time" argument to specify which of the predictor will be used as the independent variable in the plots. 

saemix.data1 <- saemixData(name.data = df_saemix, name.group=c("ID"),name.response = c("CONC"),
                           name.predictors = c("TIME_P","t2","t3", "d1", "d2", "d3", "Ti1","Ti2","Ti3","INF_NUMB"),
                           name.X="TIME_P")

saemix.model1 <- saemixModel(model=Gentainfusions,
                            psi0=matrix(c(log(0.05),log(1000)), 
                            nrow=1,dimnames=list(NULL,c("k10","V"))))


saemix.fit1 <- saemix(saemix.model1, saemix.data1, saemixControl(nbdisplay = 400,
                                                                 displayProgress = TRUE,
                                                                 save.graphs = FALSE))
saemix.fit1 <- saemix.predict(saemix.fit1)
par(mfrow=c(1,1))
lm <- lm(df_saemix$CONC ~ saemix.fit1@results@ppred)
plot(df_saemix$CONC, saemix.fit1@results@ppred,xlab="measured conc [mg/L]",type="n" ,ylab="predicted conc [mg/L]",
     main="Predicted vs. real concentrations - saemix Fit",xlim=c(0,30), ylim=c(0,30))
abline(c(0,1),col="black")
points(df_saemix$CONC, saemix.fit1@results@ppred, pch=21, bg="grey")
points(df_saemix$CONC, saemix.fit1@results@ipred, pch=21, bg="grey40")
abline(lm ,col="red")
legend("topleft",legend=c("Population prediction (PRED)","Individual prediction (IPRED)"),
       pch=c(21,21),pt.bg =c("grey","grey40"))

saemix.fit1_R2 <- cor(saemix.fit1@results@ipred,df_saemix$CONC)^2 #r-square of the fit
saemix.fit1_R2

par(mfrow=c(2,2)) 
#Individual fit
plot(saemix.fit1, plot.type="individual.fit",new=FALSE,ilist=1:4,smooth=TRUE,ylog=F, 
                       pch=1, col="lightblue4",xlab="Time in hr",ylab="Gentamicin concentrations (mg/L)")

plot(saemix.fit1, plot.type="residuals.scatter")

#possible saemix fits
#Diagnostic plot: observations versus population predictions 
#par(mfrow=c(1,1)) 
#saemix.plot.obsvspred(saemix.fit1,level=0,new=FALSE)

# LL by Importance Sampling 
#saemix.plot.llis(saemix.fit1)

# Scatter plot of residuals 
#saemix.plot.scatterresiduals(saemix.fit1)

# Boxplot of random effects 
#saemix.plot.randeff(saemix.fit1)

# Relationships between parameters and covariates
#saemix.plot.parcov(saemix.fit_CORR_V)

#Relationships between parameters and covariates, on the same page 
#par(mfrow=c(3,2)) 
#saemix.plot.parcov(saemix.fit_CORR_V,new=FALSE)




#CORR V ---------------------------------------------------------------------------------------------------
saemix.data_CORR_V <- saemixData(name.data = df_saemix, name.group=c("ID"),name.response = c("CONC"),
                                   name.predictors = c("TIME_P","t2","t3", "d1", "d2", "d3", "Ti1","Ti2","Ti3","INF_NUMB"),
                                   name.covariates = c("dKOF"),
                                   name.X="TIME_P")

saemix.model_CORR_V <- saemixModel(model=Gentainfusions,
                                     psi0=matrix(c(log(0.05),log(1000),0.1,0.1), ncol=2, byrow=TRUE, dimnames=list(NULL,c("k10","V"))),
                                     covariate.model=matrix(c(0,1), ncol=2, byrow=TRUE))

saemix.fit_CORR_V <- saemix(saemix.model_CORR_V, saemix.data_CORR_V, saemixControl(nbdisplay = 400,
                                                                                         displayProgress = TRUE,
                                                                                         save.graphs = FALSE))

saemix.fit_CORR_V <- saemix.predict(saemix.fit_CORR_V)

###AIC
#AIC.SaemixObject(saemix.fit_CORR_V) #KOF= 856.6569, AGE= 883.5517, WT= 856.9821, Height= 868.8157, GA= 891.6809, UREA= 892.5345, CREA= 885.8431
#AIC.SaemixObject(saemix.fit1) # = 892.3341
AIC_corr_V <- c(856.6569,883.5517,856.9821,868.8157,891.6809,892.5345,885.8431)-820
AIC0 <- rep(892.3341,7) -820
AIC_V_DIFF <- AIC0-AIC_corr_V
AIC_V_DIFF[AIC_V_DIFF<=0] <- 0
AIC_V <- matrix(c(AIC_corr_V, AIC_V_DIFF),byrow=TRUE,nrow=2)
COV <- c("BSA","AGE","WT","HEIGHT", "GA", "UREA", "CREA")

par(mfrow=c(1,1))
barplot(AIC_V, names.arg=COV, main = "AIC after correction of V with different covariates - saemix fit", 
        col = c("lightblue","grey"), ylab="AIC", ylim=c(0,100), yaxt="n")
legend("topleft",legend=c("AIC after the correction with covariates", "AIC no covariates"),
       pch=c(22,22), pt.bg=c("lightblue", "grey"))
axis(2, at=c(0,10,20,30,40,50,60,70,80), labels=c(820,830,840,850,860,870,880,890,900))

text(0.7,32, 856.6569, cex=1)
text(3.1,32, 856.9821, cex=1)

#Corr V with KOF and K10 with different covariates ---------------------------------------------------------------------------------------------------
saemix.data_CORR_K10 <- saemixData(name.data = df_saemix, name.group=c("ID"),name.response = c("CONC"),
                                 name.predictors = c("TIME_P","t2","t3", "d1", "d2", "d3", "Ti1","Ti2","Ti3","INF_NUMB"),
                                 name.covariates = c("dKOF","dAGE"),
                                 name.X="TIME_P")

saemix.model_CORR_K10 <- saemixModel(model=Gentainfusions,
                                   psi0=matrix(c(log(0.05),log(1000),0.1,0.1), ncol=2, byrow=TRUE, dimnames=list(NULL,c("k10","V"))),
                                   covariate.model=matrix(c(0,1,1,0), ncol=2, byrow=TRUE))

saemix.fit_CORR_K10 <- saemix(saemix.model_CORR_K10, saemix.data_CORR_K10, saemixControl(nbdisplay = 400,
                                                                                   displayProgress = TRUE,
                                                                                   save.graphs = FALSE))

saemix.fit_CORR_K10 <- saemix.predict(saemix.fit_CORR_K10)

###AIC
AIC.SaemixObject(saemix.fit_CORR_V) #KOF= 856.7545, AGE= 840.6228, WT= 856.9874, Height= 858.072, GA= 858.6764, UREA= 857.5065, CREA= 855.245
AIC.SaemixObject(saemix.fit1) # = 892.3341

AIC_corr_K10 <- c(856.7545, 840.6228, 856.9874, 858.072, 858.6764, 857.5065, 855.245) -820
AIC0 <- rep(892.3341,7) - 820
AICV <- rep(856.6569,7) - 820
AIC_V_DIFF <- AICV-AIC_corr_K10
AIC_V_DIFF[AIC_V_DIFF<=0] <- 0
AIC_V <- matrix(c(AIC_corr_K10, AIC_V_DIFF),byrow=TRUE,nrow=2)
COV <- c("BSA","AGE","WT","HEIGHT", "GA", "UREA", "CREA")

par(mfrow=c(1,1))
barplot(AIC_V, names.arg=COV, main = "AIC after correction of K10 and V with different covariates - saemix fit", 
        col = c("lightblue4", "lightblue"), ylab="AIC", ylim=c(0,100), yaxt="n")
legend("topleft",legend=c("V corr with BSA, K10 corr with different covariates", "V corr BSA"),
       pch=c(22,22), pt.bg=c("lightblue4", "lightblue"))
axis(2, at=c(0,10,20,30,40,50,60,70,80), labels=c(820,30,840,850,860,870,880,890,900))
abline(h=AICV, col="lightblue4",)

##Difference in AIC-------------------------------------------------------------------------------------
AICmin <- min(AIC_corr_V)

function.deltai <- function(AIC_corr_V) {
  deltai<- (AIC_corr_V - AICmin)
  return(deltai)
}

function.deltai(AIC_corr_V)
#0.0000, 26.8948,  0.3252, 12.1588, 35.0240, 35.8776, 29.1862

# deltai: i <2: hen there is substantial support for the i-th model (or the evidence against it is worth only a bare mention), and the proposition that it is a proper description is highly probable;
# deltai: <2 i >4: then there is strong support for the i-th model;
# deltai <4 i >7: then there is considerably less support for the i-th model;
# deltai i >10: have essentially no support.

###BIC----------------------------------------------------------------------------------------------------------
BIC.SaemixObject(saemix.fit_CORR_V) #KOF= 870.3997, AGE= 897.2945, WT= 870.7248, Height= 882.5584, GA= 905.4236, UREA= 906.2773, CREA= 899.5858 
BIC.SaemixObject(saemix.fit1) # = 903.7864

BIC_corr_V <- c(870.3997, 897.2945, 870.7248, 882.5584, 905.4236, 906.2773, 899.5858)-845 
BIC0 <- rep(903.7864,7) -845
BIC_V_DIFF <- BIC0-BIC_corr_V
BIC_V_DIFF[BIC_V_DIFF<=0] <- 0
BIC_V <- matrix(c(BIC_corr_V, BIC_V_DIFF),byrow=TRUE,nrow=2)
COV <- c("BSA","AGE","WT","HEIGHT", "GA", "UREA", "CREA")

par(mfrow=c(1,1))
barplot(BIC_V, names.arg=COV, main = "BIC after correction V with different covariates - saemix fit", 
        col = c("grey","lightblue"), ylab="BIC", ylim=c(0,80), yaxt="n")
legend("topleft",legend=c("BIC after the correction with covariates", "BIC"),
       pch=c(22,22), pt.bg=c("grey","lightblue"))
axis(2, at=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65), labels=c(845,850,855,860,865,870,875,880,885,890,895,900,905,910))

text(0.7,23, 870.3997, cex=1)
text(3.1,23, 870.7248, cex=1)
#(*) schwach signifikant (p<0.1) * signifikant (p<0.05) ** hoch signifikant (p<0.01) *** höchst signifikant (p<0.001) 

#population predicion ----------------------------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(df_saemix$CONC, saemix.fit_CORR_V@results@ppred,xlab="measured conc, mg/l",type="n" ,ylab="predicted conc, mg/l",
     main="Predicted vs. real concentrations corrected with Covariates - saemix Fit",xlim=c(0,30), ylim=c(0,30))
abline(c(0,1),col="black")
points(df_saemix$CONC, saemix.fit1@results@ppred, pch=21, bg="grey")
points(df_saemix$CONC, saemix.fit_CORR_V@results@ppred, pch=21, bg="lightblue")
points(df_saemix$CONC, saemix.fit_CORR_K10@results@ppred, pch=21, bg="lightblue4")
legend("topleft",legend=c("Pop. prediction","Pop. prediction V corrected with BSA", "Pop. prediction V corrected with BSA, K10 with AGE"),
       pch=c(21,21),pt.bg =c("grey","lightblue", "lightblue4"))


saemix.fit_CORR_KOF_R2 <- cor(saemix.fit_CORR_V@results@ipred,df_saemix$CONC)^2
saemix.fit_CORR_KOF_R2
saemix.fit_CORR_KOF_AGE_R2 <- cor(saemix.fit_CORR_K10@results@ipred,df_saemix$CONC)^2
saemix.fit_CORR_KOF_AGE_R2

#population prediction
par(mfrow=c(1,3))

lm1 <- lm(df_saemix$CONC ~ saemix.fit1@results@ppred)
plot(df_saemix$CONC, saemix.fit1@results@ppred, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - saemix fit")
abline(c(0,1),col="black")
points(df_saemix$CONC, saemix.fit1@results@ppred, pch=21, bg="grey")
abline(lm1, col="red")
text(8, 30, "adj. R-sqared = 0.8714", cex=1, col="red")

lm2 <- lm(df_saemix$CONC ~ saemix.fit_CORR_V@results@ppred)
plot(df_saemix$CONC, saemix.fit1@results@ppred, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]",type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - saemix fit
V corrected with BSA")
abline(c(0,1),col="black")
points(df_saemix$CONC, saemix.fit_CORR_V@results@ppred, pch=21, bg="lightblue")
abline(lm2, col="red")
text(8, 30, "adj. R-sqared = 0.9115", cex=1, col="red")

lm3 <- lm(df_saemix$CONC ~ saemix.fit_CORR_K10@results@ppred)
plot(df_saemix$CONC, saemix.fit1@results@ppred, xlab="measured conc [mg/L]",ylab="predicted conc [mg/L]", type="n", xlim=c(0,30), ylim=c(0,30),
     main="Population prediction - saemix fit
V corrected with BSA
K10 corr with AGE")
abline(c(0,1),col="black")
points(df_saemix$CONC, saemix.fit_CORR_K10@results@ppred, pch=21, bg="lightblue4")
abline(lm3, col="red")
text(8, 30, "adj. R-sqared = 0.9114", cex=1, col="red")

#random effects verteilung -----------------------------------------------------------------------------------------
saemix.fit1_randoms <- eta(saemix.fit1)
saemix.fit_CORR_V_randoms <- eta(saemix.fit_CORR_V)
saemix.fit_CORR_K10_randoms <- eta(saemix.fit_CORR_K10)

dfCOV$eta.lK10 <- saemix.fit1_randoms[,1]
dfCOV$eta.lv <- saemix.fit1_randoms[,2]
dfCOV$eta.lK10_corrKOF <- saemix.fit_CORR_V_randoms[,1]
dfCOV$eta.lv_corrKOF <- saemix.fit_CORR_V_randoms[,2]
dfCOV$eta.lK10_corrKOFAGE <- saemix.fit_CORR_K10_randoms[,1]
dfCOV$eta.lv_corrKOFAGE <- saemix.fit_CORR_K10_randoms[,2]

results <- saemix.fit_CORR_K10@results@betas

saemix.Corr_K10_popMean <- as.numeric(as.character(results[1,]))
saemix.Corr_K10_COVeff <- as.numeric(as.character(results[2,]))
saemix.Corr_V_popMean <- as.numeric(as.character(results[3,]))
saemix.Corr_V_COVeff <- as.numeric(as.character(results[4,]))

saemix.Corr_k10 <- exp(saemix.Corr_K10_popMean + saemix.Corr_K10_COVeff*dfCOV$dAGE)
saemix.CorrETA_k10 <- exp(saemix.Corr_K10_popMean + dfCOV$eta.lK10_corrKOFAGE + saemix.Corr_K10_COVeff*dfCOV$dAGE)
saemix.Corr_V <- exp(saemix.Corr_V_popMean + saemix.Corr_V_COVeff*dfCOV$dAGE)
saemix.CorrETA_V <- exp(saemix.Corr_V_popMean + dfCOV$eta.lv_corrKOFAGE  + saemix.Corr_V_COVeff*dfCOV$dKOF)


df.saemix_results <- data.frame(ID = dfCOV$ID,
                                saemix.Corr_k10 = saemix.Corr_k10,
                                saemix.CorrETA_k10 = saemix.CorrETA_k10,
                                saemix.Corr_V = saemix.Corr_V,
                                saemix.CorrETA_V = saemix.CorrETA_V,
                                saemix.Corr_K10_popMean = as.numeric(as.character(results[1,])),
                                saemix.Corr_K10_COVeff = as.numeric(as.character(results[2,])),
                                saemix.Corr_V_popMean = as.numeric(as.character(results[3,])),
                                saemix.Corr_V_COVeff = as.numeric(as.character(results[4,]))
)

write.csv(df.saemix_results,"D:/BioMed IT CSV/CSV_Gentamicin\\saemix_dataframe_results.csv")

#SSR Fit1
SSRfit <- sum((predict(saemix.fit1,level=1)-df_saemix$CONC[!is.na(df_saemix$CONC)])^2,na.rm=TRUE)
print(SSRfit)

#SSR V corr with KOF
SSRfit <- sum((predict(saemix.fit_CORR_V,level=1)-df_saemix$CONC[!is.na(df_saemix$CONC)])^2,na.rm=TRUE)
print(SSRfit)

#SSR V corr with KOF and k10 corr with AGE
SSRfit <- sum((predict(saemix.fit_CORR_K10,level=1)-df_saemix$CONC[!is.na(df_saemix$CONC)])^2,na.rm=TRUE)
print(SSRfit)

logLik.SaemixObject(saemix.fit1, method="is")
logLik.SaemixObject(saemix.fit_CORR_V, method="is")
logLik.SaemixObject(saemix.fit_CORR_K10, method="is")

AIC.SaemixObject(saemix.fit1)
AIC.SaemixObject(saemix.fit_CORR_V)
AIC.SaemixObject(saemix.fit_CORR_K10)

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
par(mfrow=c(2,2)) #1000 #550
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
     paste("SSR = ",round(sum(lnk10_KOF$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dKOF,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_KOF)[2],5),sep=""),adj=0)

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
     paste("SSR = ",round(sum(lnk10_AGE$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dAGE,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_AGE)[2],5),sep=""),adj=0)

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
     paste("SSR = ",round(sum(lnk10_WT$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dWT,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_WT)[2],5),sep=""),adj=0)

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
     paste("SSR = ",round(sum(lnk10_HEIGHT$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dHEIGHT,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_HEIGHT)[2],5),sep=""),adj=0)

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
     paste("SSR = ",round(sum(lnk10_GA$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dGA,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_GA)[2],5),sep=""),adj=0)

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
     paste("SSR = ",round(sum(lnk10_CREA$residuals^2),5),sep=""),adj=0)
text(min(dfCOV$dCREA,na.rm=TRUE),0.61*max(dfCOV$eta.lK10,na.rm=TRUE),
     paste("slope = ",round(coef(lnk10_CREA)[2],5),sep=""),adj=0)

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
plot(dfCOV$dKOF, dfCOV$eta.lv ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="saemix V - BSA", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF)[1]+coef(lnv_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dKOF, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dKOF, dfCOV$eta.lv_corrKOF ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="saemix V - BSA - V corr with BSA", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF_corrKOF)[1]+coef(lnv_KOF_corrKOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF_corrKOF)[1]+coef(lnv_KOF_corrKOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dKOF, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dKOF, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(V) random effects", main ="saemix V - BSA - V corr with BSA, 
     k10 corr with AGE ", type='n')
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnv_KOF_corrKOFAGE)[1]+coef(lnv_KOF_corrKOFAGE)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnv_KOF_corrKOFAGE)[1]+coef(lnv_KOF_corrKOFAGE)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dKOF, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")

lnv_AGE <- lm(dfCOV$eta.lv ~ dfCOV$dAGE)
lnv_AGE_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dAGE)
lnv_AGE_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dAGE)
plot(dfCOV$dAGE, dfCOV$eta.lv ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="saemix V - AGE", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE)[1]+coef(lnv_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dAGE, dfCOV$eta.lv_corrKOF ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="saemix V - AGE - V corr with BSA", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE_corrKOF)[1]+coef(lnv_AGE_corrKOF)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE_corrKOF)[1]+coef(lnv_AGE_corrKOF)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dAGE, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dAGE, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(V) random effects", main ="saemix V - AGE - V corr with BSA, 
     k10 corr with AGE ", type='n')
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnv_AGE_corrKOFAGE)[1]+coef(lnv_AGE_corrKOFAGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnv_AGE_corrKOFAGE)[1]+coef(lnv_AGE_corrKOFAGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dAGE, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")

lnv_WT <- lm(dfCOV$eta.lv ~ dfCOV$dWT)
lnv_WT_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dWT)
lnv_WT_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dWT)
plot(dfCOV$dWT, dfCOV$eta.lv ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="saemix V - WT")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT)[1]+coef(lnv_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT)[1]+coef(lnv_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dWT, dfCOV$eta.lv_corrKOF ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="saemix V - WT - V corr with BSA")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT_corrKOF)[1]+coef(lnv_WT_corrKOF)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT_corrKOF)[1]+coef(lnv_WT_corrKOF)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dWT, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(V) random effects",  type="n", main ="saemix V - WT - V corr with BSA, 
     k10 corr with AGE")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnv_WT_corrKOFAGE)[1]+coef(lnv_WT_corrKOFAGE)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnv_WT_corrKOFAGE)[1]+coef(lnv_WT_corrKOFAGE)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dWT, dfCOV$eta.lv_corrKOFAGE, pch=21, bg="lightblue4")


lnv_HEIGHT <- lm(dfCOV$eta.lv ~ dfCOV$dHEIGHT)
lnv_HEIGHT_corrKOF <- lm(dfCOV$eta.lv_corrKOF ~ dfCOV$dHEIGHT)
lnv_HEIGHT_corrKOFAGE <- lm(dfCOV$eta.lv_corrKOFAGE ~ dfCOV$dHEIGHT)
plot(dfCOV$dHEIGHT, dfCOV$eta.lv ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="saemix V - HEIGHT")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT)[1]+coef(lnv_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dHEIGHT, dfCOV$eta.lv, pch=21, bg="grey80")
plot(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOF ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="saemix V - HEIGHT - V corr with BSA")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnv_HEIGHT_corrKOF)[1]+coef(lnv_HEIGHT_corrKOF)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnv_HEIGHT_corrKOF)[1]+coef(lnv_HEIGHT_corrKOF)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dHEIGHT, dfCOV$eta.lv_corrKOFAGE ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(V) random effects",  type="n", main ="saemix V - HEIGHT - V corr with BSA,
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
     main ="saemix K10 - BSA")
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF)[1]+coef(lnk10_KOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dKOF, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dKOF, dfCOV$eta.lK10_corrKOF ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - BSA - V corr with BSA")
points(c(min(dfCOV$dKOF,na.rm=TRUE),max(dfCOV$dKOF,na.rm=TRUE)),
       c(coef(lnk10_KOF_corrKOF)[1]+coef(lnk10_KOF_corrKOF)[2]*min(dfCOV$dKOF,na.rm=TRUE),
         coef(lnk10_KOF_corrKOF)[1]+coef(lnk10_KOF_corrKOF)[2]*max(dfCOV$dKOF,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dKOF, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dKOF, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(KOF) - mean(ln(KOF))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - BSA - V corr with BSA,
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
     main ="saemix K10 - AGE")
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE)[1]+coef(lnk10_AGE)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dAGE, dfCOV$eta.lK10_corrKOF ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - AGE - V corr with BSA")
points(c(min(dfCOV$dAGE,na.rm=TRUE),max(dfCOV$dAGE,na.rm=TRUE)),
       c(coef(lnk10_AGE_corrKOF)[1]+coef(lnk10_AGE_corrKOF)[2]*min(dfCOV$dAGE,na.rm=TRUE),
         coef(lnk10_AGE_corrKOF)[1]+coef(lnk10_AGE_corrKOF)[2]*max(dfCOV$dAGE,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dAGE, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue")
plot(dfCOV$dAGE, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(AGE) - mean(ln(AGE))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - AGE - V corr with BSA,
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
     main ="saemix K10 - WT")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT)[1]+coef(lnk10_WT)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dWT, dfCOV$eta.lK10_corrKOF ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - WT - V corr with BSA")
points(c(min(dfCOV$dWT,na.rm=TRUE),max(dfCOV$dWT,na.rm=TRUE)),
       c(coef(lnk10_WT_corrKOF)[1]+coef(lnk10_WT_corrKOF)[2]*min(dfCOV$dWT,na.rm=TRUE),
         coef(lnk10_WT_corrKOF)[1]+coef(lnk10_WT_corrKOF)[2]*max(dfCOV$dWT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dWT, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(WT) - mean(ln(WT))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - WT - V corr with BSA,
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
     main ="saemix K10 - HEIGHT")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT)[1]+coef(lnk10_HEIGHT)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(dfCOV$dWT, dfCOV$eta.lK10, pch=21, bg="grey80")
plot(dfCOV$dHEIGHT, dfCOV$eta.lK10_corrKOF ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - HEIGHT - V corr with BSA")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT_corrKOF)[1]+coef(lnk10_HEIGHT_corrKOF)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT_corrKOF)[1]+coef(lnk10_HEIGHT_corrKOF)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOF, pch=21, bg="lightblue")
plot(dfCOV$dHEIGHT, dfCOV$eta.lK10_corrKOFAGE ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="ln(k10) random effects",  type="n",
     main ="saemix K10 - HEIGHT - V corr with BSA,
     K10 corr with AGE")
points(c(min(dfCOV$dHEIGHT,na.rm=TRUE),max(dfCOV$dHEIGHT,na.rm=TRUE)),
       c(coef(lnk10_HEIGHT_corrKOFAGE)[1]+coef(lnk10_HEIGHT_corrKOFAGE)[2]*min(dfCOV$dHEIGHT,na.rm=TRUE),
         coef(lnk10_HEIGHT_corrKOFAGE)[1]+coef(lnk10_HEIGHT_corrKOFAGE)[2]*max(dfCOV$dHEIGHT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(dfCOV$dWT, dfCOV$eta.lK10_corrKOFAGE, pch=21, bg="lightblue4")

#-----------------------------------------------------------------------------------------------------------------------
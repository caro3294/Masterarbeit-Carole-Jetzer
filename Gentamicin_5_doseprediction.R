# Gentamicin Dose Prediction
# autor: Carole Jetzer
# next run Gentamicin_6_Patient_Characterisics.R

rm(list=ls())
Corr_data <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\Gentamicin_CORRdata.csv", sep=",")
dfCOV <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\Covariates_dataframe.csv", sep=",")
nlmixr_data <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\nlmixr_dataframe_results.csv", sep=",")
nlme_data <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\nlme_dataframe_results.csv", sep=",")
saemix_data <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\saemix_dataframe_results.csv", sep=",")

Genta.Data <- data.frame(
                  #Patient parameter:
                  ID=Corr_data$ï..Studien.ID, #Patient number 
                  AGE=as.numeric(as.character(Corr_data$AgeDose0)),      # [d]
                  BSA=as.numeric(as.character(Corr_data$KOF.berechnet)), # [m^2] Body surface area
                  HEIGHT=as.numeric(as.character(Corr_data$GROESSE)),    # [cm]
                  WEIGHT=as.numeric(as.character(Corr_data$GEWICHT)),    # [kg]
                  CREA=as.numeric(as.character(Corr_data$KREATININ)),    # [umol/L]
                  UREA=as.numeric(as.character(Corr_data$HARNSTOFF)),    # [mmol/L]
                  GA=as.numeric(as.character(Corr_data$GA..Weeks.)),     # [weeks] week of gestation
                  SEX=Corr_data$Sex, 
                  DIAGNOSE=Corr_data$DIAGNOSE_GROUPED,
                  #Dosing parameter       
                  INF1=Corr_data$TimeDose1*24, #Time point of infusion [h]
                  INF2=Corr_data$TimeDose2*24,
                  INF3=Corr_data$TimeDose3*24,
                  TIME_I1=Corr_data$APPLIKATIONSDAUER/60,   #Duration of the infusion 1 [h]
                  TIME_I2=Corr_data$APPLIKATIONSDAUER.2/60, #Durarion of the infusion 2
                  TIME_I3=Corr_data$APPLIKATIONSDAUER.3/60, #Durarion of the infusion 3
                  DOSE1=Corr_data$DOSIERUNG.1, #Dose 1 [mg]
                  DOSE2=Corr_data$DOSIERUNG.2, #Dose 2
                  DOSE3=Corr_data$DOSIERUNG.3, #Dose 3
                  TIME_M1=Corr_data$TimeGenta1*24, # Time of the measurement 1 [h]
                  TIME_M2=Corr_data$TimeGenta2*24, # Time of the measurement 2
                  TIME_M3=Corr_data$TimeGenta3*24, # Time of the measurement 3
                  CONCM1=Corr_data$Gentamicin.Spiegel.nach.30.min, #measured conc 1 [mg/L]
                  CONCM2=Corr_data$Gentamicin.Spiegel.nach.4h, #measured conc 1 [mg/L]
                  CONCM3=Corr_data$Gentamicin.Spiegel.nach.24h #measured conc 1 [mg/L]
                  )

Genta.Data$CREA[Genta.Data$CREA==0] <- NA
Genta.Data$UREA[Genta.Data$UREA==0] <- NA

Genta.Results <- data.frame(ID=Corr_data$ï..Studien.ID,
                            nlme_K10_corr = nlme_data$nlme.Corr_k10,
                            nlme_K10_corrETA = nlme_data$nlme.CorrETA_k10,
                            nlme_V_corr = nlme_data$nlme.Corr_V,
                            nlme_V_corrETA = nlme_data$nlme.CorrETA_V,
                            saemix_K10_corr = saemix_data$saemix.Corr_k10,
                            saemix_K10_corrETA = saemix_data$saemix.CorrETA_k10,
                            saemix_V_corr = saemix_data$saemix.Corr_V,
                            saemix_V_corrETA = saemix_data$saemix.CorrETA_V,
                            nlmixr_K10_corrETA = nlmixr_data$nlmixr.CorrETA_k10,
                            nlmixr_V_corrETA = nlmixr_data$nlmixr.CorrETA_V
                            )

write.csv(Genta.Results,"C:/Users/carol/BioMed IT CSV/CSV\\Genta_results.csv")

#Dose:
c_max <- 12 #12mg/L
Tinf <- 2 #min
#Poster Formel: Di = Cmax × ki × Vi × Ti / [1-exp(-ki × Ti)]


Genta.Results$Dose2minPred_sameix <- c_max * (Genta.Results$saemix_K10_corr) * (Genta.Results$saemix_V_corr) * Tinf /
  (1-exp(-(Genta.Results$saemix_K10_corr) * Tinf))/1000
Genta.Results$Dose2minPred_sameixETA <- c_max * (Genta.Results$saemix_K10_corrETA) * (Genta.Results$saemix_V_corrETA) * Tinf /
  (1-exp(-(Genta.Results$saemix_K10_corrETA) * Tinf))/1000

Genta.Results$Dose2minPred_nlmixr <- c_max * (Genta.Results$nlmixr_K10_corr) * (Genta.Results$nlmixr_V_corr) * Tinf /
  (1-exp(-(Genta.Results$nlmixr_K10_corr) * Tinf)) 
Genta.Results$Dose2minPred_nlmixrETA <- c_max * (Genta.Results$nlmixr_K10_corrETA) * (Genta.Results$nlmixr_V_corrETA) * Tinf /
  (1-exp(-(Genta.Results$nlmixr_K10_corrETA) * Tinf)) 

Genta.Results$Dose2minPred_nlme <- c_max * (Genta.Results$nlme_K10_corr) * (Genta.Results$nlme_V_corr) * Tinf /
  (1-exp(-(Genta.Results$nlme_K10_corr) * Tinf))/1000
Genta.Results$Dose2minPred_nlmeETA <- c_max * (Genta.Results$nlme_K10_corrETA) * (Genta.Results$nlme_V_corrETA) * Tinf /
  (1-exp(-(Genta.Results$nlme_K10_corrETA) * Tinf))/1000


#Predicted vs. gived dose (WEIGHT)
lnDose_Weight <- lm(Genta.Data$DOSE1 ~ Genta.Data$WEIGHT)
lnDose_Weight_saemix <- lm(Genta.Results$Dose2minPred_sameix ~ Genta.Data$WEIGHT)
lnDose_Weight_saemixETA <- lm(Genta.Results$Dose2minPred_sameixETA ~ Genta.Data$WEIGHT)
lnDose_Weight_nlme <- lm(Genta.Results$Dose2minPred_nlme ~ Genta.Data$WEIGHT)
lnDose_Weight_nlmeETA <- lm(Genta.Results$Dose2minPred_nlmeETA ~ Genta.Data$WEIGHT)
lnDose_Weight_nlmixr <- lm(Genta.Results$Dose2minPred_nlmixr ~ Genta.Data$WEIGHT)
lnDose_Weight_nlmixrETA <- lm(Genta.Results$Dose2minPred_nlmixrETA ~ Genta.Data$WEIGHT)

par(mfrow=c(1,1))
plot(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_sameix ,xlab="Weight [kg]",ylab="Dose [mg]",  type="n",
     main ="Predicted and given doses vs. body weight - sameix fit")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_saemixETA)[1]+coef(lnDose_Weight_saemixETA)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_saemixETA)[1]+coef(lnDose_Weight_saemixETA)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_saemix)[1]+coef(lnDose_Weight_saemix)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_saemix)[1]+coef(lnDose_Weight_saemix)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(Genta.Data$WEIGHT, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_sameixETA, pch=21, bg="lightblue")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_sameix, pch=21, bg="lightblue4")
legend("topleft",  legend=c("given dose", "predicted dose with ETAs", "predicted dose without ETAs"), pch=21, 
       pt.bg=c("grey80","lightblue", "lightblue4"))

text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.86*max(Genta.Results$Dose2minPred_sameix,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_saemixETA)[1],4),sep=""),adj=0, col="lightblue")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.83*max(Genta.Results$Dose2minPred_sameix,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_saemixETA)[2],4),sep=""),adj=0, col="lightblue")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.80*max(Genta.Results$Dose2minPred_sameix,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_saemix)[1],4),sep=""),adj=0, col="lightblue4")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.77*max(Genta.Results$Dose2minPred_sameix,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_saemix)[2],4),sep=""),adj=0, col="lightblue4")


plot(Genta.Data$WEIGHT,  Genta.Results$Dose2minPred_nlme ,xlab="Weight [kg]",ylab="Dose [mg]",  type="n",
     main ="Predicted and given doses vs. body weight - nlme fit")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_nlmeETA)[1]+coef(lnDose_Weight_nlmeETA)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_nlmeETA)[1]+coef(lnDose_Weight_nlmeETA)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_nlme)[1]+coef(lnDose_Weight_nlme)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_nlme)[1]+coef(lnDose_Weight_nlme)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(Genta.Data$WEIGHT, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlmeETA, pch=21, bg="lightblue")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlme, pch=21, bg="lightblue4")
legend("topleft",  legend=c("given dose", "predicted dose with ETAs", "predicted dose without ETAs"), pch=21, 
       pt.bg=c("grey80","lightblue", "lightblue4"))
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.86*max(Genta.Results$Dose2minPred_nlme,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_nlmeETA)[1],4),sep=""),adj=0, col="lightblue")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.83*max(Genta.Results$Dose2minPred_nlme,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_nlmeETA)[2],4),sep=""),adj=0, col="lightblue")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.80*max(Genta.Results$Dose2minPred_nlme,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_nlme)[1],4),sep=""),adj=0, col="lightblue4")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.77*max(Genta.Results$Dose2minPred_nlme,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_nlme)[2],4),sep=""),adj=0, col="lightblue4")

plot(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlmixr ,xlab="Weight [kg]",ylab="Dose [mg]",  type="n",
     main ="Predicted and given doses vs. body weight - nlmixr fit")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_nlmixrETA)[1]+coef(lnDose_Weight_nlmixrETA)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_nlmixrETA)[1]+coef(lnDose_Weight_nlmixrETA)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_nlmixr)[1]+coef(lnDose_Weight_nlmixr)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_nlmixr)[1]+coef(lnDose_Weight_nlmixr)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue4")
points(Genta.Data$WEIGHT, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlmixrETA, pch=21, bg="lightblue")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlmixr, pch=21, bg="lightblue4")
legend("topleft",  legend=c("given dose", "predicted dose with ETAs", "predicted dose without ETAs"), pch=21, 
       pt.bg=c("grey80","lightblue", "lightblue4"))
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.86*max(Genta.Results$Dose2minPred_nlmixr,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_nlmixrETA)[1],4),sep=""),adj=0, col="lightblue")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.83*max(Genta.Results$Dose2minPred_nlmixr,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_nlmixrETA)[2],4),sep=""),adj=0, col="lightblue")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.80*max(Genta.Results$Dose2minPred_nlmixr,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_nlmixr)[1],4),sep=""),adj=0, col="lightblue4")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.77*max(Genta.Results$Dose2minPred_nlmixr,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_nlmixr)[2],4),sep=""),adj=0, col="lightblue4")

#Predicted vs. gived dose (BSA)
lnDose_BSA <- lm(Genta.Data$DOSE1 ~ Genta.Data$BSA)
lnDose_BSA_saemix <- lm(Genta.Results$Dose2minPred_sameix ~ Genta.Data$BSA)
lnDose_BSA_saemixETA <- lm(Genta.Results$Dose2minPred_sameixETA ~ Genta.Data$BSA)
lnDose_BSA_nlme <- lm(Genta.Results$Dose2minPred_nlme ~ Genta.Data$BSA)
lnDose_BSA_nlmeETA <- lm(Genta.Results$Dose2minPred_nlmeETA ~ Genta.Data$BSA)
lnDose_BSA_nlmixr <- lm(Genta.Results$Dose2minPred_nlmixr ~ Genta.Data$BSA)
lnDose_BSA_nlmixrETA <- lm(Genta.Results$Dose2minPred_nlmixrETA ~ Genta.Data$BSA)

par(mfrow=c(1,1))
plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2] ",ylab="Dose [mg]",  type="n",
     main ="Predicted and given doses vs. BSA - sameix fit")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_saemix)[1]+coef(lnDose_BSA_saemix)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_saemix)[1]+coef(lnDose_BSA_saemix)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="red")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_sameix, pch=21, bg="red")
legend("topleft",  legend=c("given dose", "predicted dose"), pch=21, pt.bg=c("grey80","red"))
text(min(Genta.Data$BSA,na.rm=TRUE),0.90*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_saemix)[1],4),sep=""),adj=0, col="red")
text(min(Genta.Data$BSA,na.rm=TRUE),0.87*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_saemix)[2],4),sep=""),adj=0, col="red")

plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2]",ylab="Dose [mg]",  type="n",
     main ="Predicted and given doses vs. BSA - nlme fit")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_nlme)[1]+coef(lnDose_BSA_nlme)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_nlme)[1]+coef(lnDose_BSA_nlme)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="red")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_nlme, pch=21, bg="red")
legend("topleft",  legend=c("given dose", "predicted dose"), pch=21, pt.bg=c("grey80","red"))
text(min(Genta.Data$BSA,na.rm=TRUE),0.90*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_nlme)[1],4),sep=""),adj=0, col="red")
text(min(Genta.Data$BSA,na.rm=TRUE),0.87*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_nlme)[2],4),sep=""),adj=0, col="red")

plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2]",ylab="Dose [mg]",  type="n",
     main ="Predicted and given doses vs. BSA - nlmixr fit")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_nlmixr)[1]+coef(lnDose_BSA_nlmixr)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_nlmixr)[1]+coef(lnDose_BSA_nlmixr)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="red")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_sameix, pch=21, bg="red")
legend("topleft",  legend=c("given dose", "predicted dose"), pch=21, pt.bg=c("grey80","red"))
text(min(Genta.Data$BSA,na.rm=TRUE),0.90*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_nlmixr)[1],4),sep=""),adj=0, col="red")
text(min(Genta.Data$BSA,na.rm=TRUE),0.87*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_nlmixr)[2],4),sep=""),adj=0, col="red")

#----------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,3))
plot(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_sameix ,xlab="Weight [kg]",ylab="Dose [mg]",  type="n",
     main ="Dose vs. Weight - saemix Prediction", ylim=c(0,60))
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_saemixETA)[1]+coef(lnDose_Weight_saemixETA)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_saemixETA)[1]+coef(lnDose_Weight_saemixETA)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$WEIGHT, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_sameixETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - saemix fit"), pch=21, 
       pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.55*max(Genta.Results$Dose2minPred_sameix,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_saemixETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.52*max(Genta.Results$Dose2minPred_sameix,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_saemixETA)[2],4),sep=""),adj=0, col="lightblue3")


plot(Genta.Data$WEIGHT,  Genta.Results$Dose2minPred_nlme ,xlab="Weight [kg]",ylab="Dose [mg]",  type="n",
     main ="Dose vs. Weight - nlme Prediction", ylim=c(0,60))
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_nlmeETA)[1]+coef(lnDose_Weight_nlmeETA)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_nlmeETA)[1]+coef(lnDose_Weight_nlmeETA)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$WEIGHT, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlmeETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - nlme fit"), pch=21, 
       pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.83*max(Genta.Results$Dose2minPred_nlme,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_nlmeETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),0.78*max(Genta.Results$Dose2minPred_nlme,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_nlmeETA)[2],4),sep=""),adj=0, col="lightblue3")


plot(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlmixr ,xlab="Weight [kg]",ylab="Dose [mg]",  type="n",
     main ="Dose vs. Weight - nlmixr Prediction", ylim=c(0,60))
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight)[1]+coef(lnDose_Weight)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$WEIGHT,na.rm=TRUE),max(Genta.Data$WEIGHT,na.rm=TRUE)),
       c(coef(lnDose_Weight_nlmixrETA)[1]+coef(lnDose_Weight_nlmixrETA)[2]*min(Genta.Data$WEIGHT,na.rm=TRUE),
         coef(lnDose_Weight_nlmixrETA)[1]+coef(lnDose_Weight_nlmixrETA)[2]*max(Genta.Data$WEIGHT,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$WEIGHT, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$WEIGHT, Genta.Results$Dose2minPred_nlmixrETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - nlmixr fit"), pch=21, 
       pt.bg=c("grey80","lightblue", "lightblue4"))
text(min(Genta.Data$WEIGHT,na.rm=TRUE),1.12*max(Genta.Results$Dose2minPred_nlmixr,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_Weight_nlmixrETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$WEIGHT,na.rm=TRUE),1.05*max(Genta.Results$Dose2minPred_nlmixr,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_Weight_nlmixrETA)[2],4),sep=""),adj=0, col="lightblue3")

#Predicted vs. gived dose (BSA)
lnDose_BSA <- lm(Genta.Data$DOSE1 ~ Genta.Data$BSA)
lnDose_BSA_saemix <- lm(Genta.Results$Dose2minPred_sameix ~ Genta.Data$BSA)
lnDose_BSA_saemixETA <- lm(Genta.Results$Dose2minPred_sameixETA ~ Genta.Data$BSA)
lnDose_BSA_nlme <- lm(Genta.Results$Dose2minPred_nlme ~ Genta.Data$BSA)
lnDose_BSA_nlmeETA <- lm(Genta.Results$Dose2minPred_nlmeETA ~ Genta.Data$BSA)
lnDose_BSA_nlmixr <- lm(Genta.Results$Dose2minPred_nlmixr ~ Genta.Data$BSA)
lnDose_BSA_nlmixrETA <- lm(Genta.Results$Dose2minPred_nlmixrETA ~ Genta.Data$BSA)

par(mfrow=c(1,3))
plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2] ",ylab="Dose [mg]",  type="n",
     main ="Dose vs. BSA - sameix Prediction")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_saemix)[1]+coef(lnDose_BSA_saemix)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_saemix)[1]+coef(lnDose_BSA_saemix)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="red")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_sameix, pch=21, bg="red")
legend("topleft",  legend=c("given dose", "predicted dose"), pch=21, pt.bg=c("grey80","red"))
text(min(Genta.Data$BSA,na.rm=TRUE),0.90*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_saemix)[1],4),sep=""),adj=0, col="red")
text(min(Genta.Data$BSA,na.rm=TRUE),0.87*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_saemix)[2],4),sep=""),adj=0, col="red")

plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2]",ylab="Dose [mg]",  type="n",
     main ="Dose vs. BSA - nlme Prediction")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_nlme)[1]+coef(lnDose_BSA_nlme)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_nlme)[1]+coef(lnDose_BSA_nlme)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="red")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_nlme, pch=21, bg="red")
legend("topleft",  legend=c("given dose", "predicted dose"), pch=21, pt.bg=c("grey80","red"))
text(min(Genta.Data$BSA,na.rm=TRUE),0.90*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_nlme)[1],4),sep=""),adj=0, col="red")
text(min(Genta.Data$BSA,na.rm=TRUE),0.87*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_nlme)[2],4),sep=""),adj=0, col="red")

plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2]",ylab="Dose [mg]",  type="n",
     main ="Dose vs. BSA - nlmixr Prediction")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_nlmixr)[1]+coef(lnDose_BSA_nlmixr)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_nlmixr)[1]+coef(lnDose_BSA_nlmixr)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="red")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_sameix, pch=21, bg="red")
legend("topleft",  legend=c("given dose", "predicted dose"), pch=21, pt.bg=c("grey80","red"))
text(min(Genta.Data$BSA,na.rm=TRUE),0.90*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_nlmixr)[1],4),sep=""),adj=0, col="red")
text(min(Genta.Data$BSA,na.rm=TRUE),0.87*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_nlmixr)[2],4),sep=""),adj=0, col="red")

#----------------------------------------------------------------------------------------------------------------
#Predicted vs. gived dose (BSA)
lnDose_BSA <- lm(Genta.Data$DOSE1 ~ Genta.Data$BSA)
lnDose_BSA_saemix <- lm(Genta.Results$Dose2minPred_sameix ~ Genta.Data$BSA)
lnDose_BSA_saemixETA <- lm(Genta.Results$Dose2minPred_sameixETA ~ Genta.Data$BSA)
lnDose_BSA_nlme <- lm(Genta.Results$Dose2minPred_nlme ~ Genta.Data$BSA)
lnDose_BSA_nlmeETA <- lm(Genta.Results$Dose2minPred_nlmeETA ~ Genta.Data$BSA)
lnDose_BSA_nlmixr <- lm(Genta.Results$Dose2minPred_nlmixr ~ Genta.Data$BSA)
lnDose_BSA_nlmixrETA <- lm(Genta.Results$Dose2minPred_nlmixrETA ~ Genta.Data$BSA)

par(mfrow=c(1,3))
plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2] ",ylab="Dose [mg]",  type="n",
     main ="Dose vs. BSA - saemix Prediction", ylim=c(0,60))
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_saemixETA)[1]+coef(lnDose_BSA_saemixETA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_saemixETA)[1]+coef(lnDose_BSA_saemixETA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_sameixETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - saemix fit"), pch=21, pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$BSA,na.rm=TRUE),1.1*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_saemixETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$BSA,na.rm=TRUE),1.03*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_saemixETA)[2],4),sep=""),adj=0, col="lightblue3")

plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2] ",ylab="Dose [mg]",  type="n",
     main ="Dose vs. BSA - nlme Prediction", ylim=c(0,60))
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_nlmeETA)[1]+coef(lnDose_BSA_nlmeETA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_nlmeETA)[1]+coef(lnDose_BSA_nlmeETA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_nlmeETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - nlme fit"), pch=21, pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$BSA,na.rm=TRUE),1.1*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_nlmeETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$BSA,na.rm=TRUE),1.03*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_nlmeETA)[2],4),sep=""),adj=0, col="lightblue3")

plot(Genta.Data$BSA, Genta.Data$DOSE1 ,xlab="BSA [m^2] ",ylab="Dose [mg]",  type="n",
     main ="Dose vs. BSA - nlmixr Prediction", ylim=c(0,60))
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA)[1]+coef(lnDose_BSA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$BSA,na.rm=TRUE),max(Genta.Data$BSA,na.rm=TRUE)),
       c(coef(lnDose_BSA_nlmixrETA)[1]+coef(lnDose_BSA_nlmixrETA)[2]*min(Genta.Data$BSA,na.rm=TRUE),
         coef(lnDose_BSA_nlmixrETA)[1]+coef(lnDose_BSA_nlmixrETA)[2]*max(Genta.Data$BSA,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$BSA, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$BSA, Genta.Results$Dose2minPred_nlmixrETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - nlmixr fit"), pch=21, pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$BSA,na.rm=TRUE),1.1*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_BSA_nlmixrETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$BSA,na.rm=TRUE),1.03*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_BSA_nlmixrETA)[2],4),sep=""),adj=0, col="lightblue3")

#-----------------------------------------------------------------------------------------------------------------------
#Predicted vs. gived dose (AGE)
lnDose_AGE <- lm(Genta.Data$DOSE1 ~ Genta.Data$AGE)
lnDose_AGE_saemix <- lm(Genta.Results$Dose2minPred_sameix ~ Genta.Data$AGE)
lnDose_AGE_saemixETA <- lm(Genta.Results$Dose2minPred_sameixETA ~ Genta.Data$AGE)
lnDose_AGE_nlme <- lm(Genta.Results$Dose2minPred_nlme ~ Genta.Data$AGE)
lnDose_AGE_nlmeETA <- lm(Genta.Results$Dose2minPred_nlmeETA ~ Genta.Data$AGE)
lnDose_AGE_nlmixrETA <- lm(Genta.Results$Dose2minPred_nlmixrETA ~ Genta.Data$AGE)

par(mfrow=c(1,3))
plot(Genta.Data$AGE, Genta.Data$DOSE1 ,xlab="AGE [d] ",ylab="Dose [mg]",  type="n",
     main ="Dose vs. AGE - saemix Prediction", ylim=c(0,60))
points(c(min(Genta.Data$AGE,na.rm=TRUE),max(Genta.Data$AGE,na.rm=TRUE)),
       c(coef(lnDose_AGE)[1]+coef(lnDose_AGE)[2]*min(Genta.Data$AGE,na.rm=TRUE),
         coef(lnDose_AGE)[1]+coef(lnDose_AGE)[2]*max(Genta.Data$AGE,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$AGE,na.rm=TRUE),max(Genta.Data$AGE,na.rm=TRUE)),
       c(coef(lnDose_AGE_saemixETA)[1]+coef(lnDose_AGE_saemixETA)[2]*min(Genta.Data$AGE,na.rm=TRUE),
         coef(lnDose_AGE_saemixETA)[1]+coef(lnDose_AGE_saemixETA)[2]*max(Genta.Data$AGE,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$AGE, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$AGE, Genta.Results$Dose2minPred_sameixETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - saemix fit"), pch=21, pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$AGE,na.rm=TRUE),1.1*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_AGE_saemixETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$AGE,na.rm=TRUE),1.03*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_AGE_saemixETA)[2],4),sep=""),adj=0, col="lightblue3")

plot(Genta.Data$AGE, Genta.Data$DOSE1 ,xlab="AGE [d] ",ylab="Dose [mg]",  type="n",
     main ="Dose vs. AGE - nlme Prediction", ylim=c(0,60))
points(c(min(Genta.Data$AGE,na.rm=TRUE),max(Genta.Data$AGE,na.rm=TRUE)),
       c(coef(lnDose_AGE)[1]+coef(lnDose_AGE)[2]*min(Genta.Data$AGE,na.rm=TRUE),
         coef(lnDose_AGE)[1]+coef(lnDose_AGE)[2]*max(Genta.Data$AGE,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$AGE,na.rm=TRUE),max(Genta.Data$AGE,na.rm=TRUE)),
       c(coef(lnDose_AGE_nlmeETA)[1]+coef(lnDose_AGE_nlmeETA)[2]*min(Genta.Data$AGE,na.rm=TRUE),
         coef(lnDose_AGE_nlmeETA)[1]+coef(lnDose_AGE_nlmeETA)[2]*max(Genta.Data$AGE,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$AGE, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$AGE, Genta.Results$Dose2minPred_nlmeETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - nlme fit"), pch=21, pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$AGE,na.rm=TRUE),1.1*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_AGE_nlmeETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$AGE,na.rm=TRUE),1.03*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_AGE_nlmeETA)[2],4),sep=""),adj=0, col="lightblue3")

plot(Genta.Data$AGE, Genta.Data$DOSE1 ,xlab="AGE [d] ",ylab="Dose [mg]",  type="n",
     main ="Dose vs. AGE - nlmixr Prediction", ylim=c(0,60))
points(c(min(Genta.Data$AGE,na.rm=TRUE),max(Genta.Data$AGE,na.rm=TRUE)),
       c(coef(lnDose_AGE)[1]+coef(lnDose_AGE)[2]*min(Genta.Data$AGE,na.rm=TRUE),
         coef(lnDose_AGE)[1]+coef(lnDose_AGE)[2]*max(Genta.Data$AGE,na.rm=TRUE)),
       type='l', col="grey60")
points(c(min(Genta.Data$AGE,na.rm=TRUE),max(Genta.Data$AGE,na.rm=TRUE)),
       c(coef(lnDose_AGE_nlmixrETA)[1]+coef(lnDose_AGE_nlmixrETA)[2]*min(Genta.Data$AGE,na.rm=TRUE),
         coef(lnDose_AGE_nlmixrETA)[1]+coef(lnDose_AGE_nlmixrETA)[2]*max(Genta.Data$AGE,na.rm=TRUE)),
       type='l', col="lightblue")
points(Genta.Data$AGE, Genta.Data$DOSE1, pch=21, bg="grey80")
points(Genta.Data$AGE, Genta.Results$Dose2minPred_nlmixrETA, pch=21, bg="lightblue")
legend("topleft",  legend=c("Given Dose", "Predicted Dose - nlmixr fit"), pch=21, pt.bg=c("grey80","lightblue"))
text(min(Genta.Data$AGE,na.rm=TRUE),1.1*max(Genta.Data$DOSE1,na.rm=TRUE), 
     paste("intercept = ",round(coef(lnDose_AGE_nlmixrETA)[1],4),sep=""),adj=0, col="lightblue3")
text(min(Genta.Data$AGE,na.rm=TRUE),1.03*max(Genta.Data$DOSE1,na.rm=TRUE),
     paste("slope = ",round(coef(lnDose_AGE_nlmixrETA)[2],4),sep=""),adj=0, col="lightblue3")




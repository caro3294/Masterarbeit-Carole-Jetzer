# Voriconazole Covariates
# autor: Carole Jetzer

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

library(hms)
if(!("hms" %in% rownames(installed.packages()))){
  install.packages("hms")
}

library(scatterplot3d)
if(!("scatterplot3d" %in% rownames(installed.packages()))){
  install.packages("scatterplot3d")
}
library(ggplot2)
if(!("ggplot2" %in% rownames(installed.packages()))){
  install.packages("ggplot2")
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


rm(list=ls()) 
par.default <- par(no.readonly = TRUE) 
setwd("C:/Users/carol/OneDrive/Studium/ETH/Master/Masterarbeit") 

Corr_data <- read.csv("D:/BioMed IT CSV/CSV_Voriconazol\\Voriconazol_CORRdata.csv", sep=",")
dfCOV <- read.csv("D:/BioMed IT CSV/CSV_Voriconazol\\Voriconazole_Covariates.csv", sep=",")
df.Voridata <- read.csv("D:/BioMed IT CSV/CSV_Voriconazol\\Voriconazol_data2.csv", sep=",")

Corr.Voridata <- data.frame(ID=Corr_data$ID,
                       CONC=Corr_data$CONC,
                       TIME_M=Corr_data$TIME_M,
                       TIME_A=Corr_data$TIME_A,
                       ROA=Corr_data$ROA,
                       DOSE=Corr_data$DOSE,
                       INT=dfCOV$INT,
                       AGE=dfCOV$AGE,
                       BSA=dfCOV$BSA,
                       HEIGHT=dfCOV$HEIGHT,
                       WT=dfCOV$WT,
                       lAGE=dfCOV$lAGE,
                       lHEIGHT=dfCOV$lHEIGHT,
                       lWT=dfCOV$lWT,
                       dAGE=dfCOV$dAGE,
                       dBSA=dfCOV$dBSA,
                       dHEIGHT=dfCOV$dHEIGHT,
                       dWT=dfCOV$dWT
                       )

Voridata <- df.Voridata

df.Voridata$AgeYears <- as.numeric(as.character(df.Voridata$AGE/365))
MeanAge <- mean(df.Voridata$AgeYears)
sdAge <- sqrt(var(df.Voridata$AgeYears))

par(mfrow=c(1,1))
hist(df.Voridata$AgeYears,
     main="Age distribution",
     xlab="Age [years]",
     xlim=c(0,18),
     ylim=c(0,60),
     breaks=20,
     border="lightblue4", 
     col="lightblue",
     freq=TRUE,
     xaxt='n'
)
axis(side=1, at=seq(0,18 ,3), labels=seq(0,18, 3))
text(2.5,57, "mean Age = 9.750 years ", cex=1)
text(2.2,47, "sd Age = 5.851 ", cex=1)

#polygon for therapeutic window (3-6)
y <- c(5.5,5.5,1,1)   # y-Werte der Eckpunkte
x <- c(-3,2,2,-3)   # x-Werte der Eckpunkte

#-------------------------------------------------------------------------------------------------------------
#Route of administration (ROA), 1 = i.v. , 2 = p.o. , 3 = i.v.+p.o.
par(mfrow=c(1,2))
ln_AGE <- lm(Voridata$CONC ~ Voridata$dAGE)
plot(Voridata$dAGE, Voridata$CONC ,xlab="ln(AGE) - mean(ln(AGE))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="AGE", yaxt='n', ylim=c(0,10))
polygon(x,y, angle=45, col="grey90")
points(Voridata$dAGE[Voridata$ROA==1], Voridata$CONC[Voridata$ROA==1], pch=21, bg="lightblue")
points(Voridata$dAGE[Voridata$ROA==2], Voridata$CONC[Voridata$ROA==2], pch=21, bg="lightblue4")
points(Voridata$dAGE[Voridata$ROA==3], Voridata$CONC[Voridata$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(0,10 ,2), labels=seq(0,10, 2))
axis(side=4, at=seq(-1,12 ,12), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

ln_BSA <- lm(Voridata$CONC ~ Voridata$dBSA)
plot(Voridata$dBSA, Voridata$CONC ,xlab="ln(BSA) - mean(ln(BSA))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="BSA", yaxt='n', ylim=c(0,10))
polygon(x,y, angle=45, col="grey90")
points(Voridata$dBSA[Voridata$ROA==1], Voridata$CONC[Voridata$ROA==1], pch=21, bg="lightblue")
points(Voridata$dBSA[Voridata$ROA==2], Voridata$CONC[Voridata$ROA==2], pch=21, bg="lightblue4")
points(Voridata$dBSA[Voridata$ROA==3], Voridata$CONC[Voridata$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(0,10 ,2), labels=seq(0,10, 2))
axis(side=4, at=seq(-1,12 ,12), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

ln_WT <- lm(Voridata$CONC ~ Voridata$dWT)
plot(Voridata$dWT, Voridata$CONC ,xlab="ln(WT) - mean(ln(WT))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="WEIGHT" , yaxt='n', ylim=c(0,10))
polygon(x,y, angle=45, col="grey90")
points(Voridata$dWT[Voridata$ROA==1], Voridata$CONC[Voridata$ROA==1], pch=21, bg="lightblue")
points(Voridata$dWT[Voridata$ROA==2], Voridata$CONC[Voridata$ROA==2], pch=21, bg="lightblue4")
points(Voridata$dWT[Voridata$ROA==3], Voridata$CONC[Voridata$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(0,10 ,2), labels=seq(0,10, 2))
axis(side=4, at=seq(-1,12 ,12), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

ln_HEIGHT <- lm(Voridata$CONC ~ Voridata$dHEIGHT)
plot(Voridata$dHEIGHT, Voridata$CONC ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="HEIGHT" , yaxt='n', ylim=c(0,10))
polygon(x,y, angle=45, col="grey90")
points(Voridata$dHEIGHT[Voridata$ROA==1], Voridata$CONC[Voridata$ROA==1], pch=21, bg="lightblue")
points(Voridata$dHEIGHT[Voridata$ROA==2], Voridata$CONC[Voridata$ROA==2], pch=21, bg="lightblue4")
points(Voridata$dHEIGHT[Voridata$ROA==3], Voridata$CONC[Voridata$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(0,10 ,2), labels=seq(0,10, 2))
axis(side=4, at=seq(-1,12 ,12), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

#-------------------------------------------------------------------------------------------------------------
#Interaction
# 1 = CYP2C19-Hemmer (Ome)
# 2 = CYP2C19-Hemmer (Ome + Eso) + CYP3A4-Induktoren (Gluco)
# 4 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Hemmer (Ery)
# 5 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Induktoren (Gluco)
# 7 = CYP2C9-Hemmer (Bac)

Voridata <- Corr.Voridata

plot(Voridata$dAGE, Voridata$CONC ,xlab="ln(AGE) - mean(ln(AGE))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Voriconazole measured concentration vs. AGE")
polygon(x,y, angle=45, col="grey90")
points(Voridata$dAGE[is.na(Voridata$INT)], Voridata$CONC[is.na(Voridata$INT)], pch=21, bg="grey50")
points(Voridata$dAGE[is.na(dfCOV$CO.Medi)], Voridata$CONC[is.na(dfCOV$CO.Medi)], pch=21, bg="lightblue")
points(Voridata$dAGE[Voridata$INT==1], Voridata$CONC[Voridata$INT==1], pch=21, bg="green")
points(Voridata$dAGE[Voridata$INT==2], Voridata$CONC[Voridata$INT==2], pch=21, bg="blue")
points(Voridata$dAGE[Voridata$INT==4], Voridata$CONC[Voridata$INT==4], pch=21, bg="maroon")
points(Voridata$dAGE[Voridata$INT==5], Voridata$CONC[Voridata$INT==5], pch=21, bg="red")
points(Voridata$dAGE[Voridata$INT==7], Voridata$CONC[Voridata$INT==7], pch=21, bg="orchid1")
legend("topleft", legend=c( "no Co-Medication ", "no Interaction", 
"CYP2C19-Hemmer (Omeprazol)", 
"CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP3A4-Induktoren (Glucocorticoide)",
"CYP2C19-Hemmer (Omeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Hemmer (Erythromycin)",
"CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Induktoren (Glucocorticoide)",
"CYP2C9-Hemmer (Co-Trimazol)"),
       pch=c(21,21,21,21,21,21,21,21), pt.bg=c("lightblue","grey50", "maroon", "blue", "red", "chartreuse1", "orchid1"), cex=0.65)

plot(Voridata$dBSA, Voridata$CONC ,xlab="ln(BSA) - mean(ln(BSA))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Voriconazole measured concentration vs. BSA")
polygon(x,y, angle=45, col="grey90")
points(Voridata$dBSA[is.na(Voridata$INT)], Voridata$CONC[is.na(Voridata$INT)], pch=21, bg="grey50")
points(Voridata$dBSA[is.na(dfCOV$CO.Medi)], Voridata$CONC[is.na(dfCOV$CO.Medi)], pch=21, bg="lightblue")
points(Voridata$dBSA[Voridata$INT==1], Voridata$CONC[Voridata$INT==1], pch=21, bg="green")
points(Voridata$dBSA[Voridata$INT==2], Voridata$CONC[Voridata$INT==2], pch=21, bg="blue")
points(Voridata$dAGE[Voridata$INT==4], Voridata$CONC[Voridata$INT==4], pch=21, bg="maroon")
points(Voridata$dBSA[Voridata$INT==5], Voridata$CONC[Voridata$INT==5], pch=21, bg="red")
points(Voridata$dBSA[Voridata$INT==7], Voridata$CONC[Voridata$INT==7], pch=21, bg="orchid1")
legend("topleft", legend=c( "no Co-Medication ", "no Interaction", 
                            "CYP2C19-Hemmer (Omeprazol)", 
                            "CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP3A4-Induktoren (Glucocorticoide)",
                            "CYP2C19-Hemmer (Omeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Hemmer (Erythromycin)",
                            "CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Induktoren (Glucocorticoide)",
                            "CYP2C9-Hemmer (Co-Trimazol)"),
       pch=c(21,21,21,21,21,21,21,21), pt.bg=c("lightblue","grey50", "maroon", "blue", "red", "chartreuse1", "orchid1"), cex=0.65)

plot(Voridata$dWT, Voridata$CONC ,xlab="ln(WT) - mean(ln(WT))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Voriconazole measured concentration vs. WT")
polygon(x,y, angle=45, col="grey90")
points(Voridata$dWT[is.na(Voridata$INT)], Voridata$CONC[is.na(Voridata$INT)], pch=21, bg="grey50")
points(Voridata$dWT[is.na(dfCOV$CO.Medi)], Voridata$CONC[is.na(dfCOV$CO.Medi)], pch=21, bg="lightblue")
points(Voridata$dWT[Voridata$INT==1], Voridata$CONC[Voridata$INT==1], pch=21, bg="green")
points(Voridata$dWT[Voridata$INT==2], Voridata$CONC[Voridata$INT==2], pch=21, bg="blue")
points(Voridata$dWT[Voridata$INT==4], Voridata$CONC[Voridata$INT==4], pch=21, bg="maroon")
points(Voridata$dWT[Voridata$INT==5], Voridata$CONC[Voridata$INT==5], pch=21, bg="red")
points(Voridata$dWT[Voridata$INT==7], Voridata$CONC[Voridata$INT==7], pch=21, bg="orchid1")
legend("topleft", legend=c( "no Co-Medication ", "no Interaction", 
                            "CYP2C19-Hemmer (Omeprazol)", 
                            "CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP3A4-Induktoren (Glucocorticoide)",
                            "CYP2C19-Hemmer (Omeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Hemmer (Erythromycin)",
                            "CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Induktoren (Glucocorticoide)",
                            "CYP2C9-Hemmer (Co-Trimazol)"),
       pch=c(21,21,21,21,21,21,21,21), pt.bg=c("lightblue","grey50", "maroon", "blue", "red", "chartreuse1", "orchid1"), cex=0.65)

plot(Voridata$dHEIGHT, Voridata$CONC ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Voriconazole measured concentration vs. HEIGHT")
polygon(x,y, angle=45, col="grey90")
points(Voridata$dHEIGHT[is.na(Voridata$INT)], Voridata$CONC[is.na(Voridata$INT)], pch=21, bg="grey50")
points(Voridata$dHEIGHT[is.na(dfCOV$CO.Medi)], Voridata$CONC[is.na(dfCOV$CO.Medi)], pch=21, bg="lightblue")
points(Voridata$dHEIGHT[Voridata$INT==1], Voridata$CONC[Voridata$INT==1], pch=21, bg="green")
points(Voridata$dHEIGHT[Voridata$INT==2], Voridata$CONC[Voridata$INT==2], pch=21, bg="blue")
points(Voridata$dHEIGHT[Voridata$INT==4], Voridata$CONC[Voridata$INT==4], pch=21, bg="maroon")
points(Voridata$dHEIGHT[Voridata$INT==5], Voridata$CONC[Voridata$INT==5], pch=21, bg="red")
points(Voridata$dHEIGHT[Voridata$INT==7], Voridata$CONC[Voridata$INT==7], pch=21, bg="orchid1")
legend("topleft", legend=c( "no Co-Medication ", "no Interaction", 
                            "CYP2C19-Hemmer (Omeprazol)", 
                            "CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP3A4-Induktoren (Glucocorticoide)",
                            "CYP2C19-Hemmer (Omeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Hemmer (Erythromycin)",
                            "CYP2C19-Hemmer (Omeprazol, Esomeprazol), CYP2C9-Hemmer (Co-Trimazol), CYP3A4-Induktoren (Glucocorticoide)",
                            "CYP2C9-Hemmer (Co-Trimazol)"),
       pch=c(21,21,21,21,21,21,21,21), pt.bg=c("lightblue","grey50", "maroon", "blue", "red", "chartreuse1", "orchid1"), cex=0.65)



#----------------------------------------------------------------------------------------------------------------
#dataframe for patients unter 2 years
df.under2 <- df.Voridata
Voridata <- df.Voridata

#remova all data when age >730days (2years)
remove <- c()
for(i in seq(1,length(df.under2$ID))){
  if(all(df.under2[i,24]>=730)){
    remove <- c(remove,i)
  }
}
df.under2 <- df.under2[-remove,]

mean(df.under2$CONC, na.rm=TRUE)
sqrt(var(df.under2$CONC, na.rm=TRUE))


y2 <- c(5.5,5.5,1,1)   # y-Werte der Eckpunkte
x2 <- c(0,650,650,0)   # x-Werte der Eckpunkte

par(mfrow=c(1,2))
plot(Voridata$DOSE, Voridata$CONC ,xlab="Dose [mg]",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Plasma Conc vs. Dose", ylim=c(0,13), yaxt='n')
polygon(x2,y2, angle=45, col="grey90")
points(Voridata$DOSE[Voridata$ROA==1], Voridata$CONC[Voridata$ROA==1], pch=21, bg="lightblue")
points(Voridata$DOSE[Voridata$ROA==2], Voridata$CONC[Voridata$ROA==2], pch=21, bg="lightblue4")
points(Voridata$DOSE[Voridata$ROA==3], Voridata$CONC[Voridata$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(0,13 ,2), labels=seq(0,13, 2))
axis(side=4, at=seq(-1,15 ,15), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

plot(df.under2$DOSE, df.under2$CONC ,xlab="Dose [mg]",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Plasma Conc vs. Dose < 2 Years", ylim=c(0,4), yaxt='n')
polygon(x2,y2, angle=45, col="grey90")
points(df.under2$DOSE[df.under2$ROA==1], df.under2$CONC[df.under2$ROA==1], pch=21, bg="lightblue")
points(df.under2$DOSE[df.under2$ROA==2], df.under2$CONC[df.under2$ROA==2], pch=21, bg="lightblue4")
points(df.under2$DOSE[df.under2$ROA==3], df.under2$CONC[df.under2$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(-1,5.5 ,1), labels=seq(-1,5.5, 1))
axis(side=3, at=seq(0,250 ,250), labels=NA)
axis(side=4, at=seq(-1, 6, 6), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))


#----------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(df.under2$dAGE, df.under2$CONC ,xlab="ln(AGE) - mean(ln(AGE))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="< 2 Years - AGE", ylim=c(0,4), yaxt='n')
polygon(x,y, angle=45, col="grey90")
points(df.under2$dAGE[df.under2$ROA==1], df.under2$CONC[df.under2$ROA==1], pch=21, bg="lightblue")
points(df.under2$dAGE[df.under2$ROA==2], df.under2$CONC[df.under2$ROA==2], pch=21, bg="lightblue4")
points(df.under2$dAGE[df.under2$ROA==3], df.under2$CONC[df.under2$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(-1,5.5 ,1), labels=seq(-1,5.5, 1))
axis(side=3, at=seq(-5,250 ,250), labels=NA)
axis(side=4, at=seq(-1, 6, 6), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

plot(df.under2$dBSA, df.under2$CONC ,xlab="ln(BSA) - mean(ln(BSA))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="< 2 Years - BSA", ylim=c(0,4),yaxt='n')
polygon(x,y, angle=45, col="grey90")
points(df.under2$dBSA[df.under2$ROA==1], df.under2$CONC[df.under2$ROA==1], pch=21, bg="lightblue")
points(df.under2$dBSA[df.under2$ROA==2], df.under2$CONC[df.under2$ROA==2], pch=21, bg="lightblue4")
points(df.under2$dBSA[df.under2$ROA==3], df.under2$CONC[df.under2$ROA==3], pch=21, bg="maroon")
axis(side=2, at=seq(-1,5.5 ,1), labels=seq(-1,5.5, 1))
axis(side=3, at=seq(-5,250 ,250), labels=NA)
axis(side=4, at=seq(-1, 6, 6), labels=NA)
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

plot(df.under2$dWT, df.under2$CONC ,xlab="ln(WT) - mean(ln(WT))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="< 2 Years - WEIGHT", ylim=c(0,4),yaxt='n')
polygon(x,y, angle=45, col="grey90")
axis(side=2, at=seq(-1,5.5 ,1), labels=seq(-1,5.5, 1))
axis(side=3, at=seq(-5,250 ,250), labels=NA)
axis(side=4, at=seq(-1, 6, 6), labels=NA)
points(df.under2$dWT[df.under2$ROA==1], df.under2$CONC[df.under2$ROA==1], pch=21, bg="lightblue")
points(df.under2$dWT[df.under2$ROA==2], df.under2$CONC[df.under2$ROA==2], pch=21, bg="lightblue4")
points(df.under2$dWT[df.under2$ROA==3], df.under2$CONC[df.under2$ROA==3], pch=21, bg="maroon")
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

plot(df.under2$dHEIGHT, df.under2$CONC ,xlab="ln(HEIGHT) - mean(ln(HEIGHT))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="< 2 Years - HEIGHT", ylim=c(0,4),yaxt='n')
polygon(x,y, angle=45, col="grey90")
axis(side=2, at=seq(-1,5.5 ,1), labels=seq(-1,5.5, 1))
axis(side=3, at=seq(-5,250 ,250), labels=NA)
axis(side=4, at=seq(-1, 6, 6), labels=NA)
points(df.under2$dHEIGHT[df.under2$ROA==1], df.under2$CONC[df.under2$ROA==1], pch=21, bg="lightblue")
points(df.under2$dHEIGHT[df.under2$ROA==2], df.under2$CONC[df.under2$ROA==2], pch=21, bg="lightblue4")
points(df.under2$dHEIGHT[df.under2$ROA==3], df.under2$CONC[df.under2$ROA==3], pch=21, bg="maroon")
legend("topleft", legend=c("ROA = i.v.", "ROA = p.o.", "ROA = i.v. + p.o."),
       pch=c(21,21,21), pt.bg=c("lightblue", "lightblue4", "maroon"))

#----------------------------------------------------------------------------------------------------------------

df.Voridata$TIME <- (df.Voridata$TIME_M-df.Voridata$AGE)*30
df.Voridata <- df.Voridata[order(df.Voridata$ID, df.Voridata$TIME),] # change of time points: new order needed

plot(df.Voridata$TIME, df.Voridata$CONC, xlab="time [months]", ylab="Measured Concentration [mg/L]",type="n")
lines(df.Voridata$TIME, df.Voridata$CONC, type = "l", col=df.Voridata$ID)
points(df.Voridata$TIME, df.Voridata$CONC, pch=21, bg=df.Voridata$ID)

par(mfrow=c(1,4))
plot(df.Voridata$TIME[df.Voridata$ID==1], df.Voridata$CONC[df.Voridata$ID==1], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #1", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==1], df.Voridata$CONC[df.Voridata$ID==1], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==2], df.Voridata$CONC[df.Voridata$ID==2], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #2", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==2], df.Voridata$CONC[df.Voridata$ID==2], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==3], df.Voridata$CONC[df.Voridata$ID==3], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #3", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==3], df.Voridata$CONC[df.Voridata$ID==3], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==5], df.Voridata$CONC[df.Voridata$ID==5], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #5", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==5], df.Voridata$CONC[df.Voridata$ID==5], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==8], df.Voridata$CONC[df.Voridata$ID==8], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #8", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==8], df.Voridata$CONC[df.Voridata$ID==8], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==9], df.Voridata$CONC[df.Voridata$ID==9], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #9", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==9], df.Voridata$CONC[df.Voridata$ID==9], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==10], df.Voridata$CONC[df.Voridata$ID==10], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #10", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==10], df.Voridata$CONC[df.Voridata$ID==10], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==11], df.Voridata$CONC[df.Voridata$ID==11], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #11", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==11], df.Voridata$CONC[df.Voridata$ID==11], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==12], df.Voridata$CONC[df.Voridata$ID==12], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #12", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==12], df.Voridata$CONC[df.Voridata$ID==12], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==14], df.Voridata$CONC[df.Voridata$ID==14], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #14", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==14], df.Voridata$CONC[df.Voridata$ID==14], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==15], df.Voridata$CONC[df.Voridata$ID==15], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #15", ylim=c(0,25))
points(df.Voridata$TIME[df.Voridata$ID==15], df.Voridata$CONC[df.Voridata$ID==15], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==17], df.Voridata$CONC[df.Voridata$ID==17], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #17", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==17], df.Voridata$CONC[df.Voridata$ID==17], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==19], df.Voridata$CONC[df.Voridata$ID==19], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #19", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==19], df.Voridata$CONC[df.Voridata$ID==19], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==20], df.Voridata$CONC[df.Voridata$ID==20], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #20", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==20], df.Voridata$CONC[df.Voridata$ID==20], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==22], df.Voridata$CONC[df.Voridata$ID==22], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #22", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==22], df.Voridata$CONC[df.Voridata$ID==22], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==23], df.Voridata$CONC[df.Voridata$ID==23], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #23", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==23], df.Voridata$CONC[df.Voridata$ID==23], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==24], df.Voridata$CONC[df.Voridata$ID==24], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #24", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==24], df.Voridata$CONC[df.Voridata$ID==24], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==25], df.Voridata$CONC[df.Voridata$ID==25], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #25", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==25], df.Voridata$CONC[df.Voridata$ID==25], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

plot(df.Voridata$TIME[df.Voridata$ID==26], df.Voridata$CONC[df.Voridata$ID==26], xlab="time [months]", 
     ylab="Measured Concentration [mg/L]",type="l", col="lightblue", main="Subject #26", ylim=c(0,10))
points(df.Voridata$TIME[df.Voridata$ID==26], df.Voridata$CONC[df.Voridata$ID==26], pch=21, bg="lightblue")
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

df.Voridata.2Y <- df.Voridata
remove <- c()
for(i in seq(1,length(df.Voridata.2Y$ID))){
  if(all(df.Voridata.2Y[i,24]>=730)){
    remove <- c(remove,i)
  }
}
df.Voridata.2Y <- df.Voridata.2Y[-remove,]

par(mfrow=c(1,1))
plot(df.Voridata.2Y$TIME, df.Voridata.2Y$CONC, xlab="time [months]", ylab="Measured Concentration [mg/L]",type="n")
lines(df.Voridata.2Y$TIME, df.Voridata.2Y$CONC, col=df.Voridata$ID)
points(df.Voridata.2Y$TIME, df.Voridata.2Y$CONC, pch=21, bg=df.Voridata$ID)


#----------------------------------------------------------------------------------------------------------------
#Route of administration vs. AGE
# oral bioavailability -> often malabsorption in children/adolescents

par(mfrow=c(1,2))
plot(Voridata$ROA, Voridata$AGE/365, xlab= "Route of administration" ,ylab="AGE [years]",  type="n",
     main ="ROA vs. AGE", xaxt='n')
points(Voridata$ROA[Voridata$ROA==1], Voridata$AGE[Voridata$ROA==1]/365, pch=21, bg="lightblue")
points(Voridata$ROA[Voridata$ROA==2], Voridata$AGE[Voridata$ROA==2]/365, pch=21, bg="lightblue4")
points(Voridata$ROA[Voridata$ROA==3], Voridata$AGE[Voridata$ROA==3]/365, pch=21, bg="maroon")
axis(side = 1, at=seq(1,3 ,1), labels = c("i.v.", "p.o", "i.v. + p.o."))

plot(df.under2$ROA, df.under2$AGE/365, xlab= "Route of administration" ,ylab="AGE [years]",  type="n",
     main ="< 2 Years - ROA vs. AGE", xaxt='n')
points(df.under2$ROA[df.under2$ROA==1], df.under2$AGE[df.under2$ROA==1]/365, pch=21, bg="lightblue")
points(df.under2$ROA[df.under2$ROA==2], df.under2$AGE[df.under2$ROA==2]/365, pch=21, bg="lightblue4")
points(df.under2$ROA[df.under2$ROA==3], df.under2$AGE[df.under2$ROA==3]/365, pch=21, bg="maroon")
axis(side = 1, at=seq(1,3 ,1), labels = c("i.v.", "p.o", "i.v. + p.o."))

par(mfrow=c(1,2))
plot(Voridata$ROA, Voridata$CONC, xlab= "Route of administration" ,ylab="Concentration [mg/L]",  type="n",
     main ="ROA vs. Conc ", xaxt='n', ylim=c(0,12))
points(Voridata$ROA[Voridata$ROA==1], Voridata$CONC[Voridata$ROA==1], pch=21, bg="lightblue")
points(Voridata$ROA[Voridata$ROA==2], Voridata$CONC[Voridata$ROA==2], pch=21, bg="lightblue4")
points(Voridata$ROA[Voridata$ROA==3], Voridata$CONC[Voridata$ROA==3], pch=21, bg="maroon")
axis(side = 1, at=seq(1,3 ,1), labels = c("i.v.", "p.o", "i.v. + p.o."))

plot(df.under2$ROA, df.under2$CONC, xlab= "Route of administration" ,ylab="Concentration [mg/L]",  type="n",
     main ="< 2 Years - ROA vs. Conc ", xaxt='n', ylim=c(0,4))
points(df.under2$ROA[df.under2$ROA==1], df.under2$CONC[df.under2$ROA==1], pch=21, bg="lightblue")
points(df.under2$ROA[df.under2$ROA==2], df.under2$CONC[df.under2$ROA==2], pch=21, bg="lightblue4")
points(df.under2$ROA[df.under2$ROA==3], df.under2$CONC[df.under2$ROA==3], pch=21, bg="maroon")
axis(side = 1, at=seq(1,3 ,1), labels = c("i.v.", "p.o", "i.v. + p.o."))

#Route of administration
df.Voridata$roa[df.Voridata$ROA==1] <- "i.v."
df.Voridata$roa[df.Voridata$ROA==2] <- "p.o."
df.Voridata$roa[df.Voridata$ROA==3] <- "i.v. + p.o."
df.Voridata$roa[is.na(df.Voridata$ROA)] <- "unknown"

par(mfrow=c(1,2))
boxplot(df.Voridata$CONC~df.Voridata$roa, xlab="",ylab="Measured Concentration [mg/L]",
        main ="Route of Administration ", 
        col=c("lightblue1", "lightblue2", "lightblue3", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

boxplot(df.Voridata$CONC~df.Voridata$roa, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,8),
        col=c("lightblue1", "lightblue2", "lightblue3", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

t.test(df.Voridata$CONC[df.Voridata$roa=="i.v."],df.Voridata$CONC[df.Voridata$roa=="i.v."])


#----------------------------------------------------------------------------------------------------------------
plot(Voridata$AGE/365, Voridata$CONC ,xlab="Age [years]",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Voriconazole measured concentration vs. AGE", ylim=c(0,18)
     )
points(Voridata$AGE[Voridata$ID==1]/365, Voridata$CONC[Voridata$ID==1], pch=21, bg="coral1")
points(Voridata$AGE[Voridata$ID==2]/365, Voridata$CONC[Voridata$ID==2], pch=21, bg="chartreuse1")
points(Voridata$AGE[Voridata$ID==3]/365, Voridata$CONC[Voridata$ID==3], pch=21, bg="cadetblue2")
points(Voridata$AGE[Voridata$ID==4]/365, Voridata$CONC[Voridata$ID==4], pch=21, bg="burlywood2")
#points(Voridata$AGE[Voridata$ID==5]/365, Voridata$CONC[Voridata$ID==5], pch=21, bg="brown2")
#points(Voridata$AGE[Voridata$ID==6]/365, Voridata$CONC[Voridata$ID==6], pch=21, bg="blue2") 
points(Voridata$AGE[Voridata$ID==7]/365, Voridata$CONC[Voridata$ID==7], pch=21, bg="azure4")
#points(Voridata$AGE[Voridata$ID==8]/365, Voridata$CONC[Voridata$ID==8], pch=21, bg="aquamarine2")
points(Voridata$AGE[Voridata$ID==9]/365, Voridata$CONC[Voridata$ID==9], pch=21, bg="dodgerblue2")
points(Voridata$AGE[Voridata$ID==10]/365, Voridata$CONC[Voridata$ID==10], pch=21, bg="deeppink1")
#points(Voridata$AGE[Voridata$ID==11]/365, Voridata$CONC[Voridata$ID==11], pch=21, bg="darkseagreen2")
points(Voridata$AGE[Voridata$ID==12]/365, Voridata$CONC[Voridata$ID==12], pch=21, bg="darkorchid2")
points(Voridata$AGE[Voridata$ID==13]/365, Voridata$CONC[Voridata$ID==13], pch=21, bg="blue2")
points(Voridata$AGE[Voridata$ID==14]/365, Voridata$CONC[Voridata$ID==14], pch=21, bg="darkolivegreen2")
#points(Voridata$AGE[Voridata$ID==15]/365, Voridata$CONC[Voridata$ID==15], pch=21, bg="darkgoldenrod2")
points(Voridata$AGE[Voridata$ID==16]/365, Voridata$CONC[Voridata$ID==16], pch=21, bg="cyan1")
points(Voridata$AGE[Voridata$ID==17]/365, Voridata$CONC[Voridata$ID==17], pch=21, bg="gold1")
#points(Voridata$AGE[Voridata$ID==18]/365, Voridata$CONC[Voridata$ID==18], pch=21, bg="firebrick2")
points(Voridata$AGE[Voridata$ID==19]/365, Voridata$CONC[Voridata$ID==19], pch=21, bg="green4")
#points(Voridata$AGE[Voridata$ID==20]/365, Voridata$CONC[Voridata$ID==20], pch=21, bg="hotpink4")
#points(Voridata$AGE[Voridata$ID==21]/365, Voridata$CONC[Voridata$ID==21], pch=21, bg="mediumblue")
points(Voridata$AGE[Voridata$ID==22]/365, Voridata$CONC[Voridata$ID==22], pch=21, bg="maroon")
points(Voridata$AGE[Voridata$ID==23]/365, Voridata$CONC[Voridata$ID==23], pch=21, bg="orchid1")
points(Voridata$AGE[Voridata$ID==24]/365, Voridata$CONC[Voridata$ID==24], pch=21, bg="navyblue")
points(Voridata$AGE[Voridata$ID==25]/365, Voridata$CONC[Voridata$ID==25], pch=21, bg="plum1")
points(Voridata$AGE[Voridata$ID==26]/365, Voridata$CONC[Voridata$ID==26], pch=21, bg="steelblue1")

plot(Voridata$dAGE, Voridata$CONC ,xlab="ln(AGE) - mean(ln(AGE))",ylab="Measured Concentration [mg/L]",  type="n",
     main ="Voriconazole measured concentration vs. AGE", ylim=c(0,18)
     )
points(Voridata$dAGE[Voridata$ID==1], Voridata$CONC[Voridata$ID==1], pch=21, bg="coral1")
points(Voridata$dAGE[Voridata$ID==2], Voridata$CONC[Voridata$ID==2], pch=21, bg="chartreuse1")
points(Voridata$dAGE[Voridata$ID==3], Voridata$CONC[Voridata$ID==3], pch=21, bg="cadetblue2")
points(Voridata$dAGE[Voridata$ID==4], Voridata$CONC[Voridata$ID==4], pch=21, bg="burlywood2")
#points(Voridata$dAGE[Voridata$ID==5], Voridata$CONC[Voridata$ID==5], pch=21, bg="brown2")
#points(Voridata$dAGE[Voridata$ID==6], Voridata$CONC[Voridata$ID==6], pch=21, bg="blue2") 
points(Voridata$dAGE[Voridata$ID==7], Voridata$CONC[Voridata$ID==7], pch=21, bg="azure4")
#points(Voridata$dAGE[Voridata$ID==8], Voridata$CONC[Voridata$ID==8], pch=21, bg="aquamarine2")
points(Voridata$dAGE[Voridata$ID==9], Voridata$CONC[Voridata$ID==9], pch=21, bg="dodgerblue2")
points(Voridata$dAGE[Voridata$ID==10], Voridata$CONC[Voridata$ID==10], pch=21, bg="deeppink1")
#points(Voridata$dAGE[Voridata$ID==11], Voridata$CONC[Voridata$ID==11], pch=21, bg="darkseagreen2")
points(Voridata$dAGE[Voridata$ID==12], Voridata$CONC[Voridata$ID==12], pch=21, bg="darkorchid2")
points(Voridata$dAGE[Voridata$ID==13], Voridata$CONC[Voridata$ID==13], pch=21, bg="blue2")
points(Voridata$dAGE[Voridata$ID==14], Voridata$CONC[Voridata$ID==14], pch=21, bg="darkolivegreen2")
#points(Voridata$dAGE[Voridata$ID==15], Voridata$CONC[Voridata$ID==15], pch=21, bg="darkgoldenrod2")
points(Voridata$dAGE[Voridata$ID==16], Voridata$CONC[Voridata$ID==16], pch=21, bg="cyan1")
points(Voridata$dAGE[Voridata$ID==17], Voridata$CONC[Voridata$ID==17], pch=21, bg="gold1")
#points(Voridata$dAGE[Voridata$ID==18], Voridata$CONC[Voridata$ID==18], pch=21, bg="firebrick2")
points(Voridata$dAGE[Voridata$ID==19], Voridata$CONC[Voridata$ID==19], pch=21, bg="green4")
#points(Voridata$dAGE[Voridata$ID==20], Voridata$CONC[Voridata$ID==20], pch=21, bg="hotpink4")
#points(Voridata$dAGE[Voridata$ID==21], Voridata$CONC[Voridata$ID==21], pch=21, bg="mediumblue")
points(Voridata$dAGE[Voridata$ID==22], Voridata$CONC[Voridata$ID==22], pch=21, bg="maroon")
points(Voridata$dAGE[Voridata$ID==23], Voridata$CONC[Voridata$ID==23], pch=21, bg="orchid1")
points(Voridata$dAGE[Voridata$ID==24], Voridata$CONC[Voridata$ID==24], pch=21, bg="navyblue")
points(Voridata$dAGE[Voridata$ID==25], Voridata$CONC[Voridata$ID==25], pch=21, bg="plum1")
points(Voridata$dAGE[Voridata$ID==26], Voridata$CONC[Voridata$ID==26], pch=21, bg="steelblue1")

#-------------------------------------------------------------------------------------------------------------------
#Boxplots

#Interaction
# 1 = CYP2C19-Hemmer (Eso)
# 2 = CYP2C19-Hemmer (Eso) + CYP3A4-Induktoren (Gluco)
# 3 = CYP2C19-Hemmer (Ome)
# 4 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Hemmer (Ery)
# 5 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Induktoren (Gluco)
# 6 = CYP2C19-Hemmer (Ome), CYP3A4-Induktoren (Gluco)
# 7 = CYP2C9-Hemmer (Bac)
# 8 = CYP2C9-Hemmer (Ser)
# 9 = CYP3A4-Induktor (Pheno)
# 10 = CYP3A4-Induktoren (Gluco)

df.Voridata$INT[is.na(df.Voridata$INT)] <- 11 #no Co-Medi

par(mfrow=c(1,2))
boxplot(df.Voridata$CONC~df.Voridata$INT, xlab="",ylab="Measured Concentration [mg/L]", 
        main ="Influence of different Co-Medications", 
        col=c("lightcyan", "lightcyan1","lightcyan2", "lightcyan3", "lightcyan4",
              "lightblue", "lightblue1", "lightblue2", "lightblue3", "lightblue4","grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
boxplot(df.Voridata$CONC~df.Voridata$INT, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,7),  
        col=c("lightcyan", "lightcyan1","lightcyan2", "lightcyan3", "lightcyan4",
              "lightblue", "lightblue1", "lightblue2", "lightblue3", "lightblue4","grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

sum(df.Voridata$INT[df.Voridata$INT==1])/1 #n = 5 
sum(df.Voridata$INT[df.Voridata$INT==2])/2 #n = 1
sum(df.Voridata$INT[df.Voridata$INT==3])/3 #n = 16 
sum(df.Voridata$INT[df.Voridata$INT==4])/4 #n = 15
sum(df.Voridata$INT[df.Voridata$INT==5])/5 #n = 1
sum(df.Voridata$INT[df.Voridata$INT==6])/6 #n = 11
sum(df.Voridata$INT[df.Voridata$INT==7])/7 #n = 2
sum(df.Voridata$INT[df.Voridata$INT==8])/8 #n = 1
sum(df.Voridata$INT[df.Voridata$INT==9])/9 #n = 2
sum(df.Voridata$INT[df.Voridata$INT==10])/10 #n = 3
sum(df.Voridata$INT[df.Voridata$INT==11])/11 #n = 194

# n Eso = 6, n Ome = 43, Co-Trim = 18, Ery = 4, Gluco = 16, Ser = 1, Pheno = 2

par(mfrow=c(1,1))
df.Voridata$OME <- NA
df.Voridata$OME[df.Voridata$INT==3] <- "1 = Omeprazole"
df.Voridata$OME[df.Voridata$INT==4] <- "1 = Omeprazole"
df.Voridata$OME[df.Voridata$INT==5] <- "1 = Omeprazole"
df.Voridata$OME[df.Voridata$INT==6] <- "1 = Omeprazole"
df.Voridata$OME[is.na(df.Voridata$OME)] <- "2 = None"

t.test(df.Voridata$CONC~df.Voridata$OME)
#t(53.638) = -1.5373, p-value = 0.1301
boxplot(df.Voridata$CONC~df.Voridata$OME, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,14),
        main ="Omeprazole", col=c("lightblue", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","Omeprazole = 43","(n = 14)", "No Omeprazole = 208","(n = 27)",
                           " ","t (53.628) = -1.5373", "p-value = 0.1301"), cex=0.8)


df.Voridata$ESO <- NA
df.Voridata$ESO[df.Voridata$INT==1] <- "1 = Esomeprazole"
df.Voridata$ESO[df.Voridata$INT==2] <- "1 = Esomeprazole"
df.Voridata$ESO[is.na(df.Voridata$ESO)] <- "2 = None"

t.test(df.Voridata$CONC~df.Voridata$ESO)
#t(4.2339) = 1.9901, p-value = 0.1135
boxplot(df.Voridata$CONC~df.Voridata$ESO, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,14),
        main ="Esomeprazole", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","Esomeprazole = 6","(n = 5)", "No Esomeprazole = 245","(n = 27)",
                           " ","t (4.2339) = 1.9901", "p-value = 0.1135"), cex=0.8)

df.Voridata$BAC <- NA
df.Voridata$BAC[df.Voridata$INT==4] <- "1 = Co-Trimazole"
df.Voridata$BAC[df.Voridata$INT==5] <- "1 = Co-Trimazole"
df.Voridata$BAC[df.Voridata$INT==7] <- "1 = Co-Trimazole"
df.Voridata$BAC[is.na(df.Voridata$BAC)] <- "2 = None"

t.test(df.Voridata$CONC~df.Voridata$BAC)
#t(41.841) = -5.5358, p-value = 1.862e-06
boxplot(df.Voridata$CONC~df.Voridata$BAC, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,14),
        main ="Co-Trimazole", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","Co-Trimazole = 18","(n = 4)", "No Co-Trimazole = 233","(n = 28)"
                           ," ","t (41.841) = -5.5358", "p-value = 1.862e-06"), cex=0.8)

df.Voridata$ERY <- NA
df.Voridata$ERY[df.Voridata$INT==4] <- "1 = Erythromycin"
df.Voridata$ERY[is.na(df.Voridata$ERY)] <- "2 = None"

t.test(df.Voridata$CONC~df.Voridata$ERY)
#t(28.776) = -4.9697, p-value = 2.813e-05
boxplot(df.Voridata$CONC~df.Voridata$ERY, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,14),
        main ="Erythromycin", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","Erythromycin = 4","(n = 2)", "No Erythromycin = 247","(n = 28)"
                           ," ","t (28.776) = -4.9697", "p-value = 2.813e-05"), cex=0.8)

df.Voridata$GLUCO <- NA
df.Voridata$GLUCO[df.Voridata$INT==5] <- "1 = Glucocort."
df.Voridata$GLUCO[df.Voridata$INT==6] <- "1 = Glucocort."
df.Voridata$GLUCO[df.Voridata$INT==10] <- "1 = Glucocort."
df.Voridata$GLUCO[is.na(df.Voridata$GLUCO)] <- "2 = None"

t.test(df.Voridata$CONC~df.Voridata$GLUCO)
#t(19.69) = -1.597, p-value = 0.1262
boxplot(df.Voridata$CONC~df.Voridata$GLUCO, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,14),
        main ="Glucocorticoides", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","Glucocort. = 16","(n = 5)", "No Glucocort. = 235","(n = 28)"
                           ," ", "t (19.69) = -1.597", "p-value = 0.1262"), cex=0.8)

df.Voridata$SER <- NA
df.Voridata$SER[df.Voridata$INT==8] <- "1 = Sertralin"
df.Voridata$SER[is.na(df.Voridata$SER)] <- "2 = None"

boxplot(df.Voridata$CONC~df.Voridata$SER, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,14),
        main ="Sertraline", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","Sertraline = 1","(n = 1)", "No Sertraline = 250","(n = 28)"
                           ," ","no t-test possible"), cex=0.8)

df.Voridata$PHENO <- NA
df.Voridata$PHENO[df.Voridata$INT==9] <- "1 = Phenobarbital"
df.Voridata$PHENO[is.na(df.Voridata$PHENO)] <- "2 = None"

boxplot(df.Voridata$CONC~df.Voridata$PHENO, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,14),
        main ="Phenobarbital", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","Phenobarbital = 2","(n = 1)", "No Phenobarbital = 249","(n = 28)"
                           ," ","no t-test possible"), cex=0.8)

# 1 = CYP2C19-Hemmer (Eso)
# 2 = CYP2C19-Hemmer (Eso) + CYP3A4-Induktoren (Gluco)
# 3 = CYP2C19-Hemmer (Ome)
# 4 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Hemmer (Ery)
# 5 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Induktoren (Gluco)
# 6 = CYP2C19-Hemmer (Ome), CYP3A4-Induktoren (Gluco)
# 7 = CYP2C9-Hemmer (Bac)
# 8 = CYP2C9-Hemmer (Ser)
# 9 = CYP3A4-Induktor (Pheno)
# 10 = CYP3A4-Induktoren (Gluco)

df.Voridata$CYP2C19 <- NA
df.Voridata$CYP2C19[df.Voridata$INT==1] <- "Inhibitor"
df.Voridata$CYP2C19[df.Voridata$INT==2] <- "Inhibitor"
df.Voridata$CYP2C19[df.Voridata$INT==3] <- "Inhibitor"
df.Voridata$CYP2C19[df.Voridata$INT==4] <- "Inhibitor"
df.Voridata$CYP2C19[df.Voridata$INT==5] <- "Inhibitor"
df.Voridata$CYP2C19[df.Voridata$INT==6] <- "Inhibitor"
df.Voridata$CYP2C19[is.na(df.Voridata$CYP2C19)] <- "None"

df.Voridata$CYP2C9 <- NA
df.Voridata$CYP2C9[df.Voridata$INT==4] <- "Inhibitor"
df.Voridata$CYP2C9[df.Voridata$INT==5] <- "Inhibitor"
df.Voridata$CYP2C9[df.Voridata$INT==7] <- "Inhibitor"
df.Voridata$CYP2C9[df.Voridata$INT==8] <- "Inhibitor"
df.Voridata$CYP2C9[is.na(df.Voridata$CYP2C9)] <- "None"

df.Voridata$CYP3A4 <- NA
df.Voridata$CYP3A4[df.Voridata$INT==4] <- "Inhibitor"
df.Voridata$CYP3A4[df.Voridata$INT==2] <- "Inductor"
df.Voridata$CYP3A4[df.Voridata$INT==5] <- "Inductor"
df.Voridata$CYP3A4[df.Voridata$INT==6] <- "Inductor"
df.Voridata$CYP3A4[df.Voridata$INT==9] <- "Inductor"
df.Voridata$CYP3A4[df.Voridata$INT==10] <- "Inductor"
df.Voridata$CYP3A4[is.na(df.Voridata$CYP3A4)] <- "None"

par(mfrow=c(1,3))
t.test(df.Voridata$CONC~df.Voridata$CYP2C19)
#t (62.216) = -0.88803, p-value = 0.3779
boxplot(df.Voridata$CONC~df.Voridata$CYP2C19, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,10),
        main ="CYP2C19 Inhibitors", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","CYP2C19-Inhibitor = 49","(n = 19)", "None = 202","(n = 28)"
                           ," ", "t (62.216) = -0.88803", "p-value = 0.3779"), cex=0.8)

t.test(df.Voridata$CONC~df.Voridata$CYP2C9)
#t (49.977) = -5.8005, p-value = 4.467e-07
boxplot(df.Voridata$CONC~df.Voridata$CYP2C9, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,10),
        main ="CYP2C9 Inhibitors", col=c("lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","CYP2C9-Inhibitor = 19","(n = 4)", "None = 232","(n = 28)"
                           ," ", "t (49.977) = -5.8005", "p-value = 4.467e-07"), cex=0.8)

#Inhibitor: t (31.306) = -5.0248, p-value = 1.952e-05
#Inductor: t (22.853) = -1.4389, p-value = 0.1637

t.test(df.Voridata$CONC[df.Voridata$CYP3A4=="Inhibitor"],df.Voridata$CONC[df.Voridata$CYP3A4=="None"])
t.test(df.Voridata$CONC[df.Voridata$CYP3A4=="Inductor"],df.Voridata$CONC[df.Voridata$CYP3A4=="None"])
boxplot(df.Voridata$CONC~df.Voridata$CYP3A4, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,10),
        main ="CYP3A4 Inhibitors/Inductors", col=c("lightblue4", "lightblue", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("Number of Measurements:","CYP3A4-Inhibitor = 15","(n = 2)", "CYP3A4-Inductor = 18","(n = 7)", "None = 218","(n = 27)"
                           ," ", "Inhibitor:","t (31.306) = -5.0248", "p-value = 1.952e-05",
                                 "Inductor:","t (22.853) = -1.4389", "p-value = 0.1637"), cex=0.8)


#---------------------------------------------------------------------------------------------------------------------
#Diagnose
#1 = pulmonale Aspergillose, 2 = invasive Aspergillose, 3 = Pneumocystes jiroveci Pneumonie, 4 = Aspergillom
#5 = Candidämie, 6 = Candida Meningitis, 7 = Varizella-Zoster-Virus (VZV) Pneumonie, 8 = Infektionsprophylaxe, 9 = other

df.Voridata$DIAG.nb[is.na(df.Voridata$DIAG.nb)] <- 10 

par(mfrow=c(1,2))
boxplot(df.Voridata$CONC~df.Voridata$DIAG.nb, xlab="",ylab="Measured Concentration [mg/L]",
        main ="Influence of differnt Diagnoses", 
        col=c("lightcyan", "lightcyan1","lightcyan2", "lightcyan3", "lightcyan4",
              "lightblue", "lightblue1", "lightblue2", "lightblue3", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("1   = pulmonary aspergillosis","2   = invasive aspergillosis","3   = Pneumocystes jiroveci pneumonia","4   = Aspergillom","5   = Candidaemia",
                           "6   = Candida meningitis","7   = VZV pneumonia", "8   = Infection prophylaxis" , "9   = Other", "10 = Diagnosis not known"), cex=0.8)

boxplot(df.Voridata$CONC~df.Voridata$DIAG.nb, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,10),
        col=c("lightcyan", "lightcyan1","lightcyan2", "lightcyan3", "lightcyan4",
              "lightblue", "lightblue1", "lightblue2", "lightblue3", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

#---------------------------------------------------------------------------------------------------------------------
#Comorbidity
#1 = ALL, 2 = AML, 3 = Betatalassämie major, 4 = CF, 5 = chronische/septische Granulomatose, 
#6 = aplastische Anämie, 7 = Hepatitis, 8 = other Cancer , 9 = other

df.Voridata$CO.DIAG.nb[is.na(df.Voridata$CO.DIAG.nb)] <- 10 

par(mfrow=c(1,2))
boxplot(df.Voridata$CONC~df.Voridata$CO.DIAG.nb, xlab="",ylab="Measured Concentration [mg/L]",
        main ="Influence of different Co-Morbidities", 
        col=c("lightcyan", "lightcyan1","lightcyan2", "lightcyan3", "lightcyan4",
              "lightblue", "lightblue1", "lightblue2", "lightblue3", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("1   = ALL","2   = AML","3   = Betatalassämie major","4   = CF","5   = chronische/septische Granulomatose",
       "6   = aplastische Anämie","7   = Hepatitis", "8   = Other Cancer" , "9   = Other", "10 = No Co-Morbidity"), cex=0.8)

boxplot(df.Voridata$CONC~df.Voridata$CO.DIAG.nb, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,8),
        col=c("lightcyan", "lightcyan1","lightcyan2", "lightcyan3", "lightcyan4",
              "lightblue", "lightblue1", "lightblue2", "lightblue3", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

#---------------------------------------------------------------------------------------------------------------------
#Transplantation
#1 = SZT, 2 = KMT, 
#Graft reaction
#1 = GvHD

df.Voridata$TRANSP <- NA
df.Voridata$TRANSP[df.Voridata$SZT==1] <- 1
df.Voridata$TRANSP[df.Voridata$GvHD==1] <- 2
df.Voridata$TRANSP[df.Voridata$KMT==1] <- 3
df.Voridata$TRANSP[df.Voridata$GvHD==2] <- 4
df.Voridata$TRANSP[is.na(df.Voridata$TRANSP)] <- 5

t.test(df.Voridata$CONC[df.Voridata$TRANSP==1],df.Voridata$CONC[df.Voridata$TRANSP==5])
t.test(df.Voridata$CONC[df.Voridata$TRANSP==2],df.Voridata$CONC[df.Voridata$TRANSP==5])
t.test(df.Voridata$CONC[df.Voridata$TRANSP==3],df.Voridata$CONC[df.Voridata$TRANSP==5])
t.test(df.Voridata$CONC[df.Voridata$TRANSP==4],df.Voridata$CONC[df.Voridata$TRANSP==5])

par(mfrow=c(1,2))
boxplot(df.Voridata$CONC~df.Voridata$TRANSP, xlab="",ylab="Measured Concentration [mg/L]",
        main ="Influence of Transplantation", 
        col=c("lightblue1", "lightblue2","lightblue3","lightblue4", "grey80"))
       
abline(h=c(1,5.5) ,col="grey50", lty="twodash")
legend("topleft", legend=c("1 = SCT","2 = SCT with GvHD","3 = BMT","4 = BMT with Engraftssyndrome","5 = no Transplantation"), cex=0.8)

boxplot(df.Voridata$CONC~df.Voridata$TRANSP, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,8), 
        col=c("lightblue1", "lightblue2","lightblue3","lightblue4", "grey80"))
abline(h=c(1,5.5) ,col="grey50", lty="twodash")


#---------------------------------------------------------------------------------------------------------------------
#IPS
#1 = ambulant, 2 = IPS
df.Voridata$IPS[df.Voridata$IPS==1] <- "ambulant"
df.Voridata$IPS[df.Voridata$IPS==2] <- "ICU"
df.Voridata$IPS[is.na(df.Voridata$IPS)] <- "unknown"

par(mfrow=c(1,2))
boxplot(df.Voridata$CONC~df.Voridata$IPS, xlab="",ylab="Measured Concentration [mg/L]",
        main ="Influence of differnet Forms of Patient Care", 
        col=c("lightblue", "lightblue4", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

boxplot(df.Voridata$CONC~df.Voridata$IPS, xlab="",ylab="Measured Concentration [mg/L]", ylim=c(0,8),
        col=c("lightblue", "lightblue4", "grey80")
)
abline(h=c(1,5.5) ,col="grey50", lty="twodash")

#---------------------------------------------------------------------------------------------------------------------
t.test(df.Voridata$CONC[df.Voridata$IPS=="ambulant"],df.Voridata$CONC[df.Voridata$IPS=="ICU"])


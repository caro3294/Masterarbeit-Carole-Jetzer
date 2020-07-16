# Voriconazole Dataframe
# autor: Carole Jetzer
# next run Voriconazole_2_Covariates.R

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

### Give some info about the data
NumberRows <- 251 # rows in csv, header not included

### Read data from csv
FileName <- "Voriconazol_Data.csv"
df <- read.csv(file.path("D:/BioMed IT CSV/CSV_Voriconazol",FileName,fsep = .Platform$file.sep),
               header=TRUE, sep=";",nrows=NumberRows,na.strings="",stringsAsFactors = FALSE)
Headers <- colnames(df)

### Correct DATES
#format="%d.%m.%Y --> the format the date is entered in the data
df$GEBURTSTAG <- as.Date(df$Geburtsdatum, format="%d.%m.%Y")
df$ENTNAHMETAG <- as.Date(df$Entnahmedatum, format="%d.%m.%Y")
df$APPLIKATIONSTAG <- as.Date(df$Datum.Applikation, format="%d.%m.%Y")

### Correct TIMES: 
#Age
df$AgeMessung <- df$ENTNAHMETAG - df$GEBURTSTAG #days
df$AgeDose <- df$APPLIKATIONSTAG - df$GEBURTSTAG #days

#Zeitpunkt der Messung
df$Entnahmezeit[is.na(df$Entnahmezeit)] <- "00:00:00"
df$ENTNAHME.Time <- df$ENTNAHMETAG - df$GEBURTSTAG  +
  as.numeric(substr(df$Entnahmezeit,1,2))/24 + 
  as.numeric(substr(df$Entnahmezeit,4,5))/60/24 

#Zeitpunkt der Applikation
df$Uhrzeit.Applikation[is.na(df$Uhrzeit.Applikation)] <- "00:00:00"
df$APPLICATION.Time <- df$APPLIKATIONSTAG - df$GEBURTSTAG +
  as.numeric(substr(df$Uhrzeit.Applikation,1,2))/24 + 
  as.numeric(substr(df$Uhrzeit.Applikation,4,5))/60/24 

#Zeitpunkt der Applikation zusätzliche Dosis
df$Uhrzeit.Applikation.zus.tzliche.Dosis.po[is.na(df$Uhrzeit.Applikation.zus.tzliche.Dosis.po)] <- "00:00:00"
df$APPLICATION.Time2 <- df$APPLIKATIONSTAG - df$GEBURTSTAG +
  as.numeric(substr(df$Uhrzeit.Applikation.zus.tzliche.Dosis.po,1,2))/24 + 
  as.numeric(substr(df$Uhrzeit.Applikation.zus.tzliche.Dosis.po,4,5))/60/24 

df$APPLICATION.Time2[is.na(df$zus.tzlich.Dosis.po..mg.)] <- NA
  
#Replace the measured concentration of <0.1 with NA
df <- data.frame(lapply(df, function(x) {
  gsub("<0.1", NA, x)
}))

df2 <- data.frame(ID=as.numeric(df$Geburtsdatum),
                  INF_TIME=as.numeric(as.character(df$Infusionsdauer..h.)),      #infusion duration,
                  TIME_A=as.numeric(as.character(df$APPLICATION.Time)),          # Time of the application
                  TIME_Apo=as.numeric(as.character(df$APPLICATION.Time2)),       # Time of the application of additional po dose
                  TIME_M=as.numeric(as.character(df$ENTNAHME.Time)),             # Time of the measurement
                  TIMEdiff=(as.numeric(as.character(df$ENTNAHME.Time))) 
                  - (as.numeric(as.character(df$APPLICATION.Time))),
                  INTERVAL=df$Dosierungsintervall..h.,                           #TAU
                  DOSE=as.numeric(as.character(df$letzte.Dosis..mg.)),           #Dose
                  DOSEpo=as.numeric(as.character(df$zus.tzlich.Dosis.po..mg.)),  #additional Dose p.o.
                  CONC=as.numeric(as.character(df$ResultatNumerisch)),           #[mg/L]
                  ROA=df$ROA,                                                    #route of administration
                  AGE=(as.numeric(as.character(df$AgeMessung))),
                  AGE_M=(as.numeric(as.character(df$AgeMessung))),
                  AGE_A=(as.numeric(as.character(df$AgeDose))),
                  HEIGHT=(as.numeric(as.character(df$Gr.sse..cm.))),             #cm
                  WT=(as.numeric(as.character(df$K.rpergewicht..kg.))),          #kg
                  DIAGNOSE=df$Diagnose,
                  CO.Medi=df$Co.Medikation,
                  Interaction=df$CYP.Induktor...Inhibitor,
                  GEBURTSDATUM=df$Geburtsdatum
                  )

df.Vori <- df2
df.Vori <- setorder(df.Vori, ID, na.last=FALSE)

mean(df2$TIMEdiff)*24
sqrt(var((df2$TIMEdiff)))*24

df.GEB <- data.frame(ID=df.Vori$ID,
                     GEB=df.Vori$GEBURTSDATUM)

write.csv(df.GEB,"D:/BioMed IT CSV/CSV_Voriconazol\\Voriconazol_GEB_ID.csv", row.names = FALSE)

#------------------------------------------------------------------------------------------

#remove row with no application date
remove <- c()
for(i in seq(1,length(df2$ID))){
  if(all(is.na(df2[i,3]))){
    remove <- c(remove,i)
  }
}
df2 <- df2[-remove,]

#remove row when application time is after measurement time 
remove2 <- c()
for(i in seq(1,length(df2$ID))){
  if(all((df2[i,6]<=0))){
    remove2 <- c(remove2,i)
  }
}
df2 <- df2[-remove2,]

#remove row when ROA is NA  
remove3 <- c()
for(i in seq(1,length(df2$ID))){
  if(all(is.na(df2[i,11]))){
    remove3 <- c(remove3,i)
  }
}
df2 <- df2[-remove3,]

#remove row when ROA = i.v. and inf_Time is NA  
remove4 <- c()
for(i in seq(1,length(df2$ID))){
  if(all((df2[i,11]=="iv"))){
    if(all(is.na(df2[i,2]))){
    remove4 <- c(remove4,i)
    }
  }
}
df2 <- df2[-remove4,]

#remove row with no Dose
remove5 <- c()
for(i in seq(1,length(df2$ID))){
  if(all(is.na(df2[i,8]))){
    remove5 <- c(remove5,i)
  }
}
df2 <- df2[-remove5,]

#change the ROA to numeric
df2$ROA <- as.numeric(df2$ROA) # iv = 2 , po = 4 , iv+po =  3
df2$ROA[df2$ROA==2] <- 1 # iv = 1
df2$ROA[df2$ROA==4] <- 2 # po = 2

df2$Interaction.nb <- as.numeric(df2$Interaction) 
# 1 = CYP2C19-Hemmer (Ome)
# 2 = CYP2C19-Hemmer (Ome) + CYP3A4-Induktoren (Gluco)
# 3 = leer
# 4 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Hemmer (Ery)
# 5 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Induktoren (Gluco)
# 6 = leer
# 7 = CYP2C9-Hemmer (Bac)
# 8 = leer

df2$Interaction.nb[df2$Interaction.nb==6] <- 2
df2$Interaction.nb[df2$Interaction.nb==8] <- 7
# 1 = CYP2C19-Hemmer (Ome + Eso)
# 2 = CYP2C19-Hemmer (Ome + Eso) + CYP3A4-Induktoren (Gluco)
# 4 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Hemmer (Ery)
# 5 = CYP2C19-Hemmer (Ome), CYP2C9-Hemmer (Bac), CYP3A4-Induktoren (Gluco)
# 7 = CYP2C9-Hemmer (Bac)

df2 <- setorder(df2, ID, na.last=FALSE)

VoriDataCorr <- data.frame(ID=df2$ID,
                           INF_TIME=df2$INF_TIME,
                           TIME_A=df2$TIME_A,
                           TIME_po=df2$TIME_Apo,
                           TIME_M=df2$TIME_M,
                           INTERVAL=df2$INTERVAL,
                           DOSE=df2$DOSE,
                           DOSEpo=df2$DOSEpo, 
                           CONC=df2$CONC,
                           ROA=df2$ROA,
                           INT=df2$Interaction.nb
)

write.csv(VoriDataCorr,"D:/BioMed IT CSV/CSV_Voriconazol\\Voriconazol_CORRdata.csv", row.names = FALSE)

dfCOV <- data.frame(ID=df2$ID,
                    AGE=df2$AGE,
                    BSA=sqrt((df2$HEIGHT*df2$WT)/3600), # Mosteller formula:  BSA = sqrt (([cm] x [kg]) /3600)
                    HEIGHT=df2$HEIGHT,
                    WT=df2$WT,
                    INT=df2$Interaction.nb
)

dfCOV$lAGE <- log(dfCOV$AGE)
dfCOV$lBSA <- log(dfCOV$BSA)
dfCOV$lHEIGHT <- log(dfCOV$HEIGHT)
dfCOV$lWT <- log(dfCOV$WT)

dfCOV$dAGE <- dfCOV$lAGE - mean(dfCOV$lAGE,na.rm=T)
dfCOV$dBSA <- dfCOV$lBSA - mean(dfCOV$lBSA,na.rm=T)
dfCOV$dHEIGHT <- dfCOV$lHEIGHT - mean(dfCOV$lHEIGHT,na.rm=T)
dfCOV$dWT <- dfCOV$lWT- mean(dfCOV$lWT, na.rm=T)

VoriCOV <- data.frame(ID=df2$ID,
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
                      dWT=dfCOV$dWT,
                      DIAGNOSE=df2$DIAGNOSE,
                      CO.Medi=df2$CO.Medi,
                      INT=df2$Interaction.nb
                      
)

write.csv(VoriCOV,"D:/BioMed IT CSV/CSV_Voriconazol\\Voriconazole_Covariates.csv", row.names = FALSE)

#--------------------------------------------------------------------------------------------------------------------------
#change the ROA to numeric
df.Vori$ROA <- as.numeric(df.Vori$ROA) # iv = 2 , po = 4 , iv+po =  3
df.Vori$ROA[df.Vori$ROA==2] <- 1 # iv = 1
df.Vori$ROA[df.Vori$ROA==4] <- 2 # po = 2

#--------------------------------------------------------------------------------------------------------------------------
#change the Interaction to numeric
df.Vori$Interaction.nb <- as.numeric(df.Vori$Interaction) 
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

#--------------------------------------------------------------------------------------------------------------------------
#change IPS to numeric
df.IPS <- data.frame(ID=df.Vori$ID, 
                     IPS <- as.numeric(df$Ambulat.IPS))
# 1 = ambulant
# 2 = IPS
#--------------------------------------------------------------------------------------------------------------------------
df.Vori$DIAGNOSE.nb <- as.numeric(df.Vori$DIAGNOSE) 
# 1 = leer
# 2 = ALL, ARDS, pulmonale Aspergillose
# 3 = ALL, pulmonale Aspergillose
# 4 = ALL, invasive Aspergillose
# 5 = ALL, Mittelschwere Pneumocystes jiroveci Pneumonie
# 6 = ALL, Stammzellentransplantation (SZT), GvHD (Graft-versus-Host-Disease), pulmonale Aspergillose
# 7 = ALL, Stammzellentransplantation (SZT), GvHD (Graft-versus-Host-Disease), pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==7] <- 6
# 8 = leer
# 9 =  AML, invasive Aspergillose
# 10 = AML, multile Lungengrundherde mit Pilzbefall, septische Embolien 
# 11 = AML, Stammzellentransplantation (SZT), semi-invasive Aspergillose
# 12 = AML, Astrozytom, Aspergillom
# 13 = AML, Astrozytom, Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==13] <- 12
# 14 = AML, Astrozytom, invasive pulmonale Aspergillose
# 15 = Beta Thalassämia major, SZT, GvHD, Candidämie
# 16 = Candida Meningitis
# 17 = CF, pulmonale Aspergillose
# 18 = CF, Aspergillus-Bronchitis
# 19 = CF, Aspergillus Bronchitis, Nachweis von Aspergillus spp in BAL
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==19] <- 18
# 20 = chronische Granulomatose, SZT, pulmonale Aspergillose
# 21 = Chronische Granulomatose, SZT, GvHD
# 22 = Dermatomyositis (Muskelerkrankung)
# 23 = Diabetes, AML, Knochenmarktransplantation (KMT), Engraftssyndrom, Pneumonie (Fieber ohne Fokus)
# 24 = Hepatitis, schwere Neutropenie
# 25 = Hyper-IgM-Syndrom, Varizella-Zoster-Virus (VZV) Pneumonie
# 26 = Medulloblastom, invasive Aspergillose
# 27 = Multiorganversagen, immunologisches Multisystemgeschehen, invasive Aspergillose
# 28 = schwere aplastische Anäme, mögliche pulmonale Pilzinfektion
# 29 = schwere aplastische Anäme, SZT
# 30 = schwere aplastische Anäme, SZT, Fieber ohne Fokus, Infektionsprophylaxe
# 31 = septische Granulomatose Leberabszess, bakteriell, mykotisch
# 32 = septische Granulomatose, SZT, GvHD, pulomanel Aspergillose
# 33 = Knochenmarktransplantation (KMT),
# 34 = rezidive AML, SZT Astrozytom, Pulmonale Aspergillose im Transplantationsverlauf
# 35 = T Lymphoblastisches Lymphom, pulmonale Aspergillose
# 36 = Tumor Erkrankung, Port-Infekt
# 37 = X-linked Lympho-Proliferatives Syndrom, SZT, GvHD


#1 = pulmonale Aspergillose, 2 = invasive Aspergillose, 3 = Pneumocystes jiroveci Pneumonie, 4 = Aspergillom
#5 = Candidämie, 6 = Candida Meningitis, 7 = Varizella-Zoster-Virus (VZV) Pneumonie, 8 = Infektionsprophylaxe, 9 = other
df.Vori$DIAGNOSE.NUM <- df.Vori$DIAGNOSE.nb

df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==2]  <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==3]  <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==6]  <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==17] <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==18] <- 1 #pulmonale Aspergillose 
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==19] <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==20] <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==28] <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==32] <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==34] <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==35] <- 1 #pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==4]  <- 2 #invasive Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==9]  <- 2 #invasive Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==11] <- 2 #invasive Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==14] <- 2 #invasive pulmonale Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==26] <- 2 #invasive Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==27] <- 2 #invasive Aspergillose
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==5]  <- 3 #Pneumocystes jiroveci Pneumonie
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==12] <- 4 #Aspergillom
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==13] <- 4 #Aspergillom
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==15] <- 5 #Candidämie
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==16] <- 6 #Candida Meningitis
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==25] <- 7 #Varizella-Zoster-Virus (VZV) Pneumonie
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==30] <- 8 #Infektionsprophylaxe
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==23] <- 9 #Other
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==31] <- 9 #Other 
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==36] <- 9 #Other
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==21] <- NA #only co-Diagnosis
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==22] <- NA #only co-Diagnosis
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==24] <- NA #only co-Diagnosis
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==29] <- NA #only co-Diagnosis
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==33] <- NA #only co-Diagnosis
df.Vori$DIAGNOSE.nb[df.Vori$DIAGNOSE.nb==37] <- NA #only co-Diagnosis

#1 = ALL, 2 = AML, 3 = Betatalassämie major, 4 = CF, 5 = chronische/septische Granulomatose, 
#6 = aplastische Anämie, 7 = Hepatitis, 8 = other Cancer , 9 = other
df.Vori$CO.DIAGNOSE.nb <- as.numeric(df.Vori$DIAGNOSE) 

df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==2]  <- 1 #ALL
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==3]  <- 1 #ALL
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==4]  <- 1 #ALL
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==5]  <- 1 #ALL
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==6]  <- 1 #ALL 
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==7]  <- 1 #ALL
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==9]  <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==10] <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==11] <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==12] <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==13] <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==14] <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==23] <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==34] <- 2 #AML
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==15] <- 3 #Betatalassämie major
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==17] <- 4 #CF
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==18] <- 4 #CF
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==19] <- 4 #CF
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==20] <- 5 #chronische Granulomatose
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==21] <- 5 #chronische Granulomatose
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==31] <- 5 #chronische Granulomatose
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==32] <- 5 #chronische Granulomatose
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==28] <- 6 #aplastische Anämie
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==29] <- 6 #aplastische Anämie
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==30] <- 6 #aplastische Anämie
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==24] <- 7 #Hepatitis
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==26] <- 8 #Other Cancer 
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==35] <- 8 #Other Cancer 
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==36] <- 8 #Other Cancer 
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==25] <- 9 #Other
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==27] <- 9 #Other
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==37] <- 9 #Other
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==16] <- NA #No Comorbidity
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==22] <- NA #No Comorbidity
df.Vori$CO.DIAGNOSE.nb[df.Vori$CO.DIAGNOSE.nb==33] <- NA #No Comorbidity

dfDIAG <- data.frame(ID=df.Vori$ID,
                     CO.DIAGNOSE.nb=df.Vori$CO.DIAGNOSE.nb,
                     DIAGNOSE.nb=df.Vori$DIAGNOSE.nb
)

dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==1] <- "ALL"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==2] <- "AML"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==3] <- "Betathal"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==4] <- "CF"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==5] <- "chron.Granu"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==6] <- "apl.Anäm"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==7] <- "Hep"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==8] <- "other.Cancer"
dfDIAG$CO.DIAGNOSE.nb[dfDIAG$CO.DIAGNOSE.nb==9] <- "other"

dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==1] <- "pulm.Asp"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==2] <- "inv.Asp"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==3] <- "Pneumocystes"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==4] <- "Aspergillom"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==5] <- "Candidämie"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==6] <- "Cand.Mening"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==7] <- "VZV.Pneum"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==8] <- "Prophyl"
dfDIAG$DIAGNOSE.nb[dfDIAG$DIAGNOSE.nb==9] <- "other"

#--------------------------------------------------------------------------------------------------------------------------
df.Transplant <- data.frame(ID  =df.Vori$ID,
                            SZT =df.Vori$DIAGNOSE.NUM,
                            KMT =df.Vori$DIAGNOSE.NUM,
                            GvHD=df.Vori$DIAGNOSE.NUM)

#SZT: 6 11 15 20 21 29 30 32 34 37 
#KMT: 23 33 
#GvHD: 6 15 21 23 32 37
  
df.Transplant$SZT[df.Transplant$SZT==6]   <- 1  #SZT
df.Transplant$SZT[df.Transplant$SZT==11]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==15]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==20]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==21]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==29]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==30]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==32]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==34]  <- 1 
df.Transplant$SZT[df.Transplant$SZT==37]  <- 1
df.Transplant$SZT[df.Transplant$SZT!=1]   <- NA

df.Transplant$KMT[df.Transplant$KMT==23]  <- 1 #KMT
df.Transplant$KMT[df.Transplant$KMT==33]  <- 1 
df.Transplant$KMT[df.Transplant$KMT!=1]   <- NA

df.Transplant$GvHD[df.Transplant$GvHD==6]   <- 1 #GvHD
df.Transplant$GvHD[df.Transplant$GvHD==15]  <- 1 
df.Transplant$GvHD[df.Transplant$GvHD==21]  <- 1 
df.Transplant$GvHD[df.Transplant$GvHD==23]  <- 2 
df.Transplant$GvHD[df.Transplant$GvHD==32]  <- 1 
df.Transplant$GvHD[df.Transplant$GvHD==37]  <- 1 
df.Transplant$GvHD[df.Transplant$GvHD>2]   <- NA 

#--------------------------------------------------------------------------------------------------------------------------
df.AGE <- data.frame(ID=df.Vori$ID,
                     AGE=df.Vori$AGE)
#----------------------------------------------------------------------------------------------------------------------
VoriData <- data.frame(ID=df.Vori$ID,
                       GEB=df$GEBURTSTAG,
                       INF_TIME=df.Vori$INF_TIME,
                       TIME_A=df.Vori$TIME_A,
                       DAY_A=df$APPLIKATIONSTAG,
                       TIME_po=df.Vori$TIME_Apo,
                       TIME_M=df.Vori$TIME_M,
                       DAY_M=df$ENTNAHMETAG,
                       TIME_DIFF=df.Vori$TIMEdiff,
                       INTERVAL=df.Vori$INTERVAL,
                       DOSE=df.Vori$DOSE,
                       DOSEpo=df.Vori$DOSEpo, 
                       CONC=df.Vori$CONC,
                       ROA=df.Vori$ROA,
                       INT=df.Vori$Interaction.nb,
                       DIAG.nb=df.Vori$DIAGNOSE.nb,
                       DIAG=dfDIAG$DIAGNOSE.nb,
                       CO.DIAG.nb=df.Vori$CO.DIAGNOSE.nb,
                       CO.DIAG=dfDIAG$CO.DIAGNOSE.nb,
                       SZT=df.Transplant$SZT,
                       KMT=df.Transplant$KMT,
                       GvHD=df.Transplant$GvHD,
                       IPS=df.IPS$IPS....as.numeric.df.Ambulat.IPS.,
                       AGE=df.Vori$AGE,
                       BSA=sqrt((df.Vori$HEIGHT*df.Vori$WT)/3600), # Mosteller formula:  BSA = sqrt (([cm] x [kg]) /3600)
                       HEIGHT=df.Vori$HEIGHT,
                       WT=df.Vori$AGE
                       )

VoriData$lAGE <- log(VoriData$AGE)
VoriData$lBSA <- log(VoriData$BSA)
VoriData$lHEIGHT <- log(VoriData$HEIGHT)
VoriData$lWT <- log(VoriData$WT)

VoriData$dAGE <- VoriData$lAGE - mean(VoriData$lAGE,na.rm=T)
VoriData$dBSA <- VoriData$lBSA - mean(VoriData$lBSA,na.rm=T)
VoriData$dHEIGHT <- VoriData$lHEIGHT - mean(VoriData$lHEIGHT,na.rm=T)
VoriData$dWT <- VoriData$lWT- mean(VoriData$lWT, na.rm=T)

write.csv(VoriData,"D:/BioMed IT CSV/CSV_Voriconazol\\Voriconazol_data2.csv", row.names = FALSE)


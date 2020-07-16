# Gentamicin nlmixr
# autor: Carole Jetzer
# next run Gentamixin_2_nlmixr_k10_V.R 

#Informations about Gentamicin
"Administered intravenously as short infusion
Half-life in children: 3-6 h
Dosing interval: 24 h
Trough levels need be < 2 mg/l to reduce ear and kidney tox
Peak levels < 12 mg/l may further reduce tox
Elimination mainly / only by renal excretion of unchanged drug"
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)

## ----libraries, cache=F, message=F, warning=F, include=F-----------------

#Installiere alle nötigen Packages falls noch nicht installiert und lade sie dann
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

rm(list=ls()) 
par.default <- par(no.readonly = TRUE) #Speicher die standart Einstellungen f?r die Plots
setwd("C:/Users/carol/OneDrive/Studium/ETH/Master/Masterarbeit") #Setze den Working-directory

### Give some info about the data
NumberRows <- 73 # rows in csv, header not included
BelowLimit <- "<0.3" #data is under detection limit
BelowLimitReplacement <- 0.3 
BelowLimitMethod <- "SSR0" #  sum of squared residuals (SSR) : SSR0 for keeping below limit concentrations and setting residue to 0 if predicted is also bleow limit
# "NA" for setting below limit concentrations to NA. 

### Read data from csv
FileName <- "Gentamicin_Data_For_Publication.csv"
df <- read.csv(file.path("D:/BioMed IT CSV/CSV_Gentamicin",FileName,fsep = .Platform$file.sep),
               header=TRUE, sep=";",nrows=NumberRows,na.strings="",stringsAsFactors = FALSE)
Headers <- colnames(df)

### Corrections of the data
#replace BelowLimit with BelowLimitReplacement 
if (BelowLimitMethod == "SSR0"){
  df[,Headers[substr(Headers,1,10)=="Gentamicin"]] <-
    data.frame(lapply(df[,Headers[substr(Headers,1,10)=="Gentamicin"]], function(x) {
      gsub(BelowLimit, as.numeric(BelowLimitReplacement), x)
    }))
}

#replace BelowLimit with NA
if (BelowLimitMethod == "NA"){
  df[,Headers[substr(Headers,1,10)=="Gentamicin"]] <-
    data.frame(lapply(df[,Headers[substr(Headers,1,10)=="Gentamicin"]], function(x) {
      gsub(BelowLimit, NA, x)
    }))
}

#replace >12 with 12
df[,Headers[substr(Headers,1,10)=="Gentamicin"]] <-
  data.frame(lapply(df[,Headers[substr(Headers,1,10)=="Gentamicin"]], function(x) {
    gsub(">12", as.numeric(12), x)
  }))

#replace >9.9 with 9.9
df[,Headers[substr(Headers,1,10)=="Gentamicin"]] <-
  data.frame(lapply(df[,Headers[substr(Headers,1,10)=="Gentamicin"]], function(x) {
    gsub(">9.9", as.numeric(9.9), x)
  }))

#replace the Crea of <27umol/L with 27
df <- data.frame(lapply(df, function(x) {
  gsub("<27", as.numeric(0), x)
}))

#replace the Urea of <1.8mmol/L with 1.8
df <- data.frame(lapply(df, function(x) {
  gsub("<1.8", as.numeric(0), x)
}))


### Correct dates
#format="%d.%m.%Y --> the format the date is entered in the data
df$GEBURTSTAG <- as.Date(df$GEBURTSTAG, format="%d.%m.%Y")
df$KREATININ.DATUM <- as.Date(df$KREATININ.DATUM, format="%d.%m.%Y")
df$HARNSTOFF.DATUM <- as.Date(df$HARNSTOFF.DATUM, format="%d.%m.%Y")
df$DATUM.1 <- as.Date(df$DATUM.1, format="%d.%m.%Y")
df$DATUM.2 <- as.Date(df$DATUM.2, format="%d.%m.%Y")
df$DATUM.3 <- as.Date(df$DATUM.3, format="%d.%m.%Y")
df$DATUM.Gentanach.30.min <- as.Date(df$DATUM.Gentanach.30.min, format="%d.%m.%Y")
df$DATUM.Gentanach.4h <- as.Date(df$DATUM.Gentanach.4h, format="%d.%m.%Y")
df$DATUM.Gentanach.24h <- as.Date(df$DATUM.Gentanach.24h, format="%d.%m.%Y")

###correct time: 
#substrings in a character vector: substr(UHRZEIT.1,1,2) = first two numbers of the vector: hours --> /24 
#substr(df$UHRZEIT.1,4,5) = 4,5. number of the vector: minutes --> /60/24 
#this gives the time in days

#Alter
df$AgeDose0 <- df$DATUM.1 - df$GEBURTSTAG + as.numeric(substr(df$UHRZEIT.1,1,2))/24 +
  as.numeric(substr(df$UHRZEIT.1,4,5)) /60/24 #in days

meanAge <- mean(df$AgeDose0)

#Zeitpunkt der Dosen
df$TimeDose1 <- df$AgeDose0 - df$AgeDose0
df$TimeDose2 <- round(df$DATUM.2 - df$GEBURTSTAG + as.numeric(substr(df$UHRZEIT.2,1,2))/24 +
                        as.numeric(substr(df$UHRZEIT.2,4,5)) /60/24 - df$AgeDose0,4)
df$TimeDose3 <- round(df$DATUM.3 - df$GEBURTSTAG + as.numeric(substr(df$UHRZEIT.3,1,2))/24 +
                        as.numeric(substr(df$UHRZEIT.3,4,5)) /60/24 - df$AgeDose0,4)

#STOPZeitpunkt der Infusion
df$TimeDoseSTOP1 <- df$TimeDose1 + as.numeric(as.character(df$APPLIKATIONSDAUER))/60/24 #in days
df$TimeDoseSTOP2 <- df$TimeDose2 + as.numeric(as.character(df$APPLIKATIONSDAUER.2))/60/24
df$TimeDoseSTOP3 <- df$TimeDose3 + as.numeric(as.character(df$APPLIKATIONSDAUER.3))/60/24

#Zeitpunkt der Messungen
df$TimeGenta1 <- df$DATUM.Gentanach.30.min - df$GEBURTSTAG  +
  as.numeric(substr(df$UHRZEIT.Gentanach.30.min,1,2))/24 + 
  as.numeric(substr(df$UHRZEIT.Gentanach.30.min,4,5))/60/24 - df$AgeDose0
df$TimeGenta2 <- df$DATUM.Gentanach.4h - df$GEBURTSTAG  +
  as.numeric(substr(df$UHRZEIT.Gentanach.4h,1,2))/24 + 
  as.numeric(substr(df$UHRZEIT.Gentanach.4h,4,5))/60/24 - df$AgeDose0
df$TimeGenta3 <- df$DATUM.Gentanach.24h - df$GEBURTSTAG  +
  as.numeric(substr(df$UHRZEIT.Gentanach.24.h,1,2))/24 + 
  as.numeric(substr(df$UHRZEIT.Gentanach.24.h,4,5))/60/24 - df$AgeDose0

write.csv(df,"D:/BioMed IT CSV/CSV_Gentamicin\\Gentamicin_CORRdata.csv", row.names = FALSE)

## Repeat rows to have one row per measurement
maxTimes = 9 # 3 rows for Dosing, 3 rows stop of infusion, 3 rows for Measurement
df1 <- df[rep(seq_len(nrow(df)), each=maxTimes),]
rows = (1:nrow(df1))

## Generate individual rows for measurements
Nrow <- 1
rows1 <- rows[1:(Nrow+8)]==(Nrow+0) #Dose1
rows2 <- rows[1:(Nrow+8)]==(Nrow+1) #STOP
rows3 <- rows[1:(Nrow+8)]==(Nrow+2) #Dose2
rows4 <- rows[1:(Nrow+8)]==(Nrow+3) #STOP
rows5 <- rows[1:(Nrow+8)]==(Nrow+4) #Dose3 
rows6 <- rows[1:(Nrow+8)]==(Nrow+5) #STOP
rows7 <- rows[1:(Nrow+8)]==(Nrow+6) #Measurement1
rows8 <- rows[1:(Nrow+8)]==(Nrow+7) #Measurement2
rows9 <- rows[1:(Nrow+8)]==(Nrow+8) #Measurement3

#Time_I Infusiontime
df1$Time_I <- 0
df1$Time_I[rows[rep(rows1,nrow(df1)/maxTimes)]] <- df1$TimeDose1[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_I[rows[rep(rows2,nrow(df1)/maxTimes)]] <- df1$TimeDoseSTOP1[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_I[rows[rep(rows3,nrow(df1)/maxTimes)]] <- df1$TimeDose2[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_I[rows[rep(rows4,nrow(df1)/maxTimes)]] <- df1$TimeDoseSTOP2[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_I[rows[rep(rows5,nrow(df1)/maxTimes)]] <- df1$TimeDose3[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_I[rows[rep(rows6,nrow(df1)/maxTimes)]] <- df1$TimeDoseSTOP3[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_I[rows[rep(rows7,nrow(df1)/maxTimes)]] <- NA
df1$Time_I[rows[rep(rows8,nrow(df1)/maxTimes)]] <- NA
df1$Time_I[rows[rep(rows9,nrow(df1)/maxTimes)]] <- NA
df1$Time_I <- round(df1$Time_I*24,4) #in h

#Time_P Time of the measurement
df1$Time_P <- 0
df1$Time_P[rows[rep(rows1,nrow(df1)/maxTimes)]] <- NA
df1$Time_P[rows[rep(rows2,nrow(df1)/maxTimes)]] <- NA
df1$Time_P[rows[rep(rows3,nrow(df1)/maxTimes)]] <- NA
df1$Time_P[rows[rep(rows4,nrow(df1)/maxTimes)]] <- NA
df1$Time_P[rows[rep(rows5,nrow(df1)/maxTimes)]] <- NA
df1$Time_P[rows[rep(rows6,nrow(df1)/maxTimes)]] <- NA
df1$Time_P[rows[rep(rows7,nrow(df1)/maxTimes)]] <- df1$TimeGenta1[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_P[rows[rep(rows8,nrow(df1)/maxTimes)]] <- df1$TimeGenta2[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_P[rows[rep(rows9,nrow(df1)/maxTimes)]] <- df1$TimeGenta3[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time_P <- round(df1$Time_P*24,4)#in h

#Time combined in different rows
df1$Time <- 0
df1$Time[rows[rep(rows1,nrow(df1)/maxTimes)]] <- df1$TimeDose1[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows2,nrow(df1)/maxTimes)]] <- df1$TimeDoseSTOP1[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows3,nrow(df1)/maxTimes)]] <- df1$TimeDose2[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows4,nrow(df1)/maxTimes)]] <- df1$TimeDoseSTOP2[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows5,nrow(df1)/maxTimes)]] <- df1$TimeDose3[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows6,nrow(df1)/maxTimes)]] <- df1$TimeDoseSTOP3[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows7,nrow(df1)/maxTimes)]] <- df1$TimeGenta1[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows8,nrow(df1)/maxTimes)]] <- df1$TimeGenta2[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time[rows[rep(rows9,nrow(df1)/maxTimes)]] <- df1$TimeGenta3[rows[rep(rows1,nrow(df1)/maxTimes)]]
df1$Time <- round(df1$Time*24,4) #in h

#Conc
df1$Conc <- 0
df1$Conc[rows[rep(rows1,nrow(df1)/maxTimes)]] <- NA
df1$Conc[rows[rep(rows2,nrow(df1)/maxTimes)]] <- NA
df1$Conc[rows[rep(rows3,nrow(df1)/maxTimes)]] <- NA
df1$Conc[rows[rep(rows4,nrow(df1)/maxTimes)]] <- NA
df1$Conc[rows[rep(rows5,nrow(df1)/maxTimes)]] <- NA
df1$Conc[rows[rep(rows6,nrow(df1)/maxTimes)]] <- NA
df1$Conc[rows[rep(rows7,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$Gentamicin.Spiegel.nach.30.min[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Conc[rows[rep(rows8,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$Gentamicin.Spiegel.nach.4h[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Conc[rows[rep(rows9,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$Gentamicin.Spiegel.nach.24h[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Conc <- as.numeric(df1$Conc) #mg/L

#Dose
df1$Dose <- 0
df1$Dose[rows[rep(rows1,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$DOSIERUNG.1[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Dose[rows[rep(rows2,nrow(df1)/maxTimes)]] <- NA
df1$Dose[rows[rep(rows3,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$DOSIERUNG.2[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Dose[rows[rep(rows4,nrow(df1)/maxTimes)]] <- NA
df1$Dose[rows[rep(rows5,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$DOSIERUNG.3[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Dose[rows[rep(rows6,nrow(df1)/maxTimes)]] <- NA
df1$Dose[rows[rep(rows7,nrow(df1)/maxTimes)]] <- NA
df1$Dose[rows[rep(rows8,nrow(df1)/maxTimes)]] <- NA
df1$Dose[rows[rep(rows9,nrow(df1)/maxTimes)]] <- NA
df1$Dose <- as.numeric(format(round(as.numeric(df1$Dose), 4), nsmall = 2)) #mg

#AMT = Rate in infusions
df1$AMT <- 0 
df1$AMT[rows[rep(rows1,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$DOSIERUNG.1[rows[rep(rows1,nrow(df1)/maxTimes)]])) /
  as.numeric(as.character(df1$APPLIKATIONSDAUER[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$AMT[rows[rep(rows2,nrow(df1)/maxTimes)]] <- -1*(
  as.numeric(as.character(df1$DOSIERUNG.1[rows[rep(rows1,nrow(df1)/maxTimes)]])) /
  as.numeric(as.character(df1$APPLIKATIONSDAUER[rows[rep(rows1,nrow(df1)/maxTimes)]])))
df1$AMT[rows[rep(rows3,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$DOSIERUNG.2[rows[rep(rows1,nrow(df1)/maxTimes)]])) /
  as.numeric(as.character(df1$APPLIKATIONSDAUER.2[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$AMT[rows[rep(rows4,nrow(df1)/maxTimes)]] <- -1*(
  as.numeric(as.character(df1$DOSIERUNG.2[rows[rep(rows1,nrow(df1)/maxTimes)]])) /
  as.numeric(as.character(df1$APPLIKATIONSDAUER.2[rows[rep(rows1,nrow(df1)/maxTimes)]])))
df1$AMT[rows[rep(rows5,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$DOSIERUNG.3[rows[rep(rows1,nrow(df1)/maxTimes)]])) /
  as.numeric(as.character(df1$APPLIKATIONSDAUER.3[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$AMT[rows[rep(rows6,nrow(df1)/maxTimes)]] <- -1*(
  as.numeric(as.character(df1$DOSIERUNG.3[rows[rep(rows1,nrow(df1)/maxTimes)]])) /
  as.numeric(as.character(df1$APPLIKATIONSDAUER.3[rows[rep(rows1,nrow(df1)/maxTimes)]])))
df1$AMT[rows[rep(rows7,nrow(df1)/maxTimes)]] <- NA
df1$AMT[rows[rep(rows8,nrow(df1)/maxTimes)]] <- NA
df1$AMT[rows[rep(rows9,nrow(df1)/maxTimes)]] <- NA
df1$AMT <- round(df1$AMT*60, 4) #mg/h


#Infusiontime
df1$Inf_time <- 0
df1$Inf_time[rows[rep(rows1,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$APPLIKATIONSDAUER[rows[rep(rows1,nrow(df1)/maxTimes)]])) 
df1$Inf_time[rows[rep(rows2,nrow(df1)/maxTimes)]] <- NA
df1$Inf_time[rows[rep(rows3,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$APPLIKATIONSDAUER.2[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Inf_time[rows[rep(rows4,nrow(df1)/maxTimes)]] <- NA
df1$Inf_time[rows[rep(rows5,nrow(df1)/maxTimes)]] <- 
  as.numeric(as.character(df1$APPLIKATIONSDAUER.3[rows[rep(rows1,nrow(df1)/maxTimes)]]))
df1$Inf_time[rows[rep(rows6,nrow(df1)/maxTimes)]] <- NA
df1$Inf_time[rows[rep(rows7,nrow(df1)/maxTimes)]] <- NA
df1$Inf_time[rows[rep(rows8,nrow(df1)/maxTimes)]] <- NA
df1$Inf_time[rows[rep(rows9,nrow(df1)/maxTimes)]] <- NA 
df1$Inf_time <- as.numeric(format(round(as.numeric(df1$Inf_time)/60, 4), nsmall = 2)) #--> in min /60 = h

write.csv(df1,"D:/BioMed IT CSV/CSV_Gentamicin\\Gentamicin_dataframe.csv", row.names = FALSE)

### nlmixr
df2 <- data.frame(ID=df1$ï..Studien.ID, #Patient number 
                  INF_TIME=df1$Inf_time, #infusion duration
                  TIME_I=df1$Time_I, #Time of the infusion
                  TIME_P=df1$Time_P, # Time of the measurement
                  TIME=df1$Time, #combined time points
                  AMT=df1$AMT, #DOSE AMOUNT COLUMN: for Infusions = Rate; rate beginning and end (rate*-1)
                  DOSE=df1$Dose, #Dose
                  EVID=0, #Event identifier, --> Event will be defined later
                  DV=df1$Conc, #DEPENDENT VARIABLE: measurement = measured conc
                  lAGE=rep(log(as.numeric(as.character(df$AgeDose0))),each=9),
                  lKOF=rep(log(as.numeric(as.character(df$KOF.berechnet))), each=9),
                  lHEIGHT=rep(log(as.numeric(as.character(df$GROESSE))),each=9),
                  lWT=rep(log(as.numeric(as.character(df$GEWICHT))), each=9),
                  lCREA=rep(log(as.numeric(as.character(df$KREATININ))), each=9),
                  lUREA=rep(log(as.numeric(as.character(df$HARNSTOFF))), each=9),
                  lGA=rep(log(as.numeric(as.character(df$GA..Weeks.))), each=9),
                  SEX=rep(df$Sex, each=9)
)

df2$lCREA[df2$lCREA=="-Inf"] <- 0
df2$lUREA[df2$lUREA=="-Inf"] <- 0
df2$lGA[is.na(df2$lGA)] <- 0
df2$SEX <- as.numeric(df2$SEX) #male= 2, female=1

exam_grdat <- groupedData(DV ~ TIME | ID, df2)
plot(exam_grdat)

### nlmixr
df2$EVID[df2$CONC != 0] <- 0 # 0 if observation --> measurement
df2$EVID[df2$AMT != 0] <- 10101 #for i.v. 10000 +100*(Compartement Number) +1

#remove the NA rows from TIME --> no dosing and no measuring timepoint 
remove <- c()
for(i in seq(1,length(df2$ID))){
  if(all(is.na(df2[i,5]))){
    remove <- c(remove,i)
  }
}
df2 <- df2[-remove,]

#remove data with only 1 measurement
for(i in unique(df2$ID)){
  if(sum(is.na(df2[df2$ID==i,]$DV))<=1){
    df2 <- df2[-which(df2$ID==i),]
  }
}

#Check if time value is dublicated within the same ID 
#aggregate(df2$TIME, by=list(df2$ID), function(x) any(duplicated(x,incomparables = NA)))

#change duplicated time points --> add 0.016 h to change the time point
for(i in unique(df2$ID)){
  if(!all(duplicated(df2[df2$ID==i,]))){
    trans_dups <- duplicated(df2[df2$ID==i,]$TIME) | duplicated(df2[df2$ID==i,]$TIME, fromLast = T)
    trans_dups <- trans_dups & is.na(df2[df2$ID==i,]$AMT)
    df2[df2$ID==i,][trans_dups,]$TIME <- df2[df2$ID==i,][trans_dups,]$TIME - 0.016
  }
}

mixr_dat <- df2
mixr_dat[is.na(mixr_dat)] <- 0
mixr_dat = mixr_dat[order(mixr_dat$ID, mixr_dat$TIME),] 

# change time points which are too similar but not the same 
for(i in unique(mixr_dat$ID)){
  similar <- diff(rev(mixr_dat[mixr_dat$ID==i,]$TIME)) < 1e-03
  similar <- c(FALSE, similar)
  similar <- similar & mixr_dat[mixr_dat$ID==i,]$AMT!=0
  mixr_dat[mixr_dat$ID==i,][similar,]$TIME <- mixr_dat[df2$ID==i,][similar,]$TIME + 0.016
  
}

#abs(diff(mixr_dat$TIME))
#which.min(abs(diff(mixr_dat$TIME)))
#min(abs(diff(mixr_dat$TIME)))

mixr_dat = mixr_dat[order(mixr_dat$ID, mixr_dat$TIME),] # change of time points: new order needed

mixr_dat$AMT <- as.integer(mixr_dat$AMT)
mixr_dat$EVID <- as.integer(mixr_dat$EVID)
mixr_dat$ID <- as.numeric(as.character(mixr_dat$ID))

write.csv(mixr_dat,"D:/BioMed IT CSV/CSV_Gentamicin\\Gentamicin_dataframe.csv", row.names = FALSE)



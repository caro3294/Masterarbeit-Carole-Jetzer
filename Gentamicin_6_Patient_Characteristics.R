# Gentamicin Patient Characteristics
# autor: Carole Jetzer
# date: 21.05.2020

rm(list=ls()) 
par.default <- par(no.readonly = TRUE) 
setwd("C:/Users/carol/OneDrive/Studium/ETH/Master/Masterarbeit/R") 

df1 <- read.csv("D:/BioMed IT CSV/CSV_Gentamicin\\Gentamicin_CORRdata.csv", sep=",")

df <- data.frame(ID=df1$ï..Studien.ID, #Patient number 
                AGE=df1$AgeDose0,
                KOF=df1$KOF.berechnet,
                HEIGHT=df1$GROESSE,
                WT=df1$GEWICHT,
                CREA=df1$KREATININ,
                UREA=df1$HARNSTOFF,
                GA=df1$GA..Weeks.,
                SEX=df1$Sex,
                IND=df1$DIAGNOSE_GROUPED
                )

df$CREA[df$CREA==27] <- 1 #lower 27
df$CREA[df$CREA!=1] <- 2 #higher than 27

df$GA[df$lGA==0] <- NA
df$SEX <- as.numeric(df$SEX) #male= 2, female=1

df$IND <- as.numeric(df$IND) #Viral Infection = 4, Bacterial Infection = 1, Respiratory Disorder = 2,3

df$IND[df$IND==2] <- 3 #Respiratory Disorder
df$IND[df$IND==1] <- 2 #Bacterial Infection
df$IND[df$IND==4] <- 1 #Viral Infection

Genta.Charact <- data.frame(
                    MeanAge=mean(df$AGE),
                    sdAGE=sqrt(var(df$AGE)),
                    MeanKOF=mean(df$KOF),
                    sdKOF=sqrt(var(df$KOF)),
                    MeanHEIGHT=mean(df$HEIGHT),
                    sdHEIGHT=sqrt(var(df$HEIGHT)),
                    MeanWT=mean(df$WT),
                    sdWT=sqrt(var(df$WT)),
                    MeanUREA=mean(df$UREA, na.rm=T),
                    sdUREA=sqrt(var(df$UREA, na.rm=T)),
                    n_CREA_M27=sum(df$CREA[df$CREA==2]),
                    n_CREA_U27=sum(df$CREA[df$CREA==1]),
                    n_MALE=sum(df$SEX[df$SEX==2]),
                    n_FEMALE=sum(df$SEX[df$SEX==1]),
                    n_Viral=sum(df$IND[df$IND==1]),
                    n_Bacterial=sum(df$IND[df$IND==2]),
                    n_Resp=sum(df$IND[df$IND==3])
                    )
print(Genta.Charact)

Percent <- 42/73*100

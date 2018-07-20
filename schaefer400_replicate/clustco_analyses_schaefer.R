library(dplyr)
library(psych)
library(mgcv)
library(car)
library(ggplot2)
library(visreg)
library(sjPlot)
library(parallel)
library(RLRsim)
library(gamm4)
library(effects)
library(lm.beta)
require(gridExtra)
require(grid)
library(gridBase)
library(ggthemes)
#set ggplot theme
theme_set(theme_few())
#unless needed, will mask dplyr select
library(MASS)
library(ppcor)

# SETUP -------------------------------------------------------------------
#Cluster
setwd("/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/")
#Local personal computer
setwd("~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/")
subinfodir="~/Documents/bassett_lab/tooleyEnviNetworks/data/subjectData/"
sublistdir="~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists/"
qadir="~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/"
clustcodir="~/Dropbox (Personal)/bassett_lab/clustco_paper/"
analysis_dir="~/Documents/bassett_lab/tooleyEnviNetworks/analyses/"

#on bassett computer
setwd("~/Documents/tooleyEnviNetworks/data/rest/")
subinfodir="~/Documents/tooleyEnviNetworks/data/subjectData/"
sublistdir="~/Documents/tooleyEnviNetworks/subjectLists/"
qadir="~/Documents/tooleyEnviNetworks/data/rest/"

#subjlist<-read.csv("/data/joy/BBL/projects/tooleyEnviNetworks/subjectLists/n1015_healthT1RestExclude.csv")
subjlist<-read.csv(paste0(sublistdir,"n1015_healthT1RestExclude.csv"))
#use LTN exclude critera
#subjlist<-read.csv(paste0(sublistdir,"n885_LTNexclude.csv"))

# Get demographics to control for
file1<-read.csv(paste0(subinfodir, "n1601_demographics_go1_20161212.csv"))
file2<-read.csv(paste0(subinfodir, "n1601_go1_environment_factor_scores_tymoore_20150909.csv"))
#Get QA values to include in analyses
file3<-read.csv(paste0(qadir, "n1601_RestQAData_20170509.csv"))
#Get the network stats file with modul and clustering in schaefer
file4<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1015_sub_net_meas_schaefer_signed.csv")

#rename the second column of the network statistics file to be 'scanid' so it matches below
file4<-dplyr::rename(file4, scanid=subjlist_2)

#Join all files together
master<-right_join(file1, subjlist, by ="scanid")
master<-right_join(file2,master, by="scanid")
master<-right_join(file3, master, by= "scanid")
master<-right_join(file4, master, by ="scanid")
master_null1<-right_join(file5, master, by ="scanid")
master_null2<-right_join(file6, master, by ="scanid")

#eliminate extraneous columns
master<-master %>% dplyr::select(., -c(restExcludeVoxelwise, restNoDataExclude))
master<-master %>% dplyr::select(., -c(restRelMeanRMSMotionExclude:restRpsMapCorrectionNotApplied))
master<-master %>% dplyr::select(., -c(bblid.y, bblid.x, bblid.y.y, subjlist_1))

#make sure factor variables are factors
master$race<-factor(master$race)
master$race2<-factor(master$race2, labels=c("White", "Black", "Other"))
master$sex<-factor(master$sex, labels=c("Male", "Female"))

#make an age squared variable
master$ageatscansqdem <- (master$ageAtScan1-mean(master$ageAtScan1))^2
master$ageatscansq <- (master$ageAtScan1)^2
master$ageAtScan1yrs<-(master$ageAtScan1)/12
summary(master$envSES)
#split SES on the median 
master$envSEShigh=NA
master$envSEShigh[master$envSES >= 0.0178] <- 1
master$envSEShigh[master$envSES < 0.0178]<- 0
master$envSEShigh<-factor(master$envSEShigh, labels=c("Low", "High"), ordered=TRUE)
summary(master$envSEShigh)

#remove the outlier subj 3815
#master<-master[master$scanid!=3815,]
#should potentially be centering the age and envSES variables before looking at interactions
master$ageAtScan1cent<-(master$ageAtScan1-mean(master$ageAtScan1))
master$ageAtScan1yrscent<-(master$ageAtScan1yrs-mean(master$ageAtScan1yrs))
master$envSEScent<-(master$envSES-mean(master$envSES))
master$medu1cent=(master$medu1-mean(master$medu1[!is.na(master$medu1)]))
#create parental education variable
master$paredu1 <- master$medu1+master$fedu1
master$paredu1cent=(master$paredu1-mean(master$paredu1[!is.na(master$paredu1)]))

# Clustering ----------------------------------------------------------------
#linear age effect
l <- lm(avgclustco_both ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh, data=master)
summary(l)
lm.beta(l)
lmageplot<-visreg(l, "ageAtScan1yrs",
                  main="Average Clustering Coefficient", xlab="Age in Years (centered)", 
                  ylab="Average Clustering Coefficient (partial residuals)",line=list(col="black",lwd=8))

#linear model med split
l <- lm(avgclustco_both~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l)
lm.beta(l)
visreg(l, "ageAtScan1yrs", by="envSEShigh", main="Mean Clustering Coefficient (partial residuals)",
       xlab="Age in Years", ylab="", overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(.18, .24), legend=FALSE,
       line=list(col=c(rgb(28, 147, 255, maxColorValue = 255), rgb(255, 168, 28, maxColorValue = 255))), 
       fill=list(col=c(alpha(rgb(28, 147, 255, maxColorValue = 255), 0.7), alpha(rgb(255, 168, 28, maxColorValue = 255),0.7))),
       strip.names=c("Low SES", "High SES"))

#test for nonlinear relationship between age and clustco
ageonlyRlrtmodel<-gamm(avgclustco_both~ageAtScan1cent+sex+race2+s(avgweight)+restRelMeanRMSMotion, method='REML', data=master)$lme
l<-exactRLRT(ageonlyRlrtmodel) #no non-linear relationship
l


# Modularity --------------------------------------------------------------
#linear age effect
l <- lm(modul ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh, data=master)
summary(l)
lm.beta(l)
lmageplot<-visreg(l, "ageAtScan1yrs",
                  main="Average Clustering Coefficient", xlab="Age in Years (centered)", 
                  ylab="Average Clustering Coefficient (partial residuals)",line=list(col="black",lwd=8))

#linear model med split
l <- lm(modul~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l)
lm.beta(l)
visreg(l, "ageAtScan1yrs", by="envSEShigh", main="Mean Clustering Coefficient (partial residuals)",
       xlab="Age in Years", ylab="", overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(.18, .25), legend=FALSE,
       line=list(col=c(rgb(28, 147, 255, maxColorValue = 255), rgb(255, 168, 28, maxColorValue = 255))), 
       fill=list(col=c(alpha(rgb(28, 147, 255, maxColorValue = 255), 0.7), alpha(rgb(255, 168, 28, maxColorValue = 255),0.7))),
       strip.names=c("Low SES", "High SES"))

#include modularity in clust co and vice-versa
l <- lm(avgclustco_both ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+modul+ envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l)

l <- lm(modul~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+avgclustco_both+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l)

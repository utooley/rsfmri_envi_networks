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
library(lmtest)
#not until needed, will mask dplyr select
library(MASS)
library(ppcor)

###############
#### SETUP ####
###############
#get files
#Local personal computer
setwd("~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/")
subinfodir="~/Documents/projects/in_progress/tooleyEnviNetworks/data/subjectData/"
sublistdir="~/Documents/projects/in_progress/tooleyEnviNetworks/subjectLists/"
qadir="~/Documents/projects/in_progress/tooleyEnviNetworks/data/rest/"
clustcodir="~/Dropbox (Personal)/projects/in_progress/tooleyenvinetworks/code/clustco_paper/"
analysis_dir="~/Documents/projects/in_progress/tooleyEnviNetworks/analyses/"

#Local mackey computer
setwd("~/Documents/projects/in_progress/tooleyEnviNetworks/data/rest/")
subinfodir="~/Documents/projects/in_progress/tooleyEnviNetworks/data/subjectData/"
sublistdir="~/Documents/projects/in_progress/tooleyEnviNetworks/subjectLists/"
qadir="~/Documents/projects/in_progress/tooleyEnviNetworks/data/rest/"
clustcodir="~/Dropbox/projects/in_progress/tooleyEnviNetworks/code/clustco_paper/"
analysis_dir="~/Documents/projects/in_progress/tooleyEnviNetworks/analyses/"
reho_dir="~/Documents/projects/in_progress/tooleyEnviNetworks/data/rest/"

#get subjlist
subjlist<-read.csv(paste0(sublistdir,"n1012_healthT1RestExclude_parcels.csv"))

# Get demographics to control for
file1<-read.csv(paste0(subinfodir, "n1601_demographics_go1_20161212.csv"))
file2<-read.csv(paste0(subinfodir, "n1601_go1_environment_factor_scores_tymoore_20150909.csv"))
#Get QA values to include in analyses
file3<-read.csv(paste0(qadir, "n1601_RestQAData_20170509.csv"))
#Get the network stats file with the different clustering coefficients
file4<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1012_net_meas_normed3clustcoefs_signed.csv")
#get the original file with modularity in it
file5<-read.csv("~/Documents/projects/in_progress//tooleyEnviNetworks/analyses/n1012_sub_net_meas_signed.csv")
#rename the second column of the network statistics file to be 'scanid' so it matches the other files
file4<-dplyr::rename(file4, scanid=subjlist_2)
file5<-dplyr::rename(file5, scanid=subjlist_2)

#Join all files together
master<-right_join(file1, subjlist, by ="scanid")
master<-right_join(file2,master, by="scanid")
master<-right_join(file3, master, by= "scanid")
master<-right_join(file4, master, by ="scanid")
master<-right_join(file5, master, by ="scanid")

##############
#### Data Cleaning ###
##############

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
#should be centering the age and envSES variables before looking at interactions
master$ageAtScan1cent<-(master$ageAtScan1-mean(master$ageAtScan1))
master$ageAtScan1yrscent<-(master$ageAtScan1yrs-mean(master$ageAtScan1yrs))
#split SES on the median 
summary(master$envSES)
master$envSEShigh=NA
master$envSEShigh[master$envSES >= 0.0178] <- 1
master$envSEShigh[master$envSES < 0.0178]<- 0
master$envSEShigh<-factor(master$envSEShigh, labels=c("Low", "High"), ordered=TRUE)
#center
master$envSEScent<-(master$envSES-mean(master$envSES))
master$medu1cent=(master$medu1-mean(master$medu1[!is.na(master$medu1)]))
#create parental education variable
master$paredu1 <- master$medu1+master$fedu1
master$paredu1cent=(master$paredu1-mean(master$paredu1[!is.na(master$paredu1)]))

#diff in age between groups?
t.test(master$ageAtScan1yrs~master$envSEShigh)

##############
#### SECTION 1: CONTINUOUS ####
##############

#linear model continuous (with demeaned variables)
l <- lm(avgclustco_both ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion+envSEScent+ageAtScan1cent*envSEScent, data=master)
summary(l)
#get betas
lm.beta(l)
#compare to scaling all variables
l2 <- lm(scale(avgclustco_both) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(envSEScent)+scale(ageAtScan1cent)*scale(envSEScent), data=master)
summary(l2)
#compare to a model without the interaction
l3 <- lm(avgclustco_both ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion+envSES, data=master)
summary(l3)
visreg(l, "ageAtScan1cent", by="envSEScent", main="Linear Age*SES Int with Median Split SES")
anova(l3, l, test="Chisq")#compare the med split interaction with an age-only model
#use the likelihood ration test to compare
lrtest(l3, l)

##linear model continuous for modularity
l <- lm(modul ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion+envSEScent+ageAtScan1cent*envSEScent, data=master)
summary(l) # no significant interaction
l3 <- lm(modul ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion, data=master)
summary(l3)
anova(l3, l, test="Chisq")#compare the med split interaction with an age-only model
#use the likelihood ration test to compare
lrtest(l3, l)

##############
#### SECTION 2: MAT ED ####
##############
shapiro.test(master$envSES)
shapiro.test(master$medu1cent)
cor.test(master$envSES,master$medu1, method = "spearman", alternative = "two.sided")

#MATERNAL EDUCATION
#effects still hold with maternal ed added to the model
lmat <- lm(avgclustco_both ~ ageAtScan1cent+sex+race2+avgweight+envSEShigh+medu1cent+restRelMeanRMSMotion, data=master)
summary(lmat)
lm.beta(lmat)
lmat <- lm(avgclustco_both ~ ageAtScan1cent+sex+race2+avgweight+medu1cent+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh, data=master)
summary(lmat)
lm.beta(lmat)
#use maternal ed instead of envSES
lmat <- lm(modul ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion+medu1cent+ageAtScan1cent*medu1cent, data=master)
summary(lmat)

#############
##### STRINGENT EXCLUSION CRITERIA ###
############
# see 12_tooley_rsfmri_nets_ses_supplement_smaller_sample.R

#############
##### NEIGHBORHOOD SES MODELED CONTINUOUSLY ###
############
# see 12_tooley_rsfmri_nets_ses_supplement_cont_ses.R

#############
##### STRINGENT EXCLUSION CRITERIA + NEIGHBORHOOD SES MODELED CONTINUOUSLY ###
############
# see 12_tooley_rsfmri_nets_ses_supplement_smaller_sample_cont_ses.R

#############
##### BIVARIATE RELATIONSHIPS BETWEEN PREDICTORS ###
############
library(Hmisc)
mytable <- data.frame(master$ageAtScan1yrs, master$avgweight.x, master$restRelMeanRMSMotion, master$envSES, master$avgclustco_both.x, master$modul)
sapply(mytable, function(x) {shapiro.test(x)})
t <- rcorr(as.matrix(mytable), type="spearman")
View(t$r)

#############
##### REPLICATION IN ALTERNATIVE PARCELLATION ###
############
#see schaefer400_replicate

##############
#### THRESHOLDED NETS/REMOVAL OF NEGATIVE EDGE WEIGHTS ####
##############

# ZHANG-HORVATH FORMULA, POSITIVE WEIGHTS ONLY #
#get net stats calculated on thresholded matrices
file4<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1012_sub_net_meas_signed_zhang_horvath.csv")

# ONNELA FORMULA, POSITIVE WEIGHTS ONLY #
file4<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1012_sub_net_meas_signed_onnela.csv")

file4<-dplyr::rename(file4, scanid=subjlist_2)
master<-right_join(file1, subjlist, by ="scanid")
master<-right_join(file2,master, by="scanid")
master<-right_join(file3, master, by= "scanid")
master<-right_join(file4, master, by ="scanid")

### RUN DATA CLEANING SECTION AT THE TOP HERE ##

#linear age effect without interaction
l <- lm(avgclco_pos ~ ageAtScan1yrs+sex+race2+restRelMeanRMSMotion+avgweight+envSEShigh, data=master)
summary(l)
l.beta<- lm.beta(l)
lmageplot<-visreg(l, "ageAtScan1yrs",
                  main="Average Clustering Coefficient", xlab="Age in Years (centered)", ylab="Average Clustering Coefficient (partial residuals)")

#linear model med split
l2 <- lm(avgclco_pos ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l2)
l2.beta <- lm.beta(l2)
anova(l,l2,test="Chisq")
l2 <- lm(scale(avgclco_pos) ~ scale(ageAtScan1yrs)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+envSEShigh+scale(ageAtScan1yrs)*envSEShigh, data=master)
summary(l2)

#linear model continuous
l2 <- lm(avgclco_pos ~ ageAtScan1yrscent+sex+race2+avgweight+restRelMeanRMSMotion+envSES, data=master)
summary(l2)
l2.beta <- lm.beta(l2)
anova(l,l2,test="Chisq")
######MODULARITY######

#linear age effect without interaction
l <- lm(modul ~ ageAtScan1yrs+sex+race2+avgweight+envSEShigh+restRelMeanRMSMotion, data=master)
summary(l)
lm.beta(l)
lmageplot<-visreg(l, "ageAtScan1yrs",
                  main="Average Clustering Coefficient", xlab="Age in Years (centered)", ylab="Average Clustering Coefficient (partial residuals)")

#linear model med split
l2 <- lm(modul ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l2)
lm.beta(l2)

##############
#### NETWORK SEGREGATION PER CHAN ET AL. 2018 ####
##############
#### SIGNED MATRICES ####
#see also file 02_net_meas_for_subjs_signed.m for input

#read in file, add to rest of data
file7 <- read.csv("~/Dropbox (Personal)/projects/in_progress/tooleyenvinetworks/data/n1012_sub_net_meas_signed_w_segregation.csv")
file7 <- read.csv("~/Dropbox/projects/in_progress/tooleyenvinetworks/data/n1012_sub_net_meas_signed_w_segregation.csv")
file7<-dplyr::rename(file7, scanid=subjlist_2)
master <- master %>% dplyr::select(.,- (avgweight:modul)) #remove first set of these variables so names don't overlap, already checked that they are exactly the same
master<-right_join(file7, master, by ="scanid")

#is system segregation correlated with modularity or avg clustco?
cor.test(master$modul,master$sys_segreg,method = "spearman")
cor.test(master$avgclustco_both,master$sys_segreg,method = "spearman")

#when controlling for other variables?
temp<-cbind(master$ageAtScan1cent, as.factor(master$sex), as.factor(master$race2), master$avgweight, master$restRelMeanRMSMotion, as.factor(master$envSEShigh))
library(ppcor)
pcor.test(master$modul, master$sys_segreg, temp, method = "spearman")
pcor.test(master$avgclustco_both, master$sys_segreg, temp, method = "spearman")

#linear age effect without interaction
l <- lm(sys_segreg ~ ageAtScan1yrs+sex+race2+avgweight+envSEShigh+restRelMeanRMSMotion, data=master)
summary(l)
lm.beta(l)

#linear interaction model
l2 <- lm(sys_segreg ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l2)
lm.beta(l2)

#### THRESHOLDED MATRICES ####
#is pos system segregation correlated with modularity or avg clustco?
cor.test(master$modul,master$sys_segreg_posonly,method = "spearman")
cor.test(master$avgclustco_both,master$sys_segreg_posonly,method = "spearman")

#when controlling for other variables?
temp<-cbind(master$ageAtScan1cent, as.factor(master$sex), as.factor(master$race2), master$avgweight, master$restRelMeanRMSMotion, as.factor(master$envSEShigh))
library(ppcor)
pcor.test(master$modul, master$sys_segreg_posonly, temp, method = "spearman")
pcor.test(master$avgclustco_both, master$sys_segreg_posonly, temp, method = "spearman")

#linear age effect without interaction
l <- lm(sys_segreg_posonly ~ ageAtScan1yrs+sex+race2+avgweight+envSEShigh+restRelMeanRMSMotion, data=master)
summary(l)
lm.beta(l)

#linear interaction model med split
l2 <- lm(sys_segreg_posonly ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l2)
lm.beta(l2)

#############################################
#### PARTICIPATION COEFFICIENT ANALYSES ####
############################################

#Get the network stats file with all the different clustering coefficients
file4<-read.csv("~/Dropbox (Personal)/projects/in_progress/tooleyEnviNetworks/data/n1012_sub_net_meas_signed_partcoefficient.csv")

#rename the second column of the network statistics file to be 'scanid' so it matches below
file4<-dplyr::rename(file4, scanid=subjlist_2)

#Join all files together
master<-right_join(file1, subjlist, by ="scanid")
master<-right_join(file2,master, by="scanid")
master<-right_join(file3, master, by= "scanid")
master<-right_join(file4, master, by ="scanid")

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
#split SES on the median 
master$envSEShigh=NA
master$envSEShigh[master$envSES >= 0.0178] <- 1
master$envSEShigh[master$envSES < 0.0178]<- 0
master$envSEShigh<-factor(master$envSEShigh, labels=c("Low", "High"), ordered=TRUE)

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

#add participation coefficients for positive and negative weights together
master$yeo_part_coef_both <- master$yeo_part_coef_neg+master$yeo_part_coef_pos

temp <- master %>% select(., yeo_part_coef_both, yeo_part_coef_pos, yeo_part_coef_neg, ageAtScan1yrscent, envSEScent, medu1cent)
measures <- temp %>% names()
#look at age effect here
for (meas in measures){
  name<-paste0("lm_age_",meas)
  formula<-formula(paste0(meas, '~ ageAtScan1yrs+sex+race2+avgweight+envSES+restRelMeanRMSMotion'))
  assign(name, lm(formula, data=master))
}
summary(lm_age_yeo_part_coef_both)
summary(lm_age_yeo_part_coef_pos)
summary(lm_age_yeo_part_coef_neg)

#look at age effect in context of maternal education
for (meas in measures){
  name<-paste0("lm_age_",meas)
  formula<-formula(paste0(meas, '~ ageAtScan1yrscent+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh*ageAtScan1yrscent'))
  assign(name, lm(formula, data=master))
}
summary(lm_age_yeo_part_coef_both)
summary(lm_age_yeo_part_coef_pos)
summary(lm_age_yeo_part_coef_neg)

##############
#### SUPP FIGURE 1 ####
##############

#linear model med split
l <- lm(avgclustco_both ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l)
visreg(l, "ageAtScan1yrs", by="envSEShigh", main="Mean Clustering Coefficient (partial residuals)",
       xlab="Age in Years", ylab="", partial=TRUE, rug=FALSE, legend=FALSE,
       line=list(col=c(rgb(28, 147, 255, maxColorValue = 255), rgb(255, 168, 28, maxColorValue = 255))), 
       fill=list(col=c(alpha(rgb(28, 147, 255, maxColorValue = 255), 0.7), alpha(rgb(255, 168, 28, maxColorValue = 255),0.7))),
       points=list(col=c(rgb(28, 147, 255, maxColorValue = 255), rgb(255, 168, 28, maxColorValue = 255))), 
       strip.names=c("Low SES", "High SES"))

#alternative getting residuals with rstandard
l <- lm(avgclustco_both ~ sex+race2+avgweight+restRelMeanRMSMotion, data=master)
residdf <- master
residdf$resid <- rstandard(l)
shapiro.test(residdf$resid)
residdf$mean<-mean(residdf$avgclustco_both)
residdf$realresid<-residdf$mean+(residdf$resid)
l<-lm(realresid~ageAtScan1yrs*envSEShigh, data=residdf)
summary(l)
visreg(l, "ageAtScan1yrs",by="envSEShigh",ylab="", xlab="SES Composite", overlay=TRUE)
scatterplot(realresid~ageAtScan1yrs|envSEShigh, data=residdf, boxplots=FALSE,
            col=palette()[c(2,1)],legend.title = "SES", grid=FALSE, pch=c(19,19),
            main="", xlab="Age", ylab="ReHo", cex.main=1.5, cex.lab=1.2)

##############
#### SUPP FIGURE 2 ####
##############
full_nodewise_clustco<-read.csv(paste0(clustcodir,"n1012_clust_co_nodewise_by_subj.csv"))
full_nodewise_clustco<-dplyr::rename(full_nodewise_clustco, scanid=subjlist_2)
full_nodewise_clustco<-right_join(full_nodewise_clustco,master,by="scanid")
a<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/parcels/Glasser_to_Yeo.csv")
a$id<-gsub("^NZMod_", "", a$id)
#remove parcel R52
a<-a[-103,]
#get a dataframe with each subject with the mean clust co in each system
full_nodewise_clustco<-select(full_nodewise_clustco, -(avgweight:paredu1cent))
b<-t(full_nodewise_clustco)
colnames(b)<-b[2,]
b<-b[-(1:2),]
l<-cbind(a,b)
subject_clustco_yeo_system<-l %>% group_by(Yeo_Parcellation.7YeoPNC.nii.gz_0...) %>% summarize_all(mean)
#transpose it
subject_clustco_yeo_system<-t(subject_clustco_yeo_system)
subject_clustco_yeo_system<-subject_clustco_yeo_system[-(1:2),]
#make the index the rownames
subject_clustco_yeo_system <- cbind(scanid = rownames(subject_clustco_yeo_system), subject_clustco_yeo_system)
rownames(subject_clustco_yeo_system) <- 1:nrow(subject_clustco_yeo_system)
colnames(subject_clustco_yeo_system)<-c("scanid", "Yeo_1", "Yeo_2", "Yeo_3", "Yeo_4", "Yeo_5", "Yeo_6", "Yeo_7")
subject_clustco_yeo_system<-data.frame(subject_clustco_yeo_system)
#now can match it up to master and look at each the effect for each system
subject_clustco_yeo_system[,c(1:8)]<-sapply(subject_clustco_yeo_system[,c(1:8)], as.character)
subject_clustco_yeo_system[,c(1:8)]<-sapply(subject_clustco_yeo_system[,c(1:8)], as.numeric)
master<-right_join(subject_clustco_yeo_system, master, by ="scanid")

#FOR AGE X SES INTERACTION ACROSS COG SYSTEMS
par(mfrow=c(2,4))#graphs for Allyson

#Yeo 1
Yeo1_scaled <- lm(scale(Yeo_1) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo1_scaled)
age_ses_beta_yeo1<-lm.beta(Yeo1)$standardized.coefficients[9]
##Plot
Yeo1<- lm(Yeo_1 ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
summary(Yeo1)
visreg(Yeo1, "ageAtScan1yrs", by="envSEShigh", overlay=TRUE, xlab="Yeo_1", ylim=c(0.13,.35), ylab="", legend=FALSE)
interact_plot(Yeo1, pred = "ageAtScan1cent", modx = "envSEShigh")
plot(effect(term="ageAtScan1yrs:envSEShigh", mod=Yeo1,default.levels=20),multiline=TRUE);

#Yeo 2
Yeo2_scaled<- lm(scale(Yeo_2) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo2_scaled)
age_ses_beta_yeo2<-lm.beta(Yeo2)$standardized.coefficients[9]
Yeo2<- lm(Yeo_2 ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
visreg(Yeo2, "ageAtScan1yrs", by="envSEShigh", overlay=TRUE, xlab="Yeo_2", ylim=c(0.13,.35),ylab="",legend=FALSE)
interact_plot(Yeo2, pred = "ageAtScan1yrs", modx = "envSEShigh")
#Yeo 3
Yeo3_scaled<- lm(scale(Yeo_3) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo3_scaled)
age_ses_beta_yeo3<-lm.beta(Yeo3)$standardized.coefficients[9]
Yeo3<- lm(Yeo_3 ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
visreg(Yeo3, "ageAtScan1yrs", by="envSEShigh", overlay=TRUE, xlab="Yeo_3", ylim=c(0.13,.35),ylab="",legend=FALSE)
interact_plot(Yeo3, pred = "ageAtScan1cent", modx = "envSEShigh")
#Yeo 4
Yeo4_scaled<- lm(scale(Yeo_4) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo4_scaled)
age_ses_beta_yeo4<-lm.beta(Yeo4)$standardized.coefficients[9]
Yeo4<- lm(Yeo_4 ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
visreg(Yeo4, "ageAtScan1yrs", by="envSEShigh", overlay=TRUE, xlab="Yeo_4", ylim=c(0.13,.35),ylab="",legend=FALSE)
interact_plot(Yeo4, pred = "ageAtScan1cent", modx = "envSEShigh")
#Yeo 5
Yeo5_scaled<- lm(scale(Yeo_5) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo5_scaled)
age_ses_beta_yeo5<-lm.beta(Yeo5)$standardized.coefficients[9]
Yeo5<- lm(Yeo_5 ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
visreg(Yeo5, "ageAtScan1yrs", by="envSEShigh", overlay=TRUE, xlab="Yeo_5", ylim=c(0.13,.35),ylab="",legend=FALSE)
interact_plot(Yeo5, pred = "ageAtScan1yrs", modx = "envSEShigh")
#Yeo 6
Yeo6_scaled <- lm(scale(Yeo_6) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo6_scaled)
age_ses_beta_yeo6<-lm.beta(Yeo6)$standardized.coefficients[9]
Yeo6<- lm(Yeo_6 ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
visreg(Yeo6, "ageAtScan1yrs", by="envSEShigh", overlay=TRUE, xlab="Yeo_6", ylim=c(0.13,.35),ylab="",legend=FALSE)
interact_plot(Yeo6, pred = "ageAtScan1yrs", modx = "envSEShigh")
#Yeo 7
Yeo7_scaled <- lm(scale(Yeo_7) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo7_scaled)
age_ses_beta_yeo7<-lm.beta(Yeo7)$standardized.coefficients[9]
Yeo7<- lm(Yeo_7 ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
visreg(Yeo7, "ageAtScan1yrs", by="envSEShigh", overlay=TRUE, xlab="Yeo_7", ylim=c(0.13,.35),ylab="",legend=FALSE)
interact_plot(Yeo7, pred = "ageAtScan1yrs", modx = "envSEShigh")

##############
#### SUPP SECTION 2: REHO ####
##############
#Get the reho file by subject
file4<-read.csv(paste0(rehodir,"n1601_glasserReHoValues_20170509.csv"))
master<-right_join(file4, master, by ="scanid")
#remove parcel R52 which is missing signal
master<-select(master, -rest_glasser_reho_Right_52)
#average whole-brain reho for each subject
master <- master %>% rowwise() %>% mutate(mean_reho=mean(rest_glasser_reho_Right_V1:rest_glasser_reho_Left_p24))
summary(master$mean_reho)

l <- lm(mean_reho~ageAtScan1cent+race2+sex+restRelMeanRMSMotion+envSEShigh, data=master)
summary(l)
lm.beta(l)
l <- lm(scale(mean_reho)~scale(ageAtScan1cent)+race2+sex+scale(restRelMeanRMSMotion)+envSEShigh+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(l)
l <- lm(mean_reho~ageAtScan1cent+race2+sex+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh, data=master)
summary(l)
lm.beta(l)

###############
#### SUPP FIGURE 4 ###
###############
#plot the across regions reho by clust co correlation, controlling for subject level covariates
nodewise_clustco <- select(full_nodewise_clustco, avgclustco_both_..1:avgclustco_both_359)
nodewise_reho <- select(master, rest_glasser_reho_Right_V1:rest_glasser_reho_Left_p24)
resid_reho_perregion=numeric(359)
for (i in 1:359){
  clustco=nodewise_clustco[i]
  reho=nodewise_reho[i]
  temp=cbind(master$ageAtScan1cent,as.factor(master$sex), as.factor(master$race2), master$restRelMeanRMSMotion,master$envSEShigh, clustco, reho)
  colnames(temp)=c("ageAtScan1cent", "sex", "race2", "restRelMeanRMSMotion", "envSEShigh", "clustco", "reho")
  residualreho <- rstandard(lm(reho~ageAtScan1cent+race2+sex+restRelMeanRMSMotion+envSEShigh, data=temp))
  residualreho <- mean(temp$reho)+residualreho
  #for each region, average across subjects
  resid_reho_perregion[[i]] <- mean(residualreho)
  #correlationspvals[[i]] <- as.numeric(pcor.test(clustco, reho, temp)$p.value)
}
plot(mean_nodes_clustco,resid_reho_perregion, col="darkgreen", pch=19)
cor.test(mean_nodes_clustco, resid_reho_perregion)
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
library(stargazer)
#not until needed, will mask dplyr select
library(MASS)
library(ppcor)

###############
#### SETUP ####
###############
#get files
#Local personal computer
setwd("~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/")
subinfodir="~/Documents/bassett_lab/tooleyEnviNetworks/data/subjectData/"
sublistdir="~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists/"
qadir="~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/"
clustcodir="~/Dropbox (Personal)/bassett_lab/clustco_paper/"
analysis_dir="~/Documents/bassett_lab/tooleyEnviNetworks/analyses/"
reho_dir="~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/"

#get subjlist
subjlist<-read.csv(paste0(sublistdir,"n1012_healthT1RestExclude_parcels.csv"))

# Get demographics to control for
file1<-read.csv(paste0(subinfodir, "n1601_demographics_go1_20161212.csv"))
file2<-read.csv(paste0(subinfodir, "n1601_go1_environment_factor_scores_tymoore_20150909.csv"))
#get the raw environment values to compare across groups
envdata <- read.csv(paste0(subinfodir, "Census_for_Scores.csv"))
#Get QA values to include in analyses
file3<-read.csv(paste0(qadir, "n1601_RestQAData_20170509.csv"))
#Get the network stats file with the different clustering coefficients
file4<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1012_net_meas_normed3clustcoefs_signed.csv")
#get the original file with modularity in it
file4<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1012_sub_net_meas_signed.csv")
#rename the second column of the network statistics file to be 'scanid' so it matches the other files
file4<-dplyr::rename(file4, scanid=subjlist_2)

#Join all files together
master<-right_join(file1, subjlist, by ="scanid")
master<-right_join(file2,master, by="scanid")
master<-right_join(file3, master, by= "scanid")
master<-right_join(file4, master, by ="scanid")
master <- dplyr::rename(master,bblid=subjlist_1)
master <-right_join(envdata, master, by="bblid")

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
#should be centering the age and envSES variables before looking at interactions?
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

#########
##### Sample Characterization ####
#########
t.test(ageAtScan1yrs~envSEShigh,data=master)
t.test(PercentInPoverty~envSEShigh, data=master)
t.test(PercentMarried~envSEShigh, data=master)
t.test(MedianFamilyIncome~envSEShigh, data=master)
chisq.test(table(master$envSEShigh, master$sex))

#Is average weight higher in females or is gender related to motion?
t.test(avgweight~sex, data=master)
t.test(restRelMeanRMSMotion~sex, data=master)
l <- lm(avgweight~restRelMeanRMSMotion+ageAtScan1yrs+race2+sex, data=master)
summary(l)

l <- lm(restRelMeanRMSMotion~ageAtScan1yrs+race2+sex, data=master)
summary(l)

##########
##### Global effects on rs-fMRI topology ####
##########
## Files from 01_z_transform_FC_matrices.m and 02_net_meas_for_subjs_signed.m
###CLUSTERING####

#non-linear relationship with age?
ageonlyRlrtmodel<-gamm(avgclustco_both~s(ageAtScan1cent)+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh, method='REML', data=master)$lme
l<-exactRLRT(ageonlyRlrtmodel)

#linear age effect without interaction
l <- lm(avgclustco_both ~ ageAtScan1yrs+sex+race2+restRelMeanRMSMotion+avgweight+envSEShigh, data=master)
summary(l)
l.beta<- lm.beta(l)
lmageplot<-visreg(l, "ageAtScan1yrs",
                  main="Average Clustering Coefficient", xlab="Age in Years (centered)", ylab="Average Clustering Coefficient (partial residuals)")

#linear model med split
l2 <- lm(avgclustco_both ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l2)
l2.beta <- lm.beta(l2)
anova(l,l2,test="Chisq")

lrtest(l,l2)
l2 <- lm(scale(avgclustco_both) ~ scale(ageAtScan1yrs)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+envSEShigh+scale(ageAtScan1yrs)*envSEShigh, data=master)
summary(l2)
#make tables
#apa.reg.table(l, filename = "Supp_Table1_APA.doc", table.number = 2)
#sjt.lm(l2, show.std = TRUE, p.zero = TRUE)
stargazer(l,l2,l.beta, l2.beta, single.row = TRUE, type="latex", report = "vcstp*", digits = NA)

##### Null model networks ######
#Get the null models file, made with 02_net_meas_for_subjs_signed_nulls.m and 04_consolidate_subjs_null_models.m

null_dir <- "~/Documents/bassett_lab/tooleyEnviNetworks/analyses/null_models_subjects/"
file5<-read.csv(paste0(null_dir, "n1012_sub_null_models.csv"))
file5<-dplyr::rename(file5, scanid=subjlist_2)
master_nulls<-right_join(file5, master, by ="scanid")

#Test the difference
t.test(master$avgclustco_both, master_nulls$avgclustco_both_null1)
t.test(master$avgclustco_both, master_nulls$avgclustco_both_null2)
l <- lm(avgclustco_both_null1~ageAtScan1yrs+sex+race2+avgweight_null1+envSEShigh+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master_nulls)
summary(l) 
l <- lm(avgclustco_both_null2~ageAtScan1yrs+sex+race2+avgweight_null2+envSEShigh+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master_nulls)
summary(l)  

#plot it
temp=data.frame(master,master_nulls)
require(reshape2)
df <- melt(data.frame(master$avgclustco_both, master_nulls$avgclustco_both_null2))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center", 
                                            position="jitter", stackratio = 0.01, binwidth =0.01,dotsize = 0.4, alpha=0.5)+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                        geom = "crossbar", width = 0.8, col=rgb(101, 58, 150, maxColorValue = 255))
######MODULARITY######

#linear age effect without interaction
l <- lm(modul ~ ageAtScan1yrs+sex+race2+avgweight+envSEShigh+restRelMeanRMSMotion, data=master)
summary(l)
lm.beta(l)
lmageplot<-visreg(l, "ageAtScan1yrs",
                  main="Average Clustering Coefficient", xlab="Age in Years (centered)", ylab="Average Clustering Coefficient (partial residuals)")

#linear model med split
l2 <- lm(modul ~ ageAtScan1yrs+sex+race2+avgweight.x+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l2)
lm.beta(l2)
anova(l,l2,test="Chisq")

lrtest(l,l2)

#modularity and clustering are highly correlated
cor.test(master$modul,master$avgclustco_both,method = "spearman")

#parsimony for modularity?
l <- lm(modul ~ ageAtScan1yrs+sex+race2+avgweight.x+avgclustco_both.x+envSEShigh+restRelMeanRMSMotion, data=master)
summary(l) #no age effect
l <- lm(modul ~ ageAtScan1yrs+sex+race2+avgweight.x+avgclustco_both.x+envSEShigh+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
summary(l) # no interaction

# #parsimony for clustco?
l <- lm(avgclustco_both.x ~ ageAtScan1yrs+sex+race2+avgweight.x+modul+envSEShigh+restRelMeanRMSMotion, data=master)
summary(l)
l <- lm(avgclustco_both.x ~ ageAtScan1yrs+sex+race2+avgweight.x+modul+envSEShigh+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
summary(l)

##############
#### FIGURE 3 : YEO SYSTEMS ###
##############
#files from 05_clustco_node_wise.m

full_nodewise_clustco<-read.csv(paste0(clustcodir,"n1012_clust_co_nodewise_by_subj.csv"))
full_nodewise_clustco<-dplyr::rename(full_nodewise_clustco, scanid=subjlist_2)
full_nodewise_clustco<-right_join(full_nodewise_clustco,master,by="scanid")
a<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/parcels/Glasser_to_Yeo.csv")
a$id<-gsub("^NZMod_", "", a$id)
#remove parcel R52 which is low signal
a<-a[-103,]
#get a dataframe with each subject with the mean clust co in each system
detach("package:ppcor", unload=TRUE)
detach("package:MASS", unload=TRUE)
full_nodewise_clustco<-select(full_nodewise_clustco, -(bblid:paredu1cent))
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

### FOR AGE BETAS ####
#Yeo 1
Yeo1 <- lm(scale(Yeo_1) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+envSEShigh+scale(restRelMeanRMSMotion), data=master)
age_beta_yeo1<-lm.beta(Yeo1)$standardized.coefficients[2]
#Yeo 2
Yeo2<- lm(scale(Yeo_2) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+envSEShigh+scale(restRelMeanRMSMotion), data=master)
age_beta_yeo2<-lm.beta(Yeo2)$standardized.coefficients[2]
#Yeo 3
Yeo3<- lm(scale(Yeo_3) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+envSEShigh+scale(restRelMeanRMSMotion), data=master)
age_beta_yeo3<-lm.beta(Yeo3)$standardized.coefficients[2]
#Yeo 4
Yeo4<- lm(scale(Yeo_4) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+envSEShigh+scale(restRelMeanRMSMotion), data=master)
age_beta_yeo4<-lm.beta(Yeo4)$standardized.coefficients[2]
#Yeo 5
Yeo5 <- lm(scale(Yeo_5) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+envSEShigh+scale(restRelMeanRMSMotion), data=master)
age_beta_yeo5<-lm.beta(Yeo5)$standardized.coefficients[2]
#Yeo 6
Yeo6 <- lm(scale(Yeo_6) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+envSEShigh+scale(restRelMeanRMSMotion), data=master)
age_beta_yeo6<-lm.beta(Yeo6)$standardized.coefficients[2]
#Yeo 7
Yeo7 <- lm(scale(Yeo_7) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+envSEShigh+scale(restRelMeanRMSMotion), data=master)
age_beta_yeo7<-lm.beta(Yeo7)$standardized.coefficients[2]

#write scaled betas to outfile
models<-list(Yeo1, Yeo2, Yeo3, Yeo4, Yeo5, Yeo6, Yeo7)
age_betas <- lapply(models, function(x) { summary(x)$coefficients[2,1]})
age_se <- lapply(models,function(x) { summary(x)$coefficients[2,2]})
age_pvals <- lapply(models,function(x) { summary(x)$coefficients[2,4]})
age_scaled_yeo_betas <- as.data.frame(cbind(unlist(age_betas),unlist(age_se), unlist(age_pvals)))
colnames(age_scaled_yeo_betas) <- c("age_betas", "age_se", "age_pvals")
age_scaled_yeo_betas$age_pvals_fdr <- p.adjust(age_scaled_yeo_betas$age_pvals, method = "fdr")
age_scaled_yeo_betas$Yeonet <- 1:7
outfile <- data.frame(age_scaled_yeo_betas, agexses_scaled_yeo_betas)
write.csv(outfile, paste0(clustcodir, "yeo_network_betas_scaled.csv"))

#### FOR AGE X SES BETAS ######
#Yeo 1
Yeo1_scaled <- lm(scale(Yeo_1) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo1_scaled)
#Yeo 2
Yeo2_scaled<- lm(scale(Yeo_2) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo2_scaled)
#Yeo 3
Yeo3_scaled<- lm(scale(Yeo_3) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo3_scaled)
#Yeo 4
Yeo4_scaled<- lm(scale(Yeo_4) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo4_scaled)
#Yeo 5
Yeo5_scaled<- lm(scale(Yeo_5) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo5_scaled)
#Yeo 6
Yeo6_scaled <- lm(scale(Yeo_6) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo6_scaled)
#Yeo 7
Yeo7_scaled <- lm(scale(Yeo_7) ~ scale(ageAtScan1cent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(Yeo7_scaled)

#scaled Yeo agexses betas
models<-list(Yeo1_scaled, Yeo2_scaled, Yeo3_scaled, Yeo4_scaled, Yeo5_scaled, Yeo6_scaled, Yeo7_scaled)
agexses_betas <- lapply(models, function(x) { summary(x)$coefficients[9,1]})
agexses_se <- lapply(models,function(x) { summary(x)$coefficients[9,2]})
agexses_pvals <- lapply(models,function(x) { summary(x)$coefficients[9,4]})
agexses_scaled_yeo_betas <- data.frame(cbind(unlist(agexses_betas),unlist(agexses_se), unlist(agexses_pvals)))
colnames(agexses_scaled_yeo_betas) <- c("agexses_betas", "agexses_se", "agexses_pvals")
agexses_scaled_yeo_betas$agexses_pvalsfdr <- p.adjust(agexses_scaled_yeo_betas$agexses_pvals,method = "fdr")
agexses_scaled_yeo_betas$Yeonet <- 1:7
outfile <- data.frame(age_scaled_yeo_betas, agexses_scaled_yeo_betas)
write.csv(outfile, paste0(clustcodir, "yeo_network_betas_scaled.csv"))

###### FIGURE 3 #######
yeo_betas <- read.csv(paste0(clustcodir, "yeo_network_betas_scaled.csv"))
#make a dummy variable to rank them
par(mfrow=c(2,2))
#Make age plot
yeo_betas$effect_age <- base::rank(-yeo_betas$age_betas)
#RGB colors of Yeo brain
age_colors <- c(rgb(196, 58, 250, maxColorValue = 255), rgb(205, 62, 78, maxColorValue = 255), rgb(0, 118, 14, maxColorValue = 255), 
                rgb(230, 148, 34, maxColorValue = 255), rgb(70, 130, 180, maxColorValue = 255), 
                rgb(120, 18, 134, maxColorValue = 255), rgb(220,248,164, maxColorValue = 255))
#all Age betas are significant, need asterisks
Fig2_Age_By_Yeo_Sys <- ggplot(data=yeo_betas, aes(effect_age, age_betas)) +geom_bar(fill=age_colors, col="black",stat="identity", size=0)
Fig2_Age_By_Yeo_Sys<-Fig2_Age_By_Yeo_Sys + scale_x_continuous(breaks=1:7, labels=c("Ventral Attention","Default","Dorsal Attention","Frontoparietal","Somatomotor","Visual","Limbic"))+ 
  theme(axis.text = element_text(size= 12), plot.title = element_text(hjust = 0.5, face = "bold")) +labs(title="Age Effects Across Cognitive Systems", x="", y="Standardized Coefficients for Age Effect")+ 
  theme(panel.border = element_blank(),axis.line.y=element_line(),panel.background = element_rect(fill = "white"))
Fig2_Age_By_Yeo_Sys+ geom_errorbar(data=yeo_betas, aes(ymin=age_betas-age_se, ymax=age_betas+age_se), width=.3, size=1.5)

#Make age by SES plot
yeo_betas$effect_age_ses <- base::rank(-yeo_betas$agexses_betas)
age_ses_colors <- c(rgb(220,248,164, maxColorValue = 255), rgb(70, 130, 180, maxColorValue = 255), rgb(196, 58, 250, maxColorValue = 255), 
                    rgb(230, 148, 34, maxColorValue = 255), rgb(205, 62, 78, maxColorValue = 255), 
                    rgb(0, 118, 14, maxColorValue = 255), rgb(120, 18, 134, maxColorValue = 255))
#all Agexses betas are significant except visual (last), doesn't need asterisks
Fig2_Age_SES_By_Yeo_Sys <- ggplot(data=yeo_betas, aes(effect_age_ses, agexses_betas)) + geom_bar(size=0, fill=age_ses_colors, col="black",stat="identity",position=position_dodge())
Fig2_Age_SES_By_Yeo_Sys<-Fig2_Age_SES_By_Yeo_Sys + scale_x_continuous(breaks=1:7, labels=c("Limbic","Somatomotor","Ventral Attention","Frontoparietal","Default","Dorsal Attention","Visual"))+ theme(axis.text = element_text(size= 12), plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Age x SES Effects Across Cognitive Systems", x="", y="Standardized Coefficients for Age x SES Effect")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white"), panel.border= element_blank())
Fig2_Age_SES_By_Yeo_Sys+ geom_errorbar(data=yeo_betas, aes(ymin=agexses_betas-agexses_se, ymax=agexses_betas+agexses_se), width=.3, size=1.5)

#############
##### FIGURE 4: REGIONAL EFFECTS #######
##################
#from file 05_clustco_node_wise.m 

#import data on clustering coef for all nodes for all subjects
clustcodir="~/Dropbox (Personal)/bassett_lab/clustco_paper/"
full_nodewise_clustco<-read.csv(paste0(clustcodir,"n1012_clust_co_nodewise_by_subj.csv"))

full_nodewise_clustco<-dplyr::rename(full_nodewise_clustco, scanid=subjlist_2)
#Join all needed datasets together
full_nodewise_clustco<-left_join(full_nodewise_clustco, file1, by ="scanid")
full_nodewise_clustco<-left_join(full_nodewise_clustco, file2, by ="scanid")
full_nodewise_clustco<-left_join(full_nodewise_clustco, file3, by ="scanid")
full_nodewise_clustco<-left_join(full_nodewise_clustco, file4, by ="scanid")

#eliminate extraneous columns
full_nodewise_clustco<-full_nodewise_clustco %>% dplyr::select(., -c(restExcludeVoxelwise, restNoDataExclude))
full_nodewise_clustco<-full_nodewise_clustco %>% dplyr::select(., -c(restRelMeanRMSMotionExclude:restRpsMapCorrectionNotApplied))
full_nodewise_clustco<-full_nodewise_clustco %>% dplyr::select(., -c(bblid.y, bblid.x))

#make sure factor variables are factors
full_nodewise_clustco$race<-factor(full_nodewise_clustco$race)
full_nodewise_clustco$race2<-factor(full_nodewise_clustco$race2, labels=c("White", "Black", "Other"))
full_nodewise_clustco$sex<-factor(full_nodewise_clustco$sex, labels=c("Male", "Female"))

#make an age squared variable
full_nodewise_clustco$ageatscansqdem <- (full_nodewise_clustco$ageAtScan1-mean(full_nodewise_clustco$ageAtScan1))^2
full_nodewise_clustco$ageatscansq <- (full_nodewise_clustco$ageAtScan1)^2
full_nodewise_clustco$ageAtScan1yrs<-(full_nodewise_clustco$ageAtScan1)/12
#split SES on the median 
full_nodewise_clustco$envSEShigh=NA
full_nodewise_clustco$envSEShigh[full_nodewise_clustco$envSES >= 0.0178] <- 1
full_nodewise_clustco$envSEShigh[full_nodewise_clustco$envSES < 0.0178]<- 0
full_nodewise_clustco$envSEShigh<-factor(full_nodewise_clustco$envSEShigh, labels=c("Low", "High"), ordered=TRUE)

#should potentially be centering the age and envSES variables before looking at interactions
full_nodewise_clustco$ageAtScan1cent<-(full_nodewise_clustco$ageAtScan1-mean(full_nodewise_clustco$ageAtScan1))
full_nodewise_clustco$ageAtScan1yrscent<-(full_nodewise_clustco$ageAtScan1yrs-mean(full_nodewise_clustco$ageAtScan1yrs))
full_nodewise_clustco$envSEScent<-(full_nodewise_clustco$envSES-mean(full_nodewise_clustco$envSES))

#LOOK AT EACH NODE FOR AGE EFFECT IN THE LINEAR MODEL
covariates=" ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion"
m <- mclapply(names(full_nodewise_clustco[,3:361]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
NodeWise_Clustco_lm_Age_pvals <- mclapply(m, function(x) { summary(lm(formula = x,data=full_nodewise_clustco))$coef[2,4]},mc.cores=1)
NodeWise_Clustco_lm_Age_pvals <- as.data.frame(NodeWise_Clustco_lm_Age_pvals)
NodeWise_Clustco_lm_Age_pvals <- t(NodeWise_Clustco_lm_Age_pvals)
NodeWise_Clustco_lm_Age_pvals <- as.data.frame(NodeWise_Clustco_lm_Age_pvals)
colnames(NodeWise_Clustco_lm_Age_pvals) <- "NodeWise_Clustco_lm_Age_pvals"
NodeWise_Clustco_lm_Age_pvals$Node_index <- 1:359
#Bonferroni
Bonferronip=0.05/dim(NodeWise_Clustco_lm_Age_pvals)[1]
sig_nodes_lm_age_bonf<-NodeWise_Clustco_lm_Age_pvals[NodeWise_Clustco_lm_Age_pvals$NodeWise_Clustco_lm_Age_pvals < Bonferronip,]
#FDR
fdr_corrected<-p.adjust(NodeWise_Clustco_lm_Age_pvals$NodeWise_Clustco_lm_Age_pvals, method = "fdr")
sig_nodes_lm_age<-cbind(NodeWise_Clustco_lm_Age_pvals, fdr_corrected)
fdr_sig_nodes_lm_age<-sig_nodes_lm_age[sig_nodes_lm_age$fdr_corrected < 0.05,]
#get betas for effect estimates
NodeWise_Clustco_lm_Age_betas <- mclapply(m, function(x) { lm.beta(lm(formula = x,data=full_nodewise_clustco))$standardized.coefficients[[2]]},mc.cores=1)
NodeWise_Clustco_lm_Age_betas <- as.data.frame(NodeWise_Clustco_lm_Age_betas)
NodeWise_Clustco_lm_Age_betas <- t(NodeWise_Clustco_lm_Age_betas)
NodeWise_Clustco_lm_Age_betas <- as.data.frame(NodeWise_Clustco_lm_Age_betas)
colnames(NodeWise_Clustco_lm_Age_betas) <- "NodeWise_Clustco_lm_Age_betas"
NodeWise_Clustco_lm_Age_betas$Node_index <- 1:359

#LOOK AT EACH NODE FOR AGE*SES EFFECT IN THE LINEAR MODEL
covariates=" ~ ageAtScan1cent+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh"
m <- mclapply(names(full_nodewise_clustco[,3:361]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
NodeWise_Clustco_lm_AgexSES_pvals <- mclapply(m, function(x) { summary(lm(formula = x,data=full_nodewise_clustco))$coef[9,4]},mc.cores=1)
NodeWise_Clustco_lm_AgexSES_pvals <- as.data.frame(NodeWise_Clustco_lm_AgexSES_pvals)
NodeWise_Clustco_lm_AgexSES_pvals <- t(NodeWise_Clustco_lm_AgexSES_pvals)
NodeWise_Clustco_lm_AgexSES_pvals <- as.data.frame(NodeWise_Clustco_lm_AgexSES_pvals)
colnames(NodeWise_Clustco_lm_AgexSES_pvals) <- "NodeWise_Clustco_lm_AgexSES_pvals"
NodeWise_Clustco_lm_AgexSES_pvals$Node_index <- 1:359
#Bonferroni
Bonferronip=0.05/dim(NodeWise_Clustco_lm_AgexSES_pvals)[1]
sig_nodes_lm_AgexSES_bonf<-NodeWise_Clustco_lm_AgexSES_pvals[NodeWise_Clustco_lm_AgexSES_pvals$NodeWise_Clustco_lm_AgexSES_pvals < Bonferronip,]
#FDR
fdr_corrected<-p.adjust(NodeWise_Clustco_lm_AgexSES_pvals$NodeWise_Clustco_lm_AgexSES_pvals, method = "fdr")
sig_nodes_lm_AgexSES<-cbind(NodeWise_Clustco_lm_AgexSES_pvals, fdr_corrected)
fdr_sig_nodes_lm_AgexSES<-sig_nodes_lm_AgexSES[sig_nodes_lm_AgexSES$fdr_corrected < 0.05,]

#get betas for effect estimates
NodeWise_Clustco_lm_AgexSES_betas <- mclapply(m, function(x) { lm.beta(lm(formula = x,data=full_nodewise_clustco))$standardized.coefficients[[9]]},mc.cores=1)
NodeWise_Clustco_lm_AgexSES_betas <- as.data.frame(NodeWise_Clustco_lm_AgexSES_betas)
NodeWise_Clustco_lm_AgexSES_betas <- t(NodeWise_Clustco_lm_AgexSES_betas)
NodeWise_Clustco_lm_AgexSES_betas <- as.data.frame(NodeWise_Clustco_lm_AgexSES_betas)
colnames(NodeWise_Clustco_lm_AgexSES_betas) <- "NodeWise_Clustco_lm_AgexSES_betas"
NodeWise_Clustco_lm_AgexSES_betas$Node_index <- 1:359

#PULL OUT NODES WITH HIGHEST INTERACTION EFFECT, PLOT THESE ONLY
dim(fdr_sig_nodes_lm_AgexSES)
#make node names to pull out from the full matrix
nodes<-fdr_sig_nodes_lm_AgexSES$Node_index
nodes[1] <- "..8"
nodes[2:9]=c(".36" , ".38" , ".40" , ".41"  ,".53" , ".57" , ".60",  ".66")
node_names<-paste("avgclustco_both_",nodes, sep="")
#average across all 26 nodes to create one "significant area" per subject
full_nodewise_clustco <- full_nodewise_clustco %>% mutate(., mean_peaks=rowMeans(select(.,one_of(node_names))))
#analyze and plot
peak_areas_only<- lm(mean_peaks ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=full_nodewise_clustco)
peak_areas_only2<- lm(scale(mean_peaks) ~ scale(ageAtScan1yrscent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1yrscent)*envSEShigh, data=full_nodewise_clustco)
summary(peak_areas_only)
#regular visreg plot
visreg(peak_areas_only, "ageAtScan1yrs", by ="envSEShigh", main="Mean Clustering Coefficient of Peak Regions",
       xlab="Age in Years", overlay=TRUE, ylab=" Mean Clustering Coefficient (partial residuals)", strip.names=c("Low SES", "High SES"))

#WRITE OUT A CSV OF BETAS FOR MATLAB PLOTTING ON BRAIN
outfile<-cbind(NodeWise_Clustco_lm_Age_betas$NodeWise_Clustco_lm_Age_betas, NodeWise_Clustco_lm_AgexSES_betas$NodeWise_Clustco_lm_AgexSES_betas)
colnames(outfile)<-c("lm_age_betas", "lm_agexses_betas")
write.csv(outfile, paste0("~/Dropbox (Personal)/bassett_lab/clustco_paper/nodewise_betas_for_lms.csv"))

#WRITE OUT A CSV OF NODE PVALS TO BE READ INTO MATLAB FOR PLOTTING ON BRAIN
outfile<-cbind(sig_nodes_GAM_age$fdr_corrected, sig_nodes_GAM_sesint_medsp$fdr_corrected, 
               sig_nodes_lm_age$fdr_corrected, sig_nodes_lm_AgexSES$fdr_corrected)
colnames(outfile)<-c("fdr_pvals_GAM_age", "fdr_pvals_GAM_agexses", "fdr_pvals_lm_age", "fdr_pvals_lm_agexses")
write.csv(outfile, paste0("~/Dropbox (Personal)/bassett_lab/clustco_paper/nodewise_pvals_for_GAMs_and_lms.csv"))

#############
##### FIGURE 4: REHO #######
##################

#Get the reho file by subject
file4<-read.csv(paste0(reho_dir,"n1601_glasserReHoValues_20170509.csv"))
master<-right_join(file4, master, by ="scanid")
#remove parcel R52 which is missing signal
#unload ppcor and MASS
detach("package:ppcor", unload=TRUE)
detach("package:MASS", unload=TRUE)
master<-select(master, -rest_glasser_reho_Right_52)
#average whole-brain reho for each subject
master <- master %>% rowwise() %>% mutate(mean_reho=mean(rest_glasser_reho_Right_V1:rest_glasser_reho_Left_p24))
summary(master$mean_reho)

#are reho and clustco correlated?
cor.test(master$mean_reho, master$avgclustco_both)
#what about when controlling for other covariates?
library(ppcor) #this will mask dplyr
temp<-cbind(master$ageAtScan1cent,as.factor(master$sex), as.factor(master$race2), master$restRelMeanRMSMotion)
pcor.test(master$avgclustco_both, master$mean_reho, temp)

### WHOLE BRAIN REHO MODEL ###
l <- lm(mean_reho~ageAtScan1cent+race2+sex+restRelMeanRMSMotion+envSEShigh, data=master)
summary(l)
lm.beta(l)

#controlling for whole-brain reho when examining the clustering coefficient
l <- lm(avgclustco_both~ageAtScan1cent+race2+sex+mean_reho+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh, data=master)
summary(l)
lm.beta(l)

#remove packages as they mask dplyr
detach("package:ppcor", unload=TRUE)
detach("package:MASS", unload=TRUE)
#get the mean of each node's clustco across subjects
mean_nodes_clustco<-colMeans(select(full_nodewise_clustco, avgclustco_both_..1:avgclustco_both_359))
#get the mean of each node's reho across subjects
mean_nodes_reho<-select(master, rest_glasser_reho_Right_V1:rest_glasser_reho_Left_p24) %>% colMeans(.)
#write this to a csv
outfile<-mean_nodes_reho
write.csv(outfile, paste0("~/Dropbox (Personal)/bassett_lab/clustco_paper/avg_reho_across_nodes.csv"))
#see if they're correlated
cor.test(mean_nodes_reho, mean_nodes_clustco) #they are across regions!
#what about when controlling for other things, within each subject
#create two dataframes of clustco and reho only
nodewise_clustco <- select(full_nodewise_clustco, avgclustco_both_..1:avgclustco_both_359)
nodewise_reho <- select(master, rest_glasser_reho_Right_V1:rest_glasser_reho_Left_p24)
correlationspvals=numeric(359)
pvalsforclustco=numeric(359)
library(ppcor)
for (i in 1:359){
  clustco=nodewise_clustco[i]
  reho=nodewise_reho[i]
  #temp=cbind(master$ageAtScan1cent,as.factor(master$sex), as.factor(master$race2), master$restRelMeanRMSMotion,as.factor(master$envSEShigh), master$avgweight, clustco,reho)
  #colnames(temp)=c("ageAtScan1cent", "sex", "race2", "restRelMeanRMSMotion", "envSEShigh", "avgweight","clustco", "reho")
  #pvalsforclustco[[i]]<- summary(lm(clustco~ageAtScan1cent+avgweight+race2+sex+restRelMeanRMSMotion+envSEShigh+reho+ageAtScan1cent*envSEShigh, data=temp))$coef[9,4]
  temp=cbind(master$ageAtScan1cent,as.factor(master$sex), as.factor(master$race2), master$restRelMeanRMSMotion,master$envSEShigh)
  correlationspvals[[i]] <- as.numeric(pcor.test(clustco, reho, temp)$p.value)
}

#see if there are any nodes for which clustco and reho are not correlated significantly
ps <- p.adjust(correlationspvals, method="fdr")
ps[ps>0.05] #nope!
correlationspvals[correlationspvals<0.05]
correlationspvals <- as.data.frame(correlationspvals)
correlationspvals <- t(correlationspvals)
correlationspvals <- as.data.frame(correlationspvals)

### REGIONAL REHO MODELS ###

#LOOK AT EACH NODE FOR AGE EFFECT IN THE LINEAR MODEL
covariates=" ~ ageAtScan1cent+sex+race2+restRelMeanRMSMotion+envSEShigh"
#check that the names with this dataframe pulls rest_glasser variables
m <- mclapply(names(master[,3:361]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
NodeWise_reho_lm_Age_pvals <- mclapply(m, function(x) { summary(lm(formula = x,data=master))$coef[2,4]},mc.cores=1)
NodeWise_Reho_lm_Age_pvals <- as.data.frame(NodeWise_reho_lm_Age_pvals)
NodeWise_Reho_lm_Age_pvals <- t(NodeWise_Reho_lm_Age_pvals)
NodeWise_Reho_lm_Age_pvals <- as.data.frame(NodeWise_Reho_lm_Age_pvals)
colnames(NodeWise_Reho_lm_Age_pvals) <- "NodeWise_Reho_lm_Age_pvals"
NodeWise_Reho_lm_Age_pvals$Node_index <- 1:359
#Bonferroni
Bonferronip=0.05/dim(NodeWise_Reho_lm_Age_pvals)[1]
sig_nodes_lm_age_bonf<-NodeWise_Reho_lm_Age_pvals[NodeWise_Reho_lm_Age_pvals$NodeWise_Reho_lm_Age_pvals < Bonferronip,]
#FDR
fdr_corrected<-p.adjust(NodeWise_Reho_lm_Age_pvals$NodeWise_Reho_lm_Age_pvals, method = "fdr")
NodeWise_Reho_lm_Age_pvals$fdr_corrected <- fdr_corrected
sig_nodes_lm_age<-cbind(NodeWise_Reho_lm_Age_pvals, fdr_corrected)
fdr_sig_nodes_lm_age<-sig_nodes_lm_age[sig_nodes_lm_age$fdr_corrected < 0.05,]

#get betas for direction of the effect estimates
NodeWise_Reho_lm_Age_betas <- mclapply(m, function(x) { lm.beta(lm(formula = x,data=master))$standardized.coefficients[[2]]},mc.cores=1)
NodeWise_Reho_lm_Age_betas <- as.data.frame(NodeWise_Reho_lm_Age_betas)
NodeWise_Reho_lm_Age_betas <- t(NodeWise_Reho_lm_Age_betas)
NodeWise_Reho_lm_Age_betas <- as.data.frame(NodeWise_Reho_lm_Age_betas)
colnames(NodeWise_Reho_lm_Age_betas) <- "NodeWise_Reho_lm_Age_betas"
NodeWise_Reho_lm_Age_betas$Node_index <- 1:359

#where is ReHo not decreasing with age, these 3 regions
increases<-NodeWise_Reho_lm_Age_betas[NodeWise_Reho_lm_Age_betas$NodeWise_Reho_lm_Age_betas>0,]

#LOOK AT EACH NODE FOR AGE*SES EFFECT IN THE LINEAR MODEL
covariates=" ~ ageAtScan1cent+sex+race2+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh"
m <- mclapply(names(master[,3:361]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
NodeWise_Reho_lm_AgexSES_pvals <- mclapply(m, function(x) { summary(lm(formula = x,data=master))$coef[8,4]},mc.cores=1)
NodeWise_Reho_lm_AgexSES_pvals <- as.data.frame(NodeWise_Reho_lm_AgexSES_pvals)
NodeWise_Reho_lm_AgexSES_pvals <- t(NodeWise_Reho_lm_AgexSES_pvals)
NodeWise_Reho_lm_AgexSES_pvals <- as.data.frame(NodeWise_Reho_lm_AgexSES_pvals)
colnames(NodeWise_Reho_lm_AgexSES_pvals) <- "NodeWise_Reho_lm_AgexSES_pvals"
NodeWise_Reho_lm_AgexSES_pvals$Node_index <- 1:359
detach("package:ppcor", unload=TRUE)
detach("package:MASS", unload=TRUE)
names <- colnames(select(master, rest_glasser_reho_Right_V1:rest_glasser_reho_Left_p24))
NodeWise_Reho_lm_AgexSES_pvals$Names <- names
#Bonferroni
Bonferronip=0.05/dim(NodeWise_Reho_lm_AgexSES_pvals)[1]
sig_nodes_lm_AgexSES_bonf<-NodeWise_Reho_lm_AgexSES_pvals[NodeWise_Reho_lm_AgexSES_pvals$NodeWise_Reho_lm_AgexSES_pvals < Bonferronip,]
#FDR
fdr_corrected<-p.adjust(NodeWise_Reho_lm_AgexSES_pvals$NodeWise_Reho_lm_AgexSES_pvals, method = "fdr")
NodeWise_Reho_lm_AgexSES_pvals$fdr_corrected <- fdr_corrected
sig_nodes_reho_lm_AgexSES<-cbind(NodeWise_Reho_lm_AgexSES_pvals, fdr_corrected)
fdr_sig_nodes_reho_lm_AgexSES<-sig_nodes_reho_lm_AgexSES[sig_nodes_reho_lm_AgexSES$fdr_corrected < 0.05,]

#get betas for effect estimates
NodeWise_Reho_lm_AgexSES_betas <- mclapply(m, function(x) { lm.beta(lm(formula = x,data=master))$standardized.coefficients[[8]]},mc.cores=1)
NodeWise_Reho_lm_AgexSES_betas <- as.data.frame(NodeWise_Reho_lm_AgexSES_betas)
NodeWise_Reho_lm_AgexSES_betas <- t(NodeWise_Reho_lm_AgexSES_betas)
NodeWise_Reho_lm_AgexSES_betas <- as.data.frame(NodeWise_Reho_lm_AgexSES_betas)
colnames(NodeWise_Reho_lm_AgexSES_betas) <- "NodeWise_Reho_lm_AgexSES_betas"
NodeWise_Reho_lm_AgexSES_betas$Node_index <- 1:359

#write out a csv of nodewise betas and pvals
outfile<-cbind(NodeWise_Reho_lm_Age_betas$NodeWise_Reho_lm_Age_betas, NodeWise_Reho_lm_Age_pvals$fdr_corrected, NodeWise_Reho_lm_AgexSES_betas$NodeWise_Reho_lm_AgexSES_betas, NodeWise_Reho_lm_AgexSES_pvals$fdr_corrected)
colnames(outfile)<-c("reho_age_effect_betas", "reho_age_effect_pvals_fdr", "reho_agexses_effect_betas", "reho_agexses_effect_pvals_fdr")
write.csv(outfile, paste0("~/Dropbox (Personal)/bassett_lab/clustco_paper/nodewise_estimates_for_reho_effect_corrected.csv"))

### FIG 5D ####

#PULL OUT NODES WITH HIGHEST INTERACTION EFFECT, PLOT THESE ONLY
dim(fdr_sig_nodes_reho_lm_AgexSES)
#make node names to pull out from the full matrix
nodes<-fdr_sig_nodes_lm_AgexSES$Node_index

nodes[1] <- "..8"
nodes[2:21]=c(".32", ".36", ".37", ".38", ".39", ".40", ".41", ".43", ".53", ".55", ".56", ".57", ".58", ".59", ".60", ".63", ".65", ".76", ".88", ".94")
node_names<-paste("avgclustco_both_",nodes, sep="")
#average across all nodes to create one "significant area" per subject
master <- master %>% ungroup(.) %>%  mutate(., mean_peaks_reho=rowMeans(select(.,one_of(fdr_sig_nodes_reho_lm_AgexSES$Names))))
#analyze and plot
peak_areas_only<- lm(mean_peaks_reho ~ ageAtScan1yrs+sex+race2+restRelMeanRMSMotion+ageAtScan1yrs*envSEShigh, data=master)
peak_areas_only2<- lm(scale(mean_peaks_reho) ~ scale(ageAtScan1yrscent)+sex+race2+scale(avgweight)+scale(restRelMeanRMSMotion)+scale(ageAtScan1yrscent)*envSEShigh, data=master)

#plot it
visreg(peak_areas_only, "ageAtScan1yrs", by ="envSEShigh", main="Mean Clustering Coefficient of Peak Regions",
       xlab="Age in Years", ylab=" Mean Clustering Coefficient (partial residuals)", overlay=TRUE, partial=FALSE, rug=FALSE, 
       line=list(col=c(rgb(28, 147, 255, maxColorValue = 255), rgb(255, 168, 28, maxColorValue = 255))), 
       fill=list(col=c(alpha(rgb(28, 147, 255, maxColorValue = 255), 0.7), alpha(rgb(255, 168, 28, maxColorValue = 255),0.7))))

### REGRESS OUT REHO FROM CLUSTERING, LOOK AT THE EFFECT ACROSS REGIONS
pvalsforclustco=numeric(359)
for (i in 1:359){
  clustco=nodewise_clustco[i]
  reho=nodewise_reho[i]
  temp=cbind(master$ageAtScan1cent,as.factor(master$sex), as.factor(master$race2), master$restRelMeanRMSMotion,as.factor(master$envSEShigh), master$avgweight, clustco,reho)
  colnames(temp)=c("ageAtScan1cent", "sex", "race2", "restRelMeanRMSMotion", "envSEShigh", "avgweight","clustco", "reho")
  temp$residclustco <- resid(lm(clustco~reho, data=temp))
  #temp$residclustco <- mean(temp$clustco)+temp$residclustco
  pvalsforclustco[[i]]<- summary(lm(residclustco~ageAtScan1cent+avgweight+race2+sex+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh, data=temp))$coef[9,4]
}
#see if there are any nodes for agexses interaction remains significant when regressed out reho
sig_nodes <- pvalsforclustco[pvalsforclustco<0.05]
length(sig_nodes)
ps <- p.adjust(pvalsforclustco, method="fdr")
ps[ps>0.05] #nope!
fdr_sig_nodes <- ps[ps<0.05]
length(fdr_sig_nodes)

#WHOLE BRAIN MODEL WITH REHO DOES NOT FIT BETTER, AFFECT SIGNIFICANCE OF CLUST INTERACTION
l <- lm(scale(avgclustco_both)~scale(ageAtScan1cent)+race2+sex+scale(avgweight)+scale(restRelMeanRMSMotion)+envSEShigh+scale(ageAtScan1cent)*envSEShigh, data=master)
l2 <- lm(scale(avgclustco_both)~scale(ageAtScan1cent)+race2+sex+scale(avgweight)+scale(mean_reho)+scale(restRelMeanRMSMotion)+envSEShigh+scale(ageAtScan1cent)*envSEShigh, data=master)
summary(l2)
lrtest(l,l2)
anova(l,l2, test="Chisq")

###############
### FIGURE 5: DISTANCE DEPENDENCE ###
##############
#From files 08_net_meas_for_subjs_distbins_dependence.m and 09_net_meas_for_subjs_signed_distbins_nulls.m and 09_consolidate_subjs_nulls_distbins_dependence.m

##DISTANCE
file7<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1012_sub_net_meas_signed_distance_bins.csv")
file7<-dplyr::rename(file7, scanid=subjlist_2)
master<-right_join(file7, master, by ="scanid")
#look at each bin for distance-dependence effect of age x SES in linear model
covariates=" ~ ageAtScan1cent+sex+race2+restRelMeanRMSMotion+envSEShigh+ageAtScan1cent*envSEShigh"
m <- mclapply(names(master[,23:42]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
distance_bins_agexses_betas <- mclapply(m, function(x) { lm.beta(lm(formula = x,data=master))$standardized.coefficient[8]},mc.cores=1)
distance_bins_agexses_betas <- as.data.frame(distance_bins_agexses_betas)
distance_bins_agexses_betas <- t(distance_bins_agexses_betas)
distance_bins_agexses_betas <- as.data.frame(distance_bins_agexses_betas)
x<-cbind(distance_bins_agexses_betas, as.character(m))
colnames(x) <- c("beta_for_agexses_interaction", "formula")
write.csv(x, file= "~/Dropbox/bassett_lab/analyses/csv/distance_weight_binned_agexses_betas.csv")

x %>% filter(.,formula == starts_with("avgweight"))
avgweight_only<-x[grep("avgweight*", x$formula),]

# INCLUDE AVERAGE WEIGHT IN EACH MODEL for distance dependence 
t0to1model_clustco_distances<-lm(scale(avgclustco_both_0to1_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_0to1_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t1to2model_clustco_distances<-lm(scale(avgclustco_both_1to2_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_1to2_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t2to3model_clustco_distances<-lm(scale(avgclustco_both_2to3_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_2to3_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t3to4model_clustco_distances<-lm(scale(avgclustco_both_3to4_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_3to4_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t4to5model_clustco_distances<-lm(scale(avgclustco_both_4to5_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_4to5_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t5to6model_clustco_distances<-lm(scale(avgclustco_both_5to6_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_5to6_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t6to7model_clustco_distances<-lm(scale(avgclustco_both_6to7_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_6to7_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t7to8model_clustco_distances<-lm(scale(avgclustco_both_7to8_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_7to8_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t8to9model_clustco_distances<-lm(scale(avgclustco_both_8to9_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_8to9_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t9to10model_clustco_distances<-lm(scale(avgclustco_both_9to10_longest) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_9to10_longest)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
#make a list of models
mna<-list(t0to1model_clustco_distances, t1to2model_clustco_distances, t2to3model_clustco_distances, t3to4model_clustco_distances, t4to5model_clustco_distances, t5to6model_clustco_distances, t6to7model_clustco_distances, t7to8model_clustco_distances, t8to9model_clustco_distances, t9to10model_clustco_distances)
#get the betas for the age x SES interaction
distance_bins_agexses_betas <- lapply(mna, function(x) { summary(x)$coefficients[9,1]})
distance_bins_agexses_se <- lapply(mna, function(x) { summary(x)$coefficients[9,2]})
distance_bins_agexses_betas_table <-data.frame(cbind(distance_bins_agexses_betas, distance_bins_agexses_se))

#reverse the list so the shortest connections are first
revers_distance_bins_agexses_betas <- rev(distance_bins_agexses_betas_table$distance_bins_agexses_betas)
revers_distance_bins_agexses_se <- rev(distance_bins_agexses_betas_table$distance_bins_agexses_se)
revers_distance_bins_agexses <- data.frame(cbind(revers_distance_bins_agexses_betas, revers_distance_bins_agexses_se))

#### GET NULL MODEL CLUST CO 
nullmodels <- read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/null_models_subjects_dist_bins/n1012_sub_null_models_dist_bins.csv")
nullmodels<-dplyr::rename(nullmodels, scanid=subjlist_2)
master<-right_join(nullmodels, master, by ="scanid")
## INCLUDE AVERAGE WEIGHT IN EACH MODEL for distance dependence 
t0to1model_clustco_distances_null<-lm(scale(avgclustco_both_0to1_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_0to1_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t1to2model_clustco_distances_null<-lm(scale(avgclustco_both_1to2_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_1to2_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t2to3model_clustco_distances_null<-lm(scale(avgclustco_both_2to3_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_2to3_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t3to4model_clustco_distances_null<-lm(scale(avgclustco_both_3to4_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_3to4_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t4to5model_clustco_distances_null<-lm(scale(avgclustco_both_4to5_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_4to5_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t5to6model_clustco_distances_null<-lm(scale(avgclustco_both_5to6_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_5to6_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t6to7model_clustco_distances_null<-lm(scale(avgclustco_both_6to7_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_6to7_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t7to8model_clustco_distances_null<-lm(scale(avgclustco_both_7to8_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_7to8_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t8to9model_clustco_distances_null<-lm(scale(avgclustco_both_8to9_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_8to9_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
t9to10model_clustco_distances_null<-lm(scale(avgclustco_both_9to10_longestnull) ~ scale(ageAtScan1cent) + sex + race2 + scale(avgweight_9to10_longestnull)+ scale(restRelMeanRMSMotion) + envSEShigh + scale(ageAtScan1cent) * envSEShigh, data=master)
#make a list of models
nullmodels<-list(t0to1model_clustco_distances_null, t1to2model_clustco_distances_null, t2to3model_clustco_distances_null, t3to4model_clustco_distances_null, t4to5model_clustco_distances_null, t5to6model_clustco_distances_null, t6to7model_clustco_distances_null, t7to8model_clustco_distances_null, t8to9model_clustco_distances_null, t9to10model_clustco_distances_null)
#get the betas for the age x SES interaction
distance_bins_agexses_null_betas <- lapply(nullmodels, function(x) { summary(x)$coefficients[9,1]})
distance_bins_agexses_null_se <- lapply(nullmodels, function(x) { summary(x)$coefficients[9,2]})
distance_bins_agexses_betas_null <-data.frame(cbind(distance_bins_agexses_null_betas, distance_bins_agexses_null_se))
#reverse the list so the shortest connections are first
revers_distance_bins_agexses_null_betas <- rev(distance_bins_agexses_betas_null$distance_bins_agexses_null_betas)
revers_distance_bins_agexses_null_se <- rev(distance_bins_agexses_betas_null$distance_bins_agexses_null_se)
revers_distance_bins_agexses_null <- cbind(revers_distance_bins_agexses_null_betas, revers_distance_bins_agexses_null_se)

#compare models for each distance bin to the null model for that distance bin
lrtest(t0to1model_clustco_distances_null,t0to1model_clustco_distances)
lrtest(t1to2model_clustco_distances_null,t1to2model_clustco_distances)
lrtest(t2to3model_clustco_distances_null,t2to3model_clustco_distances)
lrtest(t3to4model_clustco_distances_null,t3to4model_clustco_distances)
lrtest(t4to5model_clustco_distances_null,t4to5model_clustco_distances)
lrtest(t5to6model_clustco_distances_null,t5to6model_clustco_distances)
lrtest(t6to7model_clustco_distances_null,t6to7model_clustco_distances)
lrtest(t7to8model_clustco_distances_null,t7to8model_clustco_distances)
lrtest(t8to9model_clustco_distances_null,t8to9model_clustco_distances)
lrtest(t9to10model_clustco_distances_null,t9to10model_clustco_distances) #all empirical data does better than null model data

#plot them
myplot <- plot(1:10,revers_distance_bins_agexses$revers_distance_bins_agexses_betas,type="p",main="Distance-Dependence of Age x SES Clustering Effect",col="blue", xlab="Distance (closest to farthest)", 
               bg="blue",ylim = c(-0.10, 0.20),ylab="Standardized Age x SES Effect on Clustering Coefficient", pch=23)
arrows(1:10, as.numeric(revers_distance_bins_agexses_betas)-as.numeric(revers_distance_bins_agexses_se), 1:10, as.numeric(revers_distance_bins_agexses_betas)+as.numeric(revers_distance_bins_agexses_se), length=0.05, angle=90, code=3)
#add null betas and SEs
par(new=TRUE)
plot(1:10,revers_distance_bins_agexses_null_betas,main="",col="black", xlab="", 
     bg="blue",ylim = c(-0.10, 0.20),ylab="", pch=23)
arrows(1:10, as.numeric(revers_distance_bins_agexses_null_betas)-as.numeric(revers_distance_bins_agexses_null_se), 1:10, as.numeric(revers_distance_bins_agexses_null_betas)+as.numeric(revers_distance_bins_agexses_null_se), length=0.05, angle=90, code=3)

### Betas of edges within and between sig nodes
# see file 10_make_all_subs_all_edges.m and 11_averagebetas_within_26_nodes_outside.m

###############
### FIGURE 6: SENSITIVITY ANALYSIS ###
##############
### IMPORT DATA
# Get demographics to control for
file1<-read.csv(paste0(subinfodir, "n1601_demographics_go1_20161212.csv"))
file2<-read.csv(paste0(subinfodir, "n1601_go1_environment_factor_scores_tymoore_20150909.csv"))
#get the raw environment values to compare across groups
envdata <- read.csv(paste0(subinfodir, "Census_for_Scores.csv"))
#Get QA values to include in analyses
file3<-read.csv(paste0(qadir, "n1601_RestQAData_20170509.csv"))
#get the original file with modularity in it
file4<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/analyses/n1012_sub_net_meas_signed.csv")
#rename the second column of the network statistics file to be 'scanid' so it matches the other files
file4<-dplyr::rename(file4, scanid=subjlist_2)

#use LTN exclude critera
subjlist<-read.csv(paste0(sublistdir,"n885_LTNexclude.csv"))
#Join all files together
master<-right_join(file1, subjlist, by ="scanid")
master<-right_join(file2,master, by="scanid")
master<-right_join(file3, master, by= "scanid")
master<-right_join(file4, master, by ="scanid")

#### DATA CLEANING

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
#split SES on the median for this sample
summary(master$envSES)
master$envSEShigh=NA
master$envSEShigh[master$envSES >= -0.06135] <- 1
master$envSEShigh[master$envSES < -0.06135]<- 0
master$envSEShigh<-factor(master$envSEShigh, labels=c("Low", "High"), ordered=TRUE)
#center
master$envSEScent<-(master$envSES-mean(master$envSES))
master$medu1cent=(master$medu1-mean(master$medu1[!is.na(master$medu1)]))

#linear age effect without interaction
l <- lm(avgclustco_both ~ ageAtScan1yrs+sex+race2+avgweight+envSEShigh+restRelMeanRMSMotion, data=master)
summary(l)
lm.beta(l)
lmageplot<-visreg(l, "ageAtScan1yrs",
                  main="Average Clustering Coefficient", xlab="Age in Years (centered)", ylab="Average Clustering Coefficient (partial residuals)")

#linear model med split
l2 <- lm(avgclustco_both ~ ageAtScan1yrs+sex+race2+avgweight+restRelMeanRMSMotion+envSEShigh+ageAtScan1yrs*envSEShigh, data=master)
summary(l2)
lm.beta(l2)
anova(l,l2,test="Chisq")

lrtest(l,l2)

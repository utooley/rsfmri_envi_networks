setwd("/data/jag/bassett-lab/tooleyEnviNetworks/data/rest/")
setwd("~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/")
library(dplyr)
library(psych)

subjlist<-read.csv("/data/joy/BBL/projects/tooleyEnviNetworks/subjectLists/n1015_healthT1RestExclude.csv")
subjlist<-read.csv("~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists/n1015_healthT1RestExclude.csv")
head(subjlist)

test<-read.table("restNetwork_GlasserPNC/GlasserPNCTimeseries/4022_GlasserPNC_ts.1D")
test2<-read.table("restNetwork_GlasserPNC/GlasserPNCNetworks/4022_GlasserPNC_network.txt")
head(test2)
dim(test2)

# I could just use the simple ReHo datafile from rest, which should have the same NA's for the parcels applied
file<-read.csv("n1601_glasserReHoValues_20170509.csv")
#try instead with the ALFF file, just to make sure it lines up the same
file<-read.csv("n1601_glasserALFFValues_20170509.csv")
file<-right_join(file, subjlist, by ="bblid")
dim(file)
head(file)
colnames(file)
summary(file, na.rm=FALSE)
col<-colMeans(file)
row<-rowMeans(file)

summary(file$rest_glasser_reho_Right_V3B)
summary(file$rest_glasser_reho_Left_V3B)
summary(file$rest_glasser_reho_Right_52)

#count the NAs in the rows (subjects)
count_na <- function(x) sum(is.na(x))

file<-file %>% mutate(means = rowMeans(., na.rm = T),
         count_na = apply(., 1, count_na))
summary(file$count_na)
sum(file$count_na)

#Count the NAs in the columns
na_count <-sapply(file, function(y) sum(length(which(is.na(y)))))
sum(na_count)

#Who are these subjects?
a<-file$bblid[is.na(file$rest_glasser_reho_Right_V3B) == TRUE]
b<-file$bblid[is.na(file$rest_glasser_reho_Left_V3B) == TRUE]

a<-file$bblid[is.na(file$rest_glasser_alff_Right_V3B) == TRUE]
b<-file$bblid[is.na(file$rest_glasser_alff_Left_V3B) == TRUE]


#make new subject list without them
subjlist<-subjlist[subjlist$bblid != a,]
subjlist<-subjlist[subjlist$bblid != b[1],]
subjlist<-subjlist[subjlist$bblid != b[2],]

write.csv(subjlist,"/data/jag/bassett-lab/tooleyEnviNetworks/subjectLists/n1012_healthT1RestExclude_parcels.csv",row.names=FALSE, quote=FALSE)
write.csv(subjlist,"~/Documents/bassett_lab/tooleyEnviNetworks/subjectLists/n1012_healthT1RestExclude_parcels.csv",row.names=FALSE, quote=FALSE)


#run through all the connectivity matrices, check if any subjects have NAs
whichsubsnarows <- matrix(NaN,nrow=400,ncol=1015)
whichsubsnacols <-matrix(NaN,nrow=1015,ncol=400)
na_count <- numeric(1015)

for (s in 1:length(subjlist[,2])){
  sub<- subjlist[s,2]
  onesub <- sprintf("~/Documents/bassett_lab/tooleyEnviNetworks/data/rest/restNetwork_schaefer400/Schaefer400Networks/%s_Schaefer400_network.txt", sub)
  df <- data.frame(read.table(file=onesub))
  count_na <- function(x) sum(is.na(x))
  file<-df %>% mutate(count_na_rows = apply(., 1, count_na), count_na_cols = apply(., 2, count_na))
  #na_count <-sum(sapply(df, function(y) sum(length(which(is.na(y))))))
  whichsubsnarows[,s] <- file$count_na_rows
  whichsubsnacols[s,] <- file$count_na_cols
  #na_count[s] <- na_count
}

subnas <- c(subjlist,whichsubsna)
#figure out whether these are subjectwise or parcelwise
summary(whichsubsna)

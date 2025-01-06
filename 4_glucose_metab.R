

# Clean up
rm(list=ls())

# Load lme4
library(lme4)
library(lmerTest)
library(gplots)
library(ggplot2)
library(mgcv)
library(tidyr)
library(twowaytests)
library(pROC)
library(plyr)
library(caTools)

# Set the wd
wd<-"/Users/alistairsenior/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/John O'Sullivan/4 diet microbiome Heart liver plasma"
setwd(wd)

# Load in the animal key
data<-read.csv("R/clean_data.csv")

# code up factors in cols 1:5
for(i in 1:5){
	data[,i]<-as.factor(data[,i])
}

# For plotting
setwd("plots")
pdf("Glucose_metabolism.pdf")

#####################################################
####################### OGTT ########################
#####################################################

# Isolate the body mass data
OGTT<-data[, c(1:5, grep("OGTT", names(data)))]

# Reshape from wide
ogtt_long<-gather(OGTT, var, ogtt, Physiological.data_OGTT_X0: Physiological.data_OGTT_X90)
ogtt_long$time<-as.numeric(unlist(lapply(strsplit(ogtt_long$var, "X"), "[[", 2)))

# Get means and SDs by time/group
ogtt_long$id<-paste0(ogtt_long$Diet.group, "_", ogtt_long$time)
mean_ogtt<-ddply(ogtt_long, .(id), summarise, mean=mean(ogtt), se=sd(ogtt)/sqrt(length(ogtt)-1), time=time[1], Diet.group=Diet.group[1])

# Plot
a<-ggplot(mean_ogtt, aes(x=time, y=mean, group=Diet.group, color=Diet.group)) +
	geom_point(size=5) + theme_bw() + xlab("Time (mins)") + ylab("Blood Glucose") +
	theme(legend.position=c(0.85, 0.85)) + labs(subtitle="A.") + geom_path()
a

# Calculate AUC for individual animals
data$AUC<-NA
individuals<-unique(ogtt_long$ID)
for(i in 1:length(individuals)){
	
	# Isolate data for ind i
	int_i<-ogtt_long[which(ogtt_long$ID == individuals[i]),]
	
	# trapz to get AUC and add in to the main dataset
	data$AUC[which(data$ID == individuals[i])]<-trapz(int_i$time, int_i$ogtt)	
	
}

b<-ggplot(data, aes(x=Diet.group, y=AUC, group=Diet.group, color=Diet.group, fill=Diet.group)) +
	geom_violin() + theme_bw() + ylab("AUC") +
	theme(legend.position="none") + labs(subtitle="B.") + geom_point(color="black")
b

# Two way test of significance
model<-lmer(AUC ~ HF * HS + (1|Cage), data=data)
summary(model)
anova(model)

#####################################################
################## Basal Glucose ####################
#####################################################


c<-ggplot(data, aes(x=Diet.group, y=Physiological.data_basal.BGL_basal.bgl, group=Diet.group, color=Diet.group, fill=Diet.group)) +
	geom_violin() + theme_bw() + ylab("Basal Glucose") +
	theme(legend.position="none") + labs(subtitle="C.") + geom_point(color="black")
c

# Two way test of significance
model<-lmer(Physiological.data_basal.BGL_basal.bgl ~ HF * HS + (1|Cage), data=data)
summary(model)
anova(model)

#####################################################
################## Liver Size ######################
#####################################################


d<-ggplot(data, aes(x=Diet.group, y=Liver_Liver.mass_liver..mg., group=Diet.group, color=Diet.group, fill=Diet.group)) +
	geom_violin() + theme_bw() + ylab("Liver Mass (mg)") +
	theme(legend.position="none") + labs(subtitle="D.") + geom_point(color="black")
d

# Two way test of significance
model<-lmer(Liver_Liver.mass_liver..mg. ~ HF * HS + (1|Cage), data=data)
summary(model)
anova(model)

dev.off()




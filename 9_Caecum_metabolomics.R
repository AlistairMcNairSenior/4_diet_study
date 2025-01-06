

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

# Isolate the metabolomics data
data<-data[, c(1:5, grep("Caecum.metabolomics", names(data)))]

# Get the metabolites
results<-data.frame(metab=unlist(lapply(strsplit(names(data)[-c(1:5)], "_"), "[[", 3)), platform=unlist(lapply(strsplit(names(data)[-c(1:5)], "_"), "[[", 2)), error=0, b_hf=NA, b_hs=NA, b_int=NA, p_Fat=NA, p_Sugar=NA, p_Int=NA)

# Change director for plotting
setwd("plots")
pdf("Caecum_metabolomics.pdf")

# Loop for metabolites
for(i in 6:length(data)){
	
	# Create ith copy for analysis
	dat_i<-data
	dat_i$y<-data[,i]
	
	# Drop missing data
	drop<-which(is.na(dat_i$y) == T)
	if(length(drop) > 0){
		dat_i<-dat_i[-drop,]
	}
	
	# Try Analyse by lmer - note already log2 transformed
	model<-try(gam(y ~ HF * HS + s(Cage, bs="re"), data=dat_i), silent=T)

	# Refit with lm minus RE for cage	
	if(class(model)[1] == "try-error"){
		results$error[i-5]<-1
		model<-lm(y ~ HF * HS, data=dat_i)
	}
	
	# Get significanxce of 2x2
	tab<-anova(model)
	
	# Save
	if(class(model)[1] == "gam"){
		results[i-5,-c(1:6)]<-tab$pTerms.table[,3]
		results[i-5,c(4:6)]<-tab$p.table[-1,1]
	}else{
		results[i-5,-c(1:6)]<-tab$Pr[c(1:3)]
		results[i-5,c(4:6)]<-summary(model)$coef[c(2:4),1]
	}
	
}

# John Suggests Analysis should be done by platform
platforms<-unique(results$platform)

# Check distribution of p-values
par(mfrow=c(1,3))
for(i in 1:length(platforms)){
	tag<-which(results$platform == platforms[i])
	hist(results$p_Fat[tag], main=platforms[i], xlab="p Fat")
	hist(results$p_Sugar[tag], main=platforms[i], xlab="p Sug")
	hist(results$p_Int[tag], main=platforms[i], xlab="p Fat:Sug")
}

# Reset plotting
par(mfrow=c(1,1))

# Screen for fdr
results$q_Fat<-1
results$q_Sugar<-1
results$q_Int<-1
for(i in 1:length(platforms)){
	tag<-which(results$platform == platforms[i])
	results$q_Fat[tag]<-p.adjust(results$p_Fat[tag], method="fdr")
	results$q_Sugar[tag]<-p.adjust(results$p_Sugar[tag], method="fdr")
	results$q_Int[tag]<-p.adjust(results$p_Int[tag], method="fdr")
}

# Get significant effects
results[which(results$q_Fat < 0.05),]
results[which(results$q_Sugar < 0.05),]
results[which(results$q_Int < 0.05),]

# Make a volcano'ish plot
a<-ggplot(results, aes(x=b_hf, y=-log(p_Fat), color=platform)) +
		geom_point(size=2) + theme_bw() + xlab("Effect Size Fat (lnRR))") + ylab("-log(p)") + labs(subtitle="A.")
a

# Make a volcano'ish plot
b<-ggplot(results, aes(x=b_hs, y=-log(p_Sugar), color=platform)) +
		geom_point(size=2) + theme_bw() + xlab("Effect Size Sugar (lnRR))") + ylab("-log(p)") + labs(subtitle="B.")
b

# Make a volcano'ish plot
c<-ggplot(results, aes(x=b_int, y=-log(p_Int), color=platform)) +
		geom_point(size=2) + theme_bw() + xlab("Effect Size Interaction (lnRR))") + ylab("-log(p)") + labs(subtitle="C.")
c


# Get significant effects
hits<-c(results[which(results$q_Fat < 0.05),"metab"], results[which(results$q_Sugar < 0.05),"metab"], results[which(results$q_Int < 0.05),"metab"])
plats<-c(results[which(results$q_Fat < 0.05),"platform"], results[which(results$q_Sugar < 0.05),"platform"], results[which(results$q_Int < 0.05),"platform"])
hits<-unique(paste0(plats, "_", hits))

# Get the data for plotting
for(h in 1:length(hits)){
	
	# Subset the data for the right metab
	gdat<-data.frame(x=data$Diet.group, y=NA)
	gdat$y<-data[,grep(hits[h], names(data))[1]]
	
	# Create string for qs
	qs<-paste(signif(results[which(paste0(results$platform, "_", results$metab) == hits[h]),c("q_Fat", "q_Sugar", "q_Int")], 3), collapse=", ")
	
	# Make the plot
	p<-ggplot(gdat, aes(x=x, y=y, col=x, fill=x)) +
		geom_violin() + geom_point(col="black") + xlab("Diet") + ylab(hits[h]) + labs(subtitle=paste0("Fat q, ", "Sugar q, ", "Int q,
", qs)) + 
		theme_bw() + theme(legend.position="none")
		
	print(p)	
}


dev.off()

# Save the results for clustering coefficients
setwd(wd)
write.table(results, file="R/Caecum_metab_effects.csv", sep=",", row.names=F, col.names=names(results))

# How do Liver Trans-HYP and Creatine look after dropping the outlier?
plot(data$Liver.metabolomics_HILIC_Creatine, data$Liver.metabolomics_HILIC_trans.HYP)

# Find the outlier
data_drop<-data[-which(data$Liver.metabolomics_HILIC_Creatine == min(data$Liver.metabolomics_HILIC_Creatine, na.rm=T)),]

# Re-fit the model for these metabolites
model<-gam(Liver.metabolomics_HILIC_trans.HYP ~ HF * HS + s(Cage, bs="re"), data=data_drop)
anova(model)
p.adjust(0.000124, n=length(which(results$platform == "HILIC")), method="fdr")
p.adjust(0.000284, n=length(which(results$platform == "HILIC")), method="fdr")
p.adjust(0.005482, n=length(which(results$platform == "HILIC")), method="fdr")
model<-gam(Liver.metabolomics_HILIC_Creatine ~ HF * HS + s(Cage, bs="re"), data=data_drop)
anova(model)
p.adjust(0.000101, n=length(which(results$platform == "HILIC")), method="fdr")
p.adjust(0.000454, n=length(which(results$platform == "HILIC")), method="fdr")
p.adjust(0.009463, n=length(which(results$platform == "HILIC")), method="fdr")

# Main effects of fat and sugar, no interaction after dropping the outlier and re-adjusting fdr.


# How does Liver thiamine look after dropping the outlier?

# Find the outlier
data_drop<-data[-which(data$Liver.metabolomics_HILIC_Thiamine == min(data$Liver.metabolomics_HILIC_Thiamine, na.rm=T)),]

# Re-fit the model for these metabolites
model<-gam(Liver.metabolomics_HILIC_Thiamine ~ HF * HS + s(Cage, bs="re"), data=data_drop)
anova(model)
p.adjust(6.45e-05, n=length(which(results$platform == "HILIC")), method="fdr")
p.adjust(7.23e-05, n=length(which(results$platform == "HILIC")), method="fdr")
p.adjust(1.18e-06, n=length(which(results$platform == "HILIC")), method="fdr")
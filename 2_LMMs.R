

# Clean up
rm(list=ls())

# Load lme4
library(lme4)
library(lmerTest)
library(gplots)

# Set the wd
wd<-"/Users/alistairsenior/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/John O'Sullivan/4 diet microbiome Heart liver plasma"
setwd(wd)

# Load in the animal key
data<-read.csv("R/clean_data.csv")
str(data)

# code up factors in cols 1:5
for(i in 1:5){
	data[,i]<-as.factor(data[,i])
}

# Make up some hists for cols 6 through k to check normality
setwd(paste0(wd, "/plots"))
pdf("hists.pdf")
for(i in 6:length(data)){
	hist(data[,i], main=names(data)[i])
}
dev.off()
# Some look a bit skewed, need to double check residuals - body weight definitely look bi-modal!

# Table to hold results from 2x2
results<-data.frame(y=names(data)[c(6:length(data))], error=0, b_hf=NA, b_hs=NA, b_int=NA, p_Fat=NA, p_Sugar=NA, p_Int=NA)

# Plots to hold model validation
pdf("resid.pdf")
par(mfrow=c(2,2))
# Loop to analyse
for(i in 6:length(data)){
	
	# Create ith copy for analysis
	dat_i<-data
	dat_i$y<-data[,i]
	
	# Drop missing data
	drop<-which(is.na(dat_i$y) == T)
	if(length(drop) > 0){
		dat_i<-dat_i[-drop,]
	}
	

	# Scale variables - helps with model fit
	dat_i$y<-(dat_i$y - mean(dat_i$y)) / sd(dat_i$y)
	
	# Try Analyse by lmer
	model<-try(lmer(y ~ HF * HS + (1|Cage), data=dat_i), silent=T)
	
	# If we get an error or singular fit refit by lm
	if(class(model) != "try-error"){
		if(isSingular(model) == TRUE){
			results$error[i-5]<-1
			model<-lm(y ~ HF * HS, data=dat_i)
		}
	}
	if(class(model) == "try-error"){
		results$error[i-5]<-1
		model<-lm(y ~ HF * HS, data=dat_i)
	}
	
	# Get significanxce of 2x2
	tab<-anova(model)
	
	# Save
	results[i-5,-c(1:5)]<-tab$Pr[c(1:3)]
	results[i-5,c(3:5)]<-summary(model)$coef[c(2:4),1]
	
	# Plot model
	plot(model, main=results[i,1])
	
}
dev.off()

# Test significance
sig<-((results[,-c(1:5)] < 0.05) + (results[,-c(1:5)] < 0.05/nrow(results)))
row.names(sig)<-results$y

# plot heat map of signficance
pdf("Heatmap_p.pdf", height=30, width=5)
heatmap.2(sig, dendrogram="none", Colv=FALSE, trace="none", key=F, lwid=c(0.5,11), lhei=c(0.5,11), margins=c(10,15), srtRow=15, cexRow=0.4)
dev.off()

# OK, so based on that, thinking strip out the body mass, composition, food intake and glucose metabolic health and profile effect of diet on those in depth - largely fat driven







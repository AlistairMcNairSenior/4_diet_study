

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

# Isolate the SCFAs data
data<-data[, c(1:5, grep("SCFA", names(data)))]
head(data)

# Get the metabolites
results<-data.frame(metab=unlist(lapply(strsplit(names(data)[-c(1:5)], "_"), "[[", 3)), sample=unlist(lapply(strsplit(names(data)[-c(1:5)], "_"), "[[", 2)), error=0, b_hf=NA, b_hs=NA, b_int=NA, p_Fat=NA, p_Sugar=NA, p_Int=NA)

# Lets make a few plots just to look at the data
setwd("plots")
pdf("SCFAs.pdf")

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
	
	# log variables - turns coefs in lnRR
	#dat_i$y<-log(dat_i$y)
	
	# Try Analyse by lmer
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

# Do analyses and corrections sample site wise
samples<-unique(results$sample)

# Check distribution of p-values
par(mfrow=c(1,3))
for(i in 1:length(samples)){
	tag<-which(results$sample == samples[i])
	hist(results$p_Fat[tag], main=samples[i], xlab="p Fat")
	hist(results$p_Sugar[tag], main=samples[i], xlab="p Sug")
	hist(results$p_Int[tag], main=samples[i], xlab="p Fat:Sug")
}

# Reset plotting
par(mfrow=c(1,1))


# Screen for fdr
results$q_Fat<-1
results$q_Sugar<-1
results$q_Int<-1
for(i in 1:length(samples)){
	tag<-which(results$sample == samples[i])
	results$q_Fat[tag]<-p.adjust(results$p_Fat[tag], method="fdr")
	results$q_Sugar[tag]<-p.adjust(results$p_Sugar[tag], method="fdr")
	results$q_Int[tag]<-p.adjust(results$p_Int[tag], method="fdr")
}


# Get significant effects - still think we might have multiple isotopes here and need to reduce to one - to discuss in the meeting
results[which(results$q_Fat < 0.05),]
results[which(results$q_Sugar < 0.05),]
results[which(results$q_Int < 0.05),]

write.table(results, file="SCFA.Effects.csv", sep=",", row.names=F, col.names=names(results))

# Make a volcano'ish plot
a<-ggplot(results, aes(x=b_hf, y=-log(p_Fat), color=sample)) +
		geom_point(size=2) + theme_bw() + xlab("Effect Size Fat (lnRR))") + ylab("-log(p)") + labs(subtitle="A.")
a

# Make a volcano'ish plot
b<-ggplot(results, aes(x=b_hs, y=-log(p_Sugar), color=sample)) +
		geom_point(size=2) + theme_bw() + xlab("Effect Size Sugar (lnRR))") + ylab("-log(p)") + labs(subtitle="B.")
b

# Make a volcano'ish plot
c<-ggplot(results, aes(x=b_int, y=-log(p_Int), color=sample)) +
		geom_point(size=2) + theme_bw() + xlab("Effect Size Interaction (lnRR))") + ylab("-log(p)") + labs(subtitle="C.")
c

# There are main effects of fat on faecal butyric acid 5, priopionate 1 13C, propionic acid 2:4, and propionic.acid.IS.3

d<-ggplot(data, aes(x=Diet.group, y=SCFAs_Faecal_GT.butyric.acid, col=Diet.group, fill=Diet.group)) +
	geom_violin() + geom_point(col="black") + xlab("Diet") + ylab("Log(Butyric Acid) [Faecal]") + labs(subtitle="D. (Main Effect of Fat)") + 
	theme_bw() + theme(legend.position="none")
 d

e<-ggplot(data, aes(x=Diet.group, y=SCFAs_Faecal_GT.propionic.acid, col=Diet.group, fill=Diet.group)) +
	geom_violin() + geom_point(col="black") + xlab("Diet") + ylab("Log(Propionic Acid) [Faecal]") + labs(subtitle="E. (Main Effect of Fat)") + 
	theme_bw() + theme(legend.position="none")
e

f<-ggplot(data, aes(x=Diet.group, y=SCFAs_Faecal_GT.valeric.acid, col=Diet.group, fill=Diet.group)) +
	 geom_violin() + geom_point(col="black") + xlab("Diet") + ylab("Log(Valeric Acid) [Faecal]") + labs(subtitle="F. (Main Effect of Fat)") + 
	 theme_bw() + theme(legend.position="none")
f


dev.off()





# Some of the key variables that seem to have been affected by diet are food intake, body mass, composition, basal glucose and L. Homoserine
# Lets check correlations between SCFAs and these

setwd(wd)

# Load in the whole animals data
data<-read.csv("R/clean_data.csv")
names(data)

mat<-cor(data[,c("Physiological.data_Body.weight_X30", "Physiological.data_Food.intake_kcal_X30","Physiological.data_.fatmass_X.fat.mass", "Liver_Liver.mass_liver..mg.", "Physiological.data_basal.BGL_basal.bgl", "Heart.metabolomics_HILIC_L.Homoserine", "SCFAs_Faecal.ratios_Valeric.acid.4.butyric.acid.IS.3")], use="pairwise.complete.obs")

setwd("plots")
pdf("SCFAs_correlations.pdf")

heatmap.2(mat,trace="none", margins=c(15,15), srtRow=15, cexRow=0.75, srtCol=75, cexCol=0.75, col=colorRampPalette(c("blue", "white", "red")))

# Only really convincing correlation for SCFA might be with fat mass????
a<-ggplot(data, aes(x=SCFAs_Faecal.ratios_Valeric.acid.4.butyric.acid.IS.3, y=Physiological.data_.fatmass_X.fat.mass, color=Diet.group)) +
	geom_point()
a

# Bit of an ecological effect going on there.
model<-try(gam(Physiological.data_.fatmass_X.fat.mass ~ SCFAs_Faecal.ratios_Valeric.acid.4.butyric.acid.IS.3 + as.factor(HF) + s(Cage, bs="re"), data=data), silent=T)
summary(model)
# Yeah nada ...

dev.off()

# Check the cage effect on Caecal SCFAs
head(data)

m1<-lmer(SCFAs_Caecum_GT.acetic.acid ~ HF * HS + (1|Cage), data=data)
summary(m1)
0.0004676 / (0.0004676 + 0.0011907)
m1_0<-lm(SCFAs_Caecum_GT.acetic.acid ~ HF * HS, data=data)
AIC(m1_0, m1)
# Not much better with cage

m2<-lmer(SCFAs_Caecum_GT.butyric.acid ~ HF * HS + (1|Cage), data=data)
summary(m2)
# ERROR - 0 cage

m3<-lmer(SCFAs_Caecum_GT.Caproic.acid ~ HF * HS + (1|Cage), data=data)
summary(m3)
 0.1344 / (0.1344 + 0.2695)
m3_0<-lm(SCFAs_Caecum_GT.Caproic.acid ~ HF * HS, data=data)
AIC(m3_0, m3)
# Not good evidence

m4<-lmer(SCFAs_Caecum_GT.propionic.acid ~ HF * HS + (1|Cage), data=data)
summary(m4)
# ERROR - 0 cage

m5<-lmer(SCFAs_Caecum_GT.valeric.acid ~ HF * HS + (1|Cage), data=data)
summary(m5)
# ERROR - 0 cage


m3<-lmer(SCFAs_Faecal_GT.Caproic.acid ~ HF * HS + (1|Cage), data=data)
summary(m3)



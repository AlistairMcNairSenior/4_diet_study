

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
pdf("Intake_body_composition.pdf")

#####################################################
##################### Body Mass #####################
#####################################################

# Isolate the body mass data
body_mass<-data[, c(1:5, grep("Body.weight", names(data))[-1])]

# Reshape from wide
body_long<-gather(body_mass, var, mass, Physiological.data_Body.weight_X0:Physiological.data_Body.weight_X30)
body_long$age<-as.numeric(unlist(lapply(strsplit(body_long$var, "X"), "[[", 2)))

# Plot
a<-ggplot(body_long, aes(x=age, y=mass, group=Diet.group, color=Diet.group)) +
	geom_point() + geom_smooth() + theme_bw() + xlab("Age (Days)") + ylab("Body Mass (g)") +
	theme(legend.position=c(0.1, 0.85)) + labs(subtitle="A.")
a

# Test
model<-gam(mass ~ HF * HS + s(age, by=Diet.group) + s(ID, bs="re") + s(Cage, bs="re"), data=body_long)
summary(model)
anova(model)
# Very clear main effect of fat - no sugar nor interaction

#####################################################
################## Body Composition #################
#####################################################

# Plots of these two
b<-ggplot(data, aes(x=Diet.group, y=Physiological.data_.fatmass_X.fat.mass, group=Diet.group, color=Diet.group, fill=Diet.group)) +
	geom_violin() + theme_bw() + xlab("Group") + ylab("Fat Mass (%)") +
	theme(legend.position="none") + geom_point(color="black") + labs(subtitle="B.")
b

# Two-way test under heteroscedasticity
fat_mass_test<-gpTwoWay(Physiological.data_.fatmass_X.fat.mass ~ HF * HS, data=data)

# Also classic LMM with cage effect
fat_mass_test2<-lmer(Physiological.data_.fatmass_X.fat.mass ~ HF * HS + (1|Cage), data=data)
anova(fat_mass_test2)

# Clear main effect of fat


#####################################################
################### Food Intake (g) #################
#####################################################

# Isolate the body mass data
intake_g<-data[, c(1:5, grep("Food.intake_g", names(data)))]

# Need to reduce to one per cage to avoid pseudo replicating
intake_g<-intake_g[match(unique(intake_g$Cage), intake_g$Cage),]

# Reshape
intake_long<-gather(intake_g, var, mass, Physiological.data_Food.intake_g_X1:Physiological.data_Food.intake_g_X30)
intake_long$age<-as.numeric(unlist(lapply(strsplit(intake_long$var, "X"), "[[", 2)))

# Ages 1 and 30, we have missing data for some animals - lets drop those ages
intake_long[which(is.na(intake_long$mass) == T),]
intake_long<-intake_long[-which(intake_long$age == 1 | intake_long$age == 30),]

# Plot
c<-ggplot(intake_long, aes(x=age, y=mass, group=Diet.group, color=Diet.group)) +
	geom_point() + geom_smooth() + theme_bw() + xlab("Age (Days)") + ylab("Food Intake (g)") + labs(subtitle="C.")
c

# Test
model<-gam(mass ~ HF * HS + s(age, by=Diet.group, k=3) + s(ID, bs="re") + s(Cage, bs="re"), data=intake_long)
summary(model)
anova(model)
# Not super convincing, but higher fat???

c2<-ggplot(intake_long, aes(x=Diet.group, y=mass, group=Diet.group, color=Diet.group, fill=Diet.group)) +
	geom_violin() + theme_bw() + ylab("Food Intake (g)") + geom_point(color='black') + labs(subtitle="C.")
c2
# Yeah maybe HF, or maybe all non-chow???

#####################################################
############## Energy Intake (kcal) #################
#####################################################

# Isolate the body mass data
intake_cal<-data[, c(1:5, grep("Food.intake_kcal", names(data)))]

# Need to reduce to one per cage to avoid pseudo replicating
intake_cal<-intake_cal[match(unique(intake_cal$Cage), intake_cal$Cage),]

# Reshape
intake_long<-gather(intake_cal, var, cal, Physiological.data_Food.intake_kcal_X1:Physiological.data_Food.intake_kcal_X30)
intake_long$age<-as.numeric(unlist(lapply(strsplit(intake_long$var, "X"), "[[", 2)))

# Ages 1 and 30, we have missing data for some animals - lets drop those ages
intake_long[which(is.na(intake_long$cal) == T),]
intake_long<-intake_long[-which(intake_long$age == 1 | intake_long$age == 30),]

# Plot
d<-ggplot(intake_long, aes(x=age, y=cal, group=Diet.group, color=Diet.group)) +
	geom_point() + geom_smooth() + theme_bw() + xlab("Age (Days)") + ylab("Energy Intake (kcal)") + labs(subtitle="D.")
d

# Test
model<-gam(cal ~ HF * HS + s(age, by=Diet.group, k=4) + s(ID, bs="re") + s(Cage, bs="re"), data=intake_long)
summary(model)
anova(model)
# Not super convincing, but higher fat???

d2<-ggplot(intake_long, aes(x=Diet.group, y=cal, group=Diet.group, color=Diet.group, fill=Diet.group)) +
	geom_violin() + theme_bw() + ylab("Energy Intake (kcal)") + geom_point(color='black') + labs(subtitle="D.")
d2
# Looks like a bit of an interaction, but does not shake out statistically 

dev.off()

# So to conclude - High Fat increases body mass and adiposity via increased energy intake and some evidence for elevated dry weight food intake. Little evidence for sugar affecting.

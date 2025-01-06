

# Clean up
rm(list=ls())

# Load packages
library(tidyr)

# Set the wd
wd<-"/Users/alistairsenior/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/John O'Sullivan/4 diet microbiome Heart liver plasma"
setwd(wd)

# Load in the animal key
animals<-read.csv("Raw_data/Animal_Key.csv")

# Code up the 2x2 (ish) design
animals$HF<-0
animals$HS<-0
animals$HF[grep("HF", animals$Diet.group)]<-1
animals$HS[grep("HS", animals$Diet.group)]<-1

# Now read in all the phenotypes
files<-dir("Raw_data")[-1]
for(i in 1:length(files)){
	
	# File path for the raw data
	file_i<-paste0("Raw_data/", files[i])
	dat<-read.csv(file_i)
	
	# Tag to the write animals
	tag<-match(animals$ID, dat[,1])
	new<-as.data.frame(dat[tag,-1])
	names(new)<-paste0(strsplit(files[i], ".csv")[[1]], "_", names(dat)[-1])
	animals<-cbind(animals, new)
	
}

# Drop redundant cols on ID etc
animals<-animals[,-grep("diet", names(animals))]

# Save the processed data in the R folder for analysis
write.table(animals, file=paste0("R/clean_data.csv"), sep=",", row.names=F, col.names=names(animals))

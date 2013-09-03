#---------------------------------------------------------------------------------
# Name - Agilent.preprocess.mRNA.r
# Desc - Pre-processing script for Agilent mRNA expression data
# Author - Robinder Gauba (rg353@georgetown.ed; robinaguba@gmail.com)
# Source code - https://github.com/rg353/MyScripts/blob/master/Agilent.preprocess.mRNA.r
#---------------------------------------------------------------------------------

# Call libraries
library('Agi4x44PreProcess')
library('limma')

# Set working directory
setwd("C:/Users/xxx/")

# Read expression files
dat<-read.maimages(files=dir(), source="agilent")
#or
targets=read.targets(infile="targets.txt") # G4112F version
dd=read.AgilentFE(targets,makePLOT=FALSE) # G4112F version

# Background Correction & Normalization 
ddNORM=BGandNorm(dat,BGmethod='normexp',NORMmethod='quantile',foreground='MeanSignal',background='BGMedianSignal',offset=50,makePLOTpre=FALSE,makePLOTpost=FALSE)

# Summarization Correction

ddPROC=summarize.probe(ddNORM,makePLOT=F, targets)
esetPROC=build.eset(ddPROC,targets,makePLOT=F,annotation.package="hgug4112a.db")

# Saving expression matrix
write.table(exprs(esetPROC),file="dataMatrix_GSE17047.txt",sep="\t")

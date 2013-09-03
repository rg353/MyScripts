#---------------------------------------------------------------------------------
# Name - Agilent.preprocess.miRNA.r
# Desc - Pre-processing script for Agilent miRNA expression data
# Author - Robinder Gauba (rg353@georgetown.ed; robinaguba@gmail.com)
# Source code - https://github.com/rg353/MyScripts/blob/master/Agilent.preprocess.miRNA.r
#---------------------------------------------------------------------------------

# Call libraries
library(limma)
library(AgiMicroRna)

# Set working directory
setwd("C:/Users/xxx/")


# Read expression files
targets.micro=readTargets(infile="targets.txt",verbose=TRUE)
#or
dd=readMicroRnaAFE(targets.micro,verbose=TRUE)

# Reading & Signal Extraction

ddTGS.AFE = tgsMicroRna(dd,half=TRUE, makePLOT = FALSE,verbose = TRUE)  

# Normalization 
ddNORM.AFE = tgsNormalization(ddTGS.AFE, "quantile", makePLOTpre = FALSE,makePLOTpost = FALSE, targets.micro, verbose = TRUE)

# Removal of Controls   
ddPROC.AFE = filterMicroRna(ddNORM.AFE, dd,  control = TRUE,IsGeneDetected = TRUE,wellaboveNEG = FALSE, limIsGeneDetected = 75, limNEG = 25,makePLOT = FALSE, targets.micro, verbose = TRUE)

# Creating eset object
esetPROC.AFE = esetMicroRna(ddPROC.AFE, targets.micro, makePLOT = FALSE, verbose = TRUE)

# Saving expression matrix
writeEset(esetPROC.AFE, ddPROC.AFE, targets.micro, verbose = TRUE)

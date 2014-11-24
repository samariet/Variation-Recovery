##Load libraries
library(lumi)
library(ggplot2)
library(gplots)
library(TeachingDemos)

setwd("~/Lab/Variation Recovery/from minal/gene expression data for Samantha/gene expression data for Samantha")

##Create object with name of data file:
data = c('exp_combined.txt')
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)

#Take only the probes that have a detection p-value<.05 in at least one individual
hist(data.norm.all@assayData$detection)
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05)
detect.ind.all <- which(detect_quant.all > 0)

#With this threshold 28,663 probes out of 47,315 are expressed
expr_quant.all <- data.norm.all@assayData$exprs[detect.ind.all,]

###Subset expression by Darren's good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID

## Convert from factor to character
probes = as.character(goodprobes$probeID)
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
dim(expr_quant.all.clean)  
dim(expr_quant.all)

# 20,304 probes of the original 32,562 "good" probes are in this data set
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)


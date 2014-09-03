##Load libraries
library(lumi)
library(ggplot2)
library(gplots)
library(TeachingDemos)
library(VennDiagram)
library(xlsx)
library(GSEAlm)
library(biomaRt)
library(piano)
library(cluster)
# hi
# lol

setwd("~/Lab/Variation Recovery")

##Create object with name of data file:
data = c('YGilad-ST-May18-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

plot(data.lumi, what='boxplot')

##LCL only
data.lumi = data.lumi[,c(1,3,5,7,9,13,15,17,19,21,23,25,27,29,31,33,35)]
stems = c(1,3,5,7,9,13,15,17,19,21,23,25,27,29,31,33,35)

### NORMALIZATION ###
# Convert raw Illumina probe intensities to expression values
# Corrects background, log2 stabiliizes variance, and quantile normalize
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)

#Take only the probes that have a detection p-value<.05 in at least one individual
hist(data.norm.all@assayData$detection)
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05)
detect.ind.all <- which(detect_quant.all > 0)

#With this threshold 28,663 probes out of 47,315 are expressed
expr_quant.all <- data.norm.all@assayData$exprs[detect.ind.all,]

###Find the column that is lumi_ID in feature data usually column 1
head(data.norm.all@featureData[[5]]) ## Should read: [1] "ILMN_1762337" "ILMN_2055271" "ILMN_1736007" "ILMN_2383229" "ILMN_1806310" "ILMN_1779670"....

###Convert expr_quant rownames to these lumi IDs
rownames(expr_quant.all)=data.norm.all@featureData[[5]][detect.ind.all]

###Subset expression by Darren's good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
probes = as.character(goodprobes$probeID) ## Convert from factor to character
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
dim(expr_quant.all.clean)
dim(expr_quant.all) #Look at this, you lose a lot of data here. Oh well.
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)

#Label your columns by sample name
samplenames = read.delim('YGilad-ST sample names switched.txt', header=TRUE)
colnames(expr_quant.all) = samplenames[(stems),1]

#Define covariates
indiv = samplenames[,2]
type = samplenames[,3]
ID = samplenames[,4]
repr_batch = samplenames[,5]
array_batch = samplenames[,6]
gender = samplenames[,7]
extr_batch = samplenames[,8]
extr_date = samplenames[,9]

#Converted categorical covariates to a factor so they are levels. This will allow you to do math on them.
indiv.f = as.factor(indiv)
type.f = as.factor(type)
ID.f=as.factor(ID)
repr_batch.f = as.factor(repr_batch)
array_batch.f = as.factor(array_batch)
gender.f=as.factor(gender)
extr_batch.f = as.factor(extr_batch)
extr_date.f = as.factor(extr_date)

#Subset stem cells and put your covariates in a list
indiv.fs = indiv.f[stems]
type.fs <- type.f[stems]
ID.fs <- ID.f[stems]
repr_batch.fs <- repr_batch.f[stems]
array_batch.fs <- array_batch.f[stems]
gender.fs <- gender.f[stems]
extr_batch.fs <- extr_batch.f[stems]
extr_date.fs <- extr_date.f[stems]
covars<-list(array_batch.fs,indiv.fs,gender.fs,extr_date.fs,extr_batch.fs) #leave out repr batch for lcls, will bug later
covars_names <- factor(c("array_batch.fs","indiv.fs","gender.fs","extr_date.fs","extr_batch.fs"))

#Subtract out mean. Make sure this is what you want to do.
expr_nice <- t(apply(expr_quant.all,1,function(x){x-mean(x)}))
expr_stem <- expr_nice

## This next section subsets your data to include only 1 probe per gene.. 

expr_genes <- expr_stem

## Convert probe IDs to gene names using Darren's file (has HGNC & ensembl)
gene_names=c()
for(i in 1:dim(expr_stem)[1]){ #what is [1] doing?
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_quant.all)[i],7])) #creates a list of gene names the same length as probe list
}
rownames(expr_genes)=gene_names
Unique_genes = unique(rownames(expr_genes))
length(Unique_genes) #this will be the length of your final data set. One data point per gene.

## This loop will give the most 3' value for multiple probes within the same gene. 
expr_gene = matrix(NA, ncol=length(stems), nrow=length(Unique_genes))
i=0
for(gene in Unique_genes){
  
  i = i+1
  
  currRows = which(Unique_genes == gene) #all probes of that gene
  if(length(currRows)>1){
    if(expr_genes[goodprobes[1],6]=="+"){ #wait so all probes we're keeping are on the same strand?
      keepRow = currRows[which.max(goodprobes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  expr_gene[i,] = expr_stem[keepRow,]
  
} 

dim(expr_gene)
rownames(expr_gene) = Unique_genes
colnames(expr_gene) = colnames(expr_stem)

#cor.q <- cor(expr_gene,method="pearson", use="complete.obs")
#heatmap.2(cor.q, main="Gene Expression Correlation", key=T, revC=T, density.info="none", trace="none")

##Regress out array batch and add the intercept back in
abatch.residual.int.g = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_stem))
rownames(abatch.residual.int.g) = rownames(expr_gene)
colnames(abatch.residual.int.g) = colnames(expr_gene)
for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ array_batch.fs)
  abatch.residual.int.g[i,] = resid(model) + model$coefficients[1]
}
cor.abatch.int.g <- cor(abatch.residual.int.g,method="pearson", use="complete.obs")
heatmap.2(cor.abatch.int.g, Rowv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))),
          Colv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))), margins=c(5,9),key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")

abatch_lcl <- abatch.residual.int.g

#Relationship between PCs and covariates for regressed data

npcs = 4
sum.PC <- prcomp(na.omit(abatch.residual.int.g))
results<-c()
for (f in covars) {
  for (i in 1:npcs)
  {
    s = summary(lm(sum.PC$rotation[,i]~f));
    results<-c(results,pf(s$fstatistic[[1]],
                          s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
               s$adj.r.squared)
  }
}
resultsM_gene_corrected <-matrix(nrow = length(covars), ncol = 2*npcs, data =
                                   results, byrow = TRUE)
rownames(resultsM_gene_corrected) = covars_names
colnames(resultsM_gene_corrected) = c("PC1 p value","PC1 R2","PC2 p value","PC2 R2","PC3 p value","PC3 R2","PC4 p value","PC4 R2")

resultsM_gene_corrected
PC_table_lcl <- resultsM_gene_corrected

#PCA for regressed data
sum.PC <- prcomp(na.omit(abatch.residual.int.g))
sumsum_gene_corrected <- summary(sum.PC)

op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
title.PC = "PCA of Gene Expression in lcls"
sum.PC <- prcomp(na.omit(abatch.residual.int.g), scale=TRUE)
sumsum <- summary(sum.PC)

#prints out plots in c(#rows, #columns)
color = indiv.fs
par(mfrow = c(1,1),oma=c(0,0,2,0)) 
plot(c(1:17),sum.PC$rotation[,1],cex=1.5,col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:17),sum.PC$rotation[,1], samplenames[,2], cex = 0.55, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC$rotation[,1], sum.PC$rotation[,i],labels=samplenames[(stems),2], cex = 0.8, pos=3)   
}  

#Making dendrograms
cor_lcl <- cor(abatch_lcl, method="pearson")
dis  <- 1-cor_lcl
distance <- as.dist(dis)
hc <- hclust(distance)
color <- indiv.fs
plot(hc, main = "Cluster Dendrogram: lcls", cex.main = 1.5 , col = "#487AA1", col.main = "#45ADA8", col.lab = "lightcoral", col.axis = "#F38630", lwd = 3, lty = 1, sub = "", hang = -1, axes = FALSE)
axis(side = 2, at = seq(0.0, 1.4, .2), col = "#F38630", labels = FALSE, lwd = 2)
# add text in margin
mtext(seq(0, 1.4, .2), side = 2, at = seq(0, 1.4, .2), line = 1, col = "#A38630", las = 2)

#Within individual pearson correlation coefficient
cor_2l <- c(cor[1,6],cor[1,12],cor[6,12])
cor_5l <- c(cor[2,7],cor[2,13],cor[7,13])
cor_6l <- c(cor[3,8],cor[8,14],cor[3,14])
cor_9l <- c(cor[4,9],cor[4,15],cor[9,15])
cor_10l <- c(cor[10,16])
cor_14l <- c(cor[5,11],cor[5,17],cor[11,17])

cor_within_l <- cbind(cor_2l,cor_5l,cor_6l,cor_9l,cor_10l,cor_14l)

cor_wmeans_l <- sapply(cor_within_l, mean)

# across individual correlation coefficients
cor_all_lcl <- c()
for (i in 1:nrow(cor_lcl)) {
  for (j in 1:ncol(cor_lcl)) {
    print(i)
    cor_all_lcl <- c(cor_all_lcl, cor_lcl[i,j])
  }
}
cor_all_lcl
boxplot(cor_all_lcl)

#Euclidean distances
dist12_ind <- c()
PC12 <- cbind(sum.PC$rotation[,1], sum.PC$rotation[,2])
#Dist by individual
ind1 <- PC12[c(1,6,12),]
ind2 <- PC12[c(2,7,13),]
ind3 <- PC12[c(3,8,14),]
ind4 <- PC12[c(4,9,15),]
ind5 <- PC12[c(10,16),]
ind6 <- PC12[c(5,11,17)]

md1 <- mean(dist(ind1))
md2 <- mean(dist(ind2))
md3 <- mean(dist(ind3))
md4 <- mean(dist(ind4))
md5 <- mean(dist(ind5))
md6 <- mean(dist(ind6))
mdL <- c(md1,md2,md3,md4,md5,md6)
mdL_mean <- mean(mdL)

#Dist all
dist12L <- dist(cbind(sum.PC$rotation[,1], sum.PC$rotation[,2]))
dist12L
mean(dist12L)

# STEMS

setwd("~/Lab/Variation Recovery")

##Create object with name of data file:
data = c('YGilad-ST-May18-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

plot(data.lumi, what='boxplot')

##LCL only
data.lumi = data.lumi[,c(2,4,6,8,10,14,16,18,20,22,24,26,28,30,32,34,36)]
stems = c(2,4,6,8,10,14,16,18,20,22,24,26,28,30,32,34,36)

### NORMALIZATION ###
# Convert raw Illumina probe intensities to expression values
# Corrects background, log2 stabiliizes variance, and quantile normalize
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)

#Take only the probes that have a detection p-value<.05 in at least one individual
hist(data.norm.all@assayData$detection)
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05)
detect.ind.all <- which(detect_quant.all > 0)

#With this threshold 28,663 probes out of 47,315 are expressed
expr_quant.all <- data.norm.all@assayData$exprs[detect.ind.all,]

###Find the column that is lumi_ID in feature data usually column 1
head(data.norm.all@featureData[[5]]) ## Should read: [1] "ILMN_1762337" "ILMN_2055271" "ILMN_1736007" "ILMN_2383229" "ILMN_1806310" "ILMN_1779670"....

###Convert expr_quant rownames to these lumi IDs
rownames(expr_quant.all)=data.norm.all@featureData[[5]][detect.ind.all]

###Subset expression by Darren's good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
probes = as.character(goodprobes$probeID) ## Convert from factor to character
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
dim(expr_quant.all.clean)
dim(expr_quant.all) #Look at this, you lose a lot of data here. Oh well.
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)

#Label your columns by sample name
samplenames = read.delim('YGilad-ST sample names switched.txt', header=TRUE)
colnames(expr_quant.all) = samplenames[(stems),1]

#Define covariates
indiv = samplenames[,2]
type = samplenames[,3]
ID = samplenames[,4]
repr_batch = samplenames[,5]
array_batch = samplenames[,6]
gender = samplenames[,7]
extr_batch = samplenames[,8]
extr_date = samplenames[,9]

#Converted categorical covariates to a factor so they are levels. This will allow you to do math on them.
indiv.f = as.factor(indiv)
type.f = as.factor(type)
ID.f=as.factor(ID)
repr_batch.f = as.factor(repr_batch)
array_batch.f = as.factor(array_batch)
gender.f=as.factor(gender)
extr_batch.f = as.factor(extr_batch)
extr_date.f = as.factor(extr_date)

#Subset stem cells and put your covariates in a list
indiv.fs = indiv.f[stems]
type.fs <- type.f[stems]
ID.fs <- ID.f[stems]
repr_batch.fs <- repr_batch.f[stems]
array_batch.fs <- array_batch.f[stems]
gender.fs <- gender.f[stems]
extr_batch.fs <- extr_batch.f[stems]
extr_date.fs <- extr_date.f[stems]
covars<-list(array_batch.fs, repr_batch.fs,indiv.fs,gender.fs,extr_date.fs,extr_batch.fs)
covars_names <- factor(c("array_batch.fs", "repr_batch.fs","indiv.fs","gender.fs","extr_date.fs","extr_batch.fs"))

#Subtract out mean. Make sure this is what you want to do.
expr_nice <- t(apply(expr_quant.all,1,function(x){x-mean(x)}))
expr_stem <- expr_nice

## This next section subsets your data to include only 1 probe per gene.. 

expr_genes <- expr_stem

## Convert probe IDs to gene names using Darren's file (has HGNC & ensembl)
gene_names=c()
for(i in 1:dim(expr_stem)[1]){ #what is [1] doing?
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_quant.all)[i],7])) #creates a list of gene names the same length as probe list
}
rownames(expr_genes)=gene_names
Unique_genes = unique(rownames(expr_genes))
length(Unique_genes) #this will be the length of your final data set. One data point per gene.

## This loop will give the most 3' value for multiple probes within the same gene. 
expr_gene = matrix(NA, ncol=length(stems), nrow=length(Unique_genes))
i=0
for(gene in Unique_genes){
  
  i = i+1
  
  currRows = which(Unique_genes == gene) #all probes of that gene
  if(length(currRows)>1){
    if(expr_genes[goodprobes[1],6]=="+"){ #wait so all probes for a gene are on the same strand?
      keepRow = currRows[which.max(goodprobes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  expr_gene[i,] = expr_stem[keepRow,]
  
} 

dim(expr_gene)
rownames(expr_gene) = Unique_genes
colnames(expr_gene) = colnames(expr_stem)

#cor.q <- cor(expr_gene,method="pearson", use="complete.obs")
#heatmap.2(cor.q, main="Gene Expression Correlation", key=T, revC=T, density.info="none", trace="none")

##Regress out array batch and add the intercept back in
abatch.residual.int.g = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_stem))
rownames(abatch.residual.int.g) = rownames(expr_gene)
colnames(abatch.residual.int.g) = colnames(expr_gene)
for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ array_batch.fs)
  abatch.residual.int.g[i,] = resid(model) + model$coefficients[1]
}
cor.abatch.int.g <- cor(abatch.residual.int.g,method="pearson", use="complete.obs")
heatmap.2(cor.abatch.int.g, Rowv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))),
          Colv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))), margins=c(5,9),key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")

abatch_stem <- abatch.residual.int.g

#Relationship between PCs and covariates for regressed data

npcs = 4
sum.PC <- prcomp(na.omit(abatch.residual.int.g))
results<-c()
for (f in covars) {
  for (i in 1:npcs)
  {
    s = summary(lm(sum.PC$rotation[,i]~f));
    results<-c(results,pf(s$fstatistic[[1]],
                          s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
               s$adj.r.squared)
  }
}
resultsM_gene_corrected <-matrix(nrow = length(covars), ncol = 2*npcs, data =
                                   results, byrow = TRUE)
rownames(resultsM_gene_corrected) = covars_names
colnames(resultsM_gene_corrected) = c("PC1 p value","PC1 R2","PC2 p value","PC2 R2","PC3 p value","PC3 R2","PC4 p value","PC4 R2")

resultsM_gene_corrected
PC_table_stem <- resultsM_gene_corrected

#PCA for regressed data
sum.PC <- prcomp(na.omit(abatch.residual.int.g))
sumsum_gene_corrected <- summary(sum.PC)

op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
title.PC = "PCA of Gene Expression in iPSCs"
sum.PC <- prcomp(na.omit(abatch.residual.int.g), scale=TRUE)
sumsum <- summary(sum.PC)

#prints out plots in c(#rows, #columns)
color = indiv.fs
par(mfrow = c(1,1),oma=c(0,0,2,0)) 
plot(c(1:17),sum.PC$rotation[,1],cex=1.5,col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:17),sum.PC$rotation[,1], samplenames[,2], cex = 0.55, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC$rotation[,1], sum.PC$rotation[,i],labels=samplenames[(stems),2], cex = 0.8, pos=3)   
}  

PC_12 <- c(sum.PC$rotation[,1], sum.PC$rotation[,2])
PC_12M <- matrix(PC_12, ncol=2)
dim(PC_12M)
plot(PC_12M)

fit <- kmeans(PC_12M, 6)
#get cluster means
aggregate(PC_12M,by=list(fit$cluster),FUN=mean)
# append cluster assignment
PC_12DF <- data.frame(PC_12M, fit$cluster) 

clusplot(PC_12M, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

PC_12 <- c(sum.PC$rotation[,1], sum.PC$rotation[,2], sum.PC$rotation[,3])
PC_12M <- matrix(PC_12, ncol=3)
dim(PC_12M)
plot(PC_12M)

fit <- kmeans(PC_12M, 6)
#get cluster means
aggregate(PC_12M,by=list(fit$cluster),FUN=mean)
# append cluster assignment
PC_12DF <- data.frame(PC_12M, fit$cluster) 

clusplot(PC_12M, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

#Making dendrograms
cor <- cor(abatch_stem, method="pearson")
dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)
color <- indiv.fs
plot(hc, main = "Cluster Dendrogram: iPSCs", cex.main = 1.5 , col = "#487AA1", col.main = "#45ADA8", col.lab = "lightcoral", col.axis = "#F38630", lwd = 3, lty = 1, sub = "", hang = -1, axes = FALSE)
axis(side = 2, at = seq(0.0, 1.4, .2), col = "#F38630", labels = FALSE, lwd = 2)
# add text in margin
mtext(seq(0, 1.4, .2), side = 2, at = seq(0, 1.4, .2), line = 1, col = "#A38630", las = 2)

#Within individual pearson correlation coefficient in stems REMEMBER 9 AND 6 ARE SWITCHED SO THESE LABELS AREN'T THE SAME AS LCLS
cor_2 <- c(cor[1,6],cor[1,12],cor[6,12])
cor_5 <- c(cor[2,7],cor[2,13],cor[7,13])
cor_6 <- c(cor[4,8],cor[4,14],cor[8,14])
cor_9 <- c(cor[3,9],cor[3,15],cor[9,15])
cor_10 <- c(cor[10,16])
cor_14 <- c(cor[5,11],cor[5,17],cor[11,17])

cor_within <- cbind(cor_2,cor_5,cor_6,cor_9,cor_10,cor_14)
cor_wmeans <- sapply(cor_within, mean)
cor_both <- cbind(cor_within, cor_within_l)
cor_bmeans <- cbind(cor_wmeans, cor_wmeans_l)
boxplot(cor_wmeans)
boxplot(cor_both)
boxplot(cor_bmeans, main="Within Individual Pearson correlation coefficients for Stem cells and LCLs")

# across individual correlation coefficients
cor_all_lcl <- c()
for (i in 1:nrow(cor_lcl)) {
  for (j in 1:ncol(cor_lcl)) {
    print(i)
    cor_all_lcl <- c(cor_all_lcl, cor_lcl[i,j])
  }
}
cor_all_lcl
boxplot(cor_all_lcl)

#pvalues for each gene association with individual

pvals_stems <- c()
for (i in 1:nrow(abatch_stem)) {
  s <- summary(lm(abatch_stem[i,]~indiv.fs));
  pvals_stems <- c(pvals_stems,pf(s$fstatistic[[1]], s$fstatistic[[2]], s$fstatistic[[3]], lower.tail = FALSE))
}

pvals_lcls <- c()
for (i in 1:nrow(abatch_lcl)) {
  s <- summary(lm(abatch_lcl[i,]~indiv.fs));
  pvals_lcls <- c(pvals_lcls,pf(s$fstatistic[[1]], s$fstatistic[[2]], s$fstatistic[[3]], lower.tail = FALSE))
}

hist(pvals_stems)
hist(pvals_lcls)

plot(stem_var_rat, pvals_stems)
plot(lcl_var_rat, pvals_lcls)

#count significant genes
length(pvals_stems[pvals_stems<.001])
length(pvals_lcls[pvals_lcls<.001])

# Euclidean distance within and between indvl of PC projections 1 & 2

dist12_ind <- c()
PC12 <- cbind(sum.PC$rotation[,1], sum.PC$rotation[,2])
#Dist by individual
ind1 <- PC12[c(1,6,12),]
ind2 <- PC12[c(2,7,13),]
ind3 <- PC12[c(4,8,14),]
ind4 <- PC12[c(3,9,15),]
ind5 <- PC12[c(10,16),]
ind6 <- PC12[c(5,11,17)]

md1 <- mean(dist(ind1))
md2 <- mean(dist(ind2))
md3 <- mean(dist(ind3))
md4 <- mean(dist(ind4))
md5 <- mean(dist(ind5))
md6 <- mean(dist(ind6))
mds <- c(md1,md2,md3,md4,md5,md6)
mds_mean <- mean(mds)

#Dist all
dist12 <- dist(cbind(sum.PC$rotation[,1], sum.PC$rotation[,2]))
dist12
mean(dist12)

t.test(dist12, mds)
t.test(dist12L, mdL)

## Variance within vs among indvls

#Between samples coefficient of variation
cv <- function(x) (sd(x)/mean(x))
CV_S <- apply(abatch_stem, 1, cv)
CV_L <- apply(abatch_lcl, 1, cv)

#Density plots of CVs
cv_all <- data.frame(coefvar=c(CV_L, CV_S), type = rep(c("LCL", "Stem"), times=c(length(CV_L),length(CV_S))))
ggplot(cv_all, aes(x=coefvar, fill=type)) + geom_density(alpha=0.5) +xlim(0,2)+xlab("Coefficients of Variation") + ggtitle("Gene Expression: Coefficients of Variation") + theme(legend.position=c(.75,.75)) +annotate(geom = "text", label="**p-value = 2.2e-16", x=1.9, y=.9)+theme(text = element_text(size=23)) 

## Variance within vs among indvls
expr_LCL <- abatch_lcl

var_1 <- apply(expr_LCL[,c(1,6,12)],1,var)
var_2 <- apply(expr_LCL[,c(2,7,13)],1,var)
var_3 <- apply(expr_LCL[,c(3,8,14)],1,var)
var_4 <- apply(expr_LCL[,c(4,9,15)],1,var)
var_5 <- apply(expr_LCL[,c(10,16)],1,var)
var_6 <- apply(expr_LCL[,c(5,11,17)],1,var)

var_win_lcl <- cbind(var_1,var_2,var_3,var_4,var_5,var_6)
var_lcl <- apply(expr_LCL,1,var)

#Calculate ratio of variance btw and within individuals for every gene
var_ratio_lcl <- vector()
for(i in 1:length(var_lcl)) {
  var_ratio_lcl <- c(var_ratio_lcl, (mean(var_win_lcl[i,]))/(var_lcl[i]))  
}
length(var_ratio_lcl)

#Do same for stem cells REMEMBER 9 AND 6 ARE SWITCHED SO THESE LABELS AREN'T THE SAME AS LCLS
expr_stem <- abatch_stem #if you've stored it this way

var_1s <- apply(expr_stem[,c(1,6,12)],1,var)
var_2s <- apply(expr_stem[,c(2,7,13)],1,var)
var_3s <- apply(expr_stem[,c(4,8,14)],1,var)
var_4s <- apply(expr_stem[,c(3,9,15)],1,var)
var_5s <- apply(expr_stem[,c(10,16)],1,var)
var_6s <- apply(expr_stem[,c(5,11,17)],1,var)

var_win <- cbind(var_1s,var_2s,var_3s,var_4s,var_5s,var_6s)
var_stem <- apply(expr_stem,1,var) #calculate variance for each gene across all individuals
#rownames(var_stem) <- rownames(expr_stem)

#Calculate ratio of variance btw and within individuals for every gene
var_ratio <- vector()
for(i in 1:length(var_stem)) {
  var_ratio <- c(var_ratio, (mean(var_win[i,]))/(var_stem[i]))  
}
length(var_ratio)

# Order by variance ratio. Plot by gene CHANGE INDIVIDUAL IDENTIFIERS!!!!

color=indiv.fs
stem_bind <- data.frame(expr_stem, var_ratio)
ordered_stem_bind <- stem_bind[order(stem_bind[,18]),]
stem_ordered <- ordered_stem_bind[,-18]

lcl_bind <- data.frame(expr_LCL, var_ratio_lcl)
ordered_lcl_bind <- lcl_bind[order(lcl_bind[,18]),]
lcl_ordered <- ordered_lcl_bind[,-18]

#vector <-unlist(c(stem_ordered[1,]))

#boxplot(vector~indiv.fs, stem_ordered)

for(i in 1:20) {
  #plot(c(1:17), ordered_stem_bind[i,-18], col=color, cex=2, pch=20, main=rownames(stem_bind[i,]), xlab=paste("Individual"), ylab=paste("Relative Gene Expression")) 
  vectori <-  unlist(stem_ordered[i,])
  boxplot(vectori~indiv.fs, stem_ordered, names=(c("ind 1", "ind 2", "ind 3", "ind 4", "ind 5", "ind 6", "")), ylab="Gene Expression", main= rownames(stem_ordered[i,]))
}

for(i in 1:20) {
  #plot(c(1:17), ordered_stem_bind[i,-18], col=color, cex=2, pch=20, main=rownames(stem_bind[i,]), xlab=paste("Individual"), ylab=paste("Relative Gene Expression")) 
  vectori <-  unlist(lcl_ordered[i,])
  boxplot(vectori~indiv.fs, lcl_ordered, names=(c("ind 1", "ind 2", "ind 3", "ind 4", "ind 5", "ind 6", "")), ylab="Gene Expression", main= rownames(lcl_ordered[i,]))
}

# #For cycle 0 lcls
# var_ratio0 <- vector()
# for(i in 1:length(var_lcl0)) {
#   var_ratio0 <- c(var_ratio0, (mean(var_win_lcl0[i,]))/(var_lcl0[i]))  
# }
# length(var_ratio0)

## Calculate variance across stem cells and variance across LCLS
stem_var_rat <- var_ratio
lcl_var_rat <- var_ratio_lcl
t.test(stem_var_rat,lcl_var_rat)

# Density plots of variance across sample type
#Ratio
var_all <- data.frame(var=c(lcl_var_rat, stem_var_rat), type = rep(c("LCL", "Stem"), times=c(length(lcl_var_rat),length(stem_var_rat))))
ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.25,2.5)+xlab("Variance Within/Variance Between") + ggtitle("Gene Expression: Variance within vs between individuals") + theme(legend.position=c(.75,.75)) + theme(text = element_text(size=23)) 

#incl cycle 0
var_all <- data.frame(var=c(var_ratio0,lcl_var_rat, stem_var_rat), type = rep(c("LCL0", "LCL7", "Stem"), times=c(length(var_lcl0),length(lcl_var_rat),length(stem_var_rat))))
ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.25,2.5)+xlab("Variance Within/Variance Between") + ggtitle("Gene Expression: Variance within vs between individuals") + theme(legend.position=c(.75,.75)) +annotate(geom = "text", label="**p-value = 2.2e-16", x=1.9, y=.9)+theme(text = element_text(size=23)) 

#Absolute
# var_all <- data.frame(var=c(var_lcl, var_stem), type = rep(c("LCL", "Stem"), times=c(length(lcl_var_rat),length(stem_var_rat))))
var_all <- data.frame(var=c(var_stem, var_lcl), type = rep(c("Stem", "LCL"), times=c(length(stem_var_rat),length(lcl_var_rat))))
ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.5) + scale_fill_manual(values=c("#9999CC", "#CC6666"))+xlim(-.01,.2)+xlab("Variance")  + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23)) #+ annotate(geom = "text", label="LCL mean = 0.923", x=2.5, y=.9) +annotate(geom = "text", label="Stem mean = 3.002", x=2.5, y=.8) +theme(text = element_text(size=23)) 

#incl cycle 0
var_all <- data.frame(var=c(var_lcl0, var_lcl, var_stem), type = rep(c("LCL0","LCL7", "Stem"), times=c(length(var_lcl0),length(lcl_var_rat),length(stem_var_rat))))
ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23)) #+ annotate(geom = "text", label="LCL mean = 0.923", x=2.5, y=.9) +annotate(geom = "text", label="Stem mean = 3.002", x=2.5, y=.8) +theme(text = element_text(size=23)) 


## Rank genes

#Top variable genes
stems_list <- data.frame(var_ratio, names(var_ratio))
lcl_list <- data.frame(var_ratio_lcl, names(var_ratio_lcl))

#Order dataframe
ordereds<- stems_list[order(stems_list[,1]),]
orderedl<- lcl_list[order(lcl_list[,1]),]  

#Reverse order
rordereds<- stems_list[order(stems_list[,1], decreasing=TRUE),]
rorderedl<- lcl_list[order(lcl_list[,1], decreasing=TRUE),]  

#Subset top and bottom
tops <- ordereds[c(1:5177),]
topl <- orderedl[c(1:2341),]
bottom_stem <- rordereds[c(1:5000),]
bottom_lcl <- rorderedl[c(1:5000),]

#Overlap of top and bottom between cell types
overlap <- intersect(tops[,2], topl[,2])
length(overlap)
roverlap <- intersect(bottom_stem[,2], bottom_lcl[,2])
length(roverlap)

#write.table(tops, "C:/Users/a a/Documents/Lab/Variation Recovery/topstem_hgnc.txt", sep="\t")
#write.table(topl, "C:/Users/a a/Documents/Lab/Variation Recovery/toplcl_hgnc.txt", sep="\t")

## eQTL enrichment

#read in gtex eQTL data
eQTL_lcl <- read.table(c("eQTL lcl overlap results top 1500.txt"), fill=TRUE)
eQTL_stem <- read.table(c("eQTL stem overlap results top 1500.txt"),  fill=TRUE)

#V9 is gene name
unique_eQTL_lcl <- unique(eQTL_lcl[,9])
unique_eQTL_stem <- unique(eQTL_stem[,9])

#Get variances for genes with eQTLs in stems and LCLs
num_eQTL_all.l <- length(unique_eQTL_lcl)
num_eQTL_all.s <- length(unique_eQTL_stem)

num_eQTL_all <- c(num_eQTL_all.l, num_eQTL_all.s)

cortex_stems <- read.table(c("stems_cortex.tab"), fill=TRUE)
cortex_lcls <- read.table(c("lcl_cortex.tab"), fill=TRUE)

unique_cortex.s <- unique(cortex_stems[,9])
unique_cortex.l <- unique(cortex_lcls[,9])
length(unique_cortex.s)
length(unique_cortex.l)

#Pritchard's eQTLs from LCLs
eQTL <- read.table(c('final_eqtl_list.txt'))

#subset our data by genes with eQTLs or without eQTLs
eQTL_overlap.s <- intersect(eQTL[,1], tops[,2])
eQTL_overlap.l <- intersect(eQTL[,1], topl[,2])
length(eQTL_overlap.s)
length(eQTL_overlap.l)
stem_leQTL <- intersect(eQTL[,1], ordereds[,2])
stem_leQTL <- unique(stem_leQTL)
lcl_leQTL <- intersect(eQTL[,1], orderedl[,2])
lcl_leQTL <- unique(lcl_leQTL)
length(stem_leQTL)
length(lcl_leQTL)
leQTL.s <- ordereds[rownames(ordereds) %in% stem_leQTL,]
nleQTL.s <- ordereds[!(rownames(ordereds) %in% stem_leQTL),]
leQTL.l <- orderedl[c(lcl_leQTL),]
nleQTL.l <- orderedl[!(rownames(ordereds) %in% stem_leQTL),]

#mean of variance and wb ratios in genes with or without eQTLs
mean(leQTL.s[,1])
mean(nleQTL.s[,1])
mean(leQTL.l[,1])
mean(nleQTL.l[,1], na.rm=TRUE)

var_leQTL.s <- var_stem[stem_leQTL]
nvar_leQTL.s <- var_stem[!(names(var_stem) %in% rownames(leQTL.s))]
mean(var_leQTL.s)
mean(nvar_leQTL.s)
varat_leQTL.s <- var_ratio[stem_leQTL]
nvarat_leQTL.s <- var_ratio[!(names(var_ratio) %in% rownames(leQTL.s))]
mean(varat_leQTL.s)
mean(nvarat_leQTL.s)

var_leQTL.l <- var_lcl[lcl_leQTL]
nvar_leQTL.l <- var_lcl[!(names(var_lcl) %in% rownames(leQTL.l))]
mean(var_leQTL.l)
mean(nvar_leQTL.l)
varat_leQTL.l <- var_ratio_lcl[lcl_leQTL]
nvarat_leQTL.l <- var_ratio_lcl[!(names(var_ratio_lcl) %in% rownames(leQTL.l))]
mean(varat_leQTL.l)
mean(nvarat_leQTL.l)

#Make table of values
eQTLs_ls <- matrix(c(mean(var_leQTL.l), mean(nvar_leQTL.l), mean(varat_leQTL.l), mean(nvarat_leQTL.l), mean(var_leQTL.s), mean(nvar_leQTL.s), mean(varat_leQTL.s), mean(nvarat_leQTL.s)), ncol=2)
colnames(eQTLs_ls) <- c("LCLs", "iPSCs")
rownames(eQTLs_ls) <- c("Variance in eQTL genes", "Variance in non-eQTL genes", "Ratio in eQTL genes", "Ratio in non-eQTL genes")
barplot(eQTLs_ls, xlim=c(0,15), ylim=c(0,1),beside=TRUE, legend.text=TRUE)

# Order
ordered_leQTL.s<- leQTL.s[order(leQTL.s[,1]),]

ordered_leQTL.l<- leQTL.l[order(leQTL.l[,1]),]

dim(ordered_leQTL.l)

num_eQTL.sl <- length(eQTL_overlap)

eQTL_overlap_l <- intersect(eQTL[,1], topl[,2])

num_eQTL.ll <- length(eQTL_overlap_l)

num_eQTL.l <- (c(num_eQTL.ll, num_eQTL.sl))

num_eQTL <- matrix(c(num_eQTL.l, num_eQTL_all), nrow=2, byrow=TRUE)
colnames(num_eQTL) <- c("LCLs","iPSCs")
rownames(num_eQTL) <- c("LCLs", "All Tissue")
barplot(num_eQTL, main="eQTLs in highly variable genes")

subset_stem <- stems_list[eQTL_overlap,]
subset_lcl <- lcl_list[eQTL_overlap_l,]

var_all <- data.frame(var=c(subset_lcl[,1], subset_stem[,1]), type = rep(c("LCL", "Stem"), times=c(nrow(subset_lcl),nrow(subset_stem))))
ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.25,2.5)+xlab("Variance Within/Variance Between") + ggtitle("Gene Expression: Variance within vs between individuals") + theme(legend.position=c(.75,.75)) + theme(text = element_text(size=23)) 

# 
# #ordereds_8000 <- ordereds[c(1:8000),]
# 
# #Var btw/ var within
# #Calculate ratio of variance btw and within individuals for every gene
# var_ratio_bw <- vector()
# for(i in 1:length(var_stem)) {
#   var_ratio_bw <- c(var_ratio_bw, (var_stem[i])/(mean(var_win[i,])))  
# }
# length(var_ratio_bw)
# 
# ordered_stem_bw <- var_ratio_bw[order(var_ratio_bw)]
# 
# var_ratio_bwl <- vector()
# for(i in 1:length(var_lcl)) {
#   var_ratio_bwl <- c(var_ratio_bwl, (var_lcl[i])/(mean(var_win_lcl[i,])))  
# }
# length(var_ratio_bwl)
# 
# ordered_lcl_bwl <- var_ratio_bw[order(var_ratio_bw)]
# write.table(ordered_lcl_bwl, "C:/Users/a a/Documents/Lab/Variation Recovery/orderedlbw.txt", sep="\t")
# 
# 
# names_stem <- rownames(ordereds)
# top1500_names_stem <- rownames(tops)
# 
# names_lcl <- rownames(orderedl)
# top1500_names_lcl <- rownames(topl)
# 
# write.table(names_lcl, "C:/Users/a a/Documents/Lab/Variation Recovery/ensembl_lcl_names.txt", sep="\t")
# write.table(top1500_names_lcl, "C:/Users/a a/Documents/Lab/Variation Recovery/ensembl_top_1500_lcl.txt", sep="\t")

indiv = samplenames[,2]
indiv.f = as.factor(indiv)
stems = c(2,4,6,8,10,14,16,18,20,22,24,26,28,30,32,34,36)
indiv.fs = indiv.f[stems]

results <- c()
for (i in 1:nrow(abatch_stem)) {
  s = summary(lm(abatch_stem[i,]~indiv.fs));
  results<-c(results,pf(s$fstatistic[[1]],
                        s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
             s$adj.r.squared)
}

resultsM_s <- matrix(nco=2, data=results, byrow=TRUE)
mean(resultsM_s[,2])
mean(resultsM_s[,1])
boxplot(resultsM_s[,1])

#LCLs
stems = c(1,3,5,7,9,13,15,17,19,21,23,25,27,29,31,33,35)
indiv = samplenames[,2]
indiv.f = as.factor(indiv)
indiv.fs = indiv.f[stems]

results <- c()
for (i in 1:nrow(abatch_lcl)) {
  s = summary(lm(abatch_lcl[i,]~indiv.fs));
  results<-c(results,pf(s$fstatistic[[1]],
                        s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
             s$adj.r.squared)
}

resultsM <- matrix(nco=2, data=results, byrow=TRUE)
mean(resultsM[,2])
mean(resultsM[,1])
boxplot(resultsM[,1], resultsM_s[,1])

boxplot(resultsM[,2], resultsM_s[,2])

plot(resultsM[,1], resultsM_s[,1], what="density")
hist(resultsM[,2])

both <- data.frame(var=c(resultsM[,2], resultsM_s[,2]), type = rep(c("LCL", "Stem"), times=c(nrow(resultsM),nrow(resultsM_s))))
ggplot(both, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.25,2.5)+xlab("p-Value") + ggtitle("Gene expression correlation with covariate individual") + theme(legend.position=c(.75,.75)) + theme(text = element_text(size=23)) 

# #Change gene annotation
# ensembl <- useMart("ensembl")
# listDatasets(ensembl)
# ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl )
# attributes = listAttributes(ensembl)
# attributes

##Create object with name of data file:
data = c('YGilad-ST-May18-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

plot(data.lumi, what='boxplot')

## All data
data.lumi = data.lumi[,c(1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]
stems = c(1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)

### NORMALIZATION ###
# Convert raw Illumina probe intensities to expression values
# Corrects background, log2 stabiliizes variance, and quantile normalize
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)

#Take only the probes that have a detection p-value<.05 in at least one individual
hist(data.norm.all@assayData$detection)
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05)
detect.ind.all <- which(detect_quant.all > 0)

#With this threshold 28,663 probes out of 47,315 are expressed
expr_quant.all <- data.norm.all@assayData$exprs[detect.ind.all,]

###Find the column that is lumi_ID in feature data usually column 1
head(data.norm.all@featureData[[5]]) ## Should read: [1] "ILMN_1762337" "ILMN_2055271" "ILMN_1736007" "ILMN_2383229" "ILMN_1806310" "ILMN_1779670"....

###Convert expr_quant rownames to these lumi IDs
rownames(expr_quant.all)=data.norm.all@featureData[[5]][detect.ind.all]

###Subset expression by Darren's good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
probes = as.character(goodprobes$probeID) ## Convert from factor to character
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
dim(expr_quant.all.clean)
dim(expr_quant.all) #Look at this, you lose a lot of data here. Oh well.
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)

#Label your columns by sample name
samplenames = read.delim('YGilad-ST sample names switched.txt', header=TRUE)
colnames(expr_quant.all) = samplenames[(stems),1]

#Define covariates
indiv = samplenames[,2]
type = samplenames[,3]
ID = samplenames[,4]
repr_batch = samplenames[,5]
array_batch = samplenames[,6]
gender = samplenames[,7]
extr_batch = samplenames[,8]
extr_date = samplenames[,9]

#Converted categorical covariates to a factor so they are levels. This will allow you to do math on them.
indiv.f = as.factor(indiv)
type.f = as.factor(type)
ID.f=as.factor(ID)
repr_batch.f = as.factor(repr_batch)
array_batch.f = as.factor(array_batch)
gender.f=as.factor(gender)
extr_batch.f = as.factor(extr_batch)
extr_date.f = as.factor(extr_date)

#Subset stem cells and put your covariates in a list
indiv.fs = indiv.f[stems]
type.fs <- type.f[stems]
ID.fs <- ID.f[stems]
repr_batch.fs <- repr_batch.f[stems]
array_batch.fs <- array_batch.f[stems]
gender.fs <- gender.f[stems]
extr_batch.fs <- extr_batch.f[stems]
extr_date.fs <- extr_date.f[stems]
covars<-list(array_batch.fs, repr_batch.fs,indiv.fs,gender.fs,extr_date.fs,extr_batch.fs)
covars_names <- factor(c("array_batch.fs", "repr_batch.fs","indiv.fs","gender.fs","extr_date.fs","extr_batch.fs"))

#Subtract out mean. Make sure this is what you want to do.
expr_nice <- t(apply(expr_quant.all,1,function(x){x-mean(x)}))
expr_stem <- expr_quant.all

## This next section subsets your data to include only 1 probe per gene.. 

expr_genes <- expr_stem

## Convert probe IDs to gene names using Darren's file (has HGNC & ensembl)
gene_names=c()
for(i in 1:dim(expr_stem)[1]){ #what is [1] doing?
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_quant.all)[i],8])) #creates a list of gene names the same length as probe list
}
rownames(expr_genes)=gene_names
Unique_genes = unique(rownames(expr_genes))
length(Unique_genes) #this will be the length of your final data set. One data point per gene.

## This loop will give the most 3' value for multiple probes within the same gene. 
expr_gene = matrix(NA, ncol=length(stems), nrow=length(Unique_genes))
i=0
for(gene in Unique_genes){
  
  i = i+1
  
  currRows = which(Unique_genes == gene) #all probes of that gene
  if(length(currRows)>1){
    if(expr_genes[goodprobes[1],6]=="+"){ #wait so all probes for a gene are on the same strand?
      keepRow = currRows[which.max(goodprobes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  expr_gene[i,] = expr_stem[keepRow,]
  
} 

dim(expr_gene)
rownames(expr_gene) = Unique_genes
colnames(expr_gene) = colnames(expr_stem)

#cor.q <- cor(expr_gene,method="pearson", use="complete.obs")
#heatmap.2(cor.q, main="Gene Expression Correlation", key=T, revC=T, density.info="none", trace="none")

##Regress out array batch and add the intercept back in
abatch.residual.int.g = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_stem))
rownames(abatch.residual.int.g) = rownames(expr_gene)
colnames(abatch.residual.int.g) = colnames(expr_gene)
for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ array_batch.fs)
  abatch.residual.int.g[i,] = resid(model) + model$coefficients[1]
}
cor.abatch.int.g <- cor(abatch.residual.int.g,method="pearson", use="complete.obs")
#heatmap.2(abatch.residual.int.g, main="All clustered heatmap- all genes", margins=c(5,9),key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")

abatch_all <- abatch.residual.int.g

#Relationship between PCs and covariates for regressed data

npcs = 4
sum.PC <- prcomp(na.omit(abatch.residual.int.g))
results<-c()
for (f in covars) {
  for (i in 1:npcs)
  {
    s = summary(lm(sum.PC$rotation[,i]~f));
    results<-c(results,pf(s$fstatistic[[1]],
                          s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
               s$adj.r.squared)
  }
}
resultsM_gene_corrected <-matrix(nrow = length(covars), ncol = 2*npcs, data =
                                   results, byrow = TRUE)
rownames(resultsM_gene_corrected) = covars_names
colnames(resultsM_gene_corrected) = c("PC1 p value","PC1 R2","PC2 p value","PC2 R2","PC3 p value","PC3 R2","PC4 p value","PC4 R2")

resultsM_gene_corrected
PC_table_lcl <- resultsM_gene_corrected

#PCA for regressed data
sum.PC <- prcomp(na.omit(abatch.residual.int.g))
sumsum_gene_corrected <- summary(sum.PC)

op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
title.PC = "PCA of Gene Expression in lcls and iPSCs"
sum.PC <- prcomp(na.omit(abatch.residual.int.g), scale=TRUE)
sumsum <- summary(sum.PC)

#prints out plots in c(#rows, #columns)
color = indiv.fs
par(mfrow = c(1,1),oma=c(0,0,2,0)) 
plot(c(1:17),sum.PC$rotation[,1],cex=1.5,col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:17),sum.PC$rotation[,1], samplenames[,2], cex = 0.55, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC$rotation[,1], sum.PC$rotation[,i],labels=samplenames[(stems),1], cex = 0.8, pos=3)   
}  

#Making dendrograms
cor <- cor(abatch.residual.int.g, method="pearson")
dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)
color <- indiv.fs
plot(hc, main = "Cluster Dendrogram", cex.main = 1.5 , col = "#487AA1", col.main = "#45ADA8", col.lab = "lightcoral", col.axis = "#F38630", lwd = 3, lty = 1, sub = "", hang = -1, axes = FALSE)
axis(side = 2, at = seq(0.0, 1.4, .2), col = "#F38630", labels = FALSE, lwd = 2)
# add text in margin
mtext(seq(0, 1.4, .2), side = 2, at = seq(0, 1.4, .2), line = 1, col = "#A38630", las = 2)

#Within individual pearson correlation coefficient
cor_2l <- c(cor[1,6],cor[1,12],cor[6,12])
cor_5l <- c(cor[2,7],cor[2,13],cor[7,13])
cor_6l <- c(cor[3,8],cor[8,14],cor[3,14])
cor_9l <- c(cor[4,9],cor[4,15],cor[9,15])
cor_10l <- c(cor[10,16])
cor_14l <- c(cor[5,11],cor[5,17],cor[11,17])

cor_within_l <- cbind(cor_2l,cor_5l,cor_6l,cor_9l,cor_10l,cor_14l)


# #Convert to ensembl gene names
# ensembl <- useMart("ensembl")
# ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
# test <- gene_names[c(1:5)]
# #getBM(attributes=c('illumina_humanht_12_v4', 'ensembl_gene_id'), filters='illumina_humanht_12_v4', values=test, mart=ensembl)
# getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters='hgnc_symbol', values=test, mart=ensembl)

setwd("~/Lab/Variation Recovery/from minal/gene expression data for Samantha/gene expression data for Samantha")

##Create object with name of data file:
data = c('expression_data_cycle0.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

## LCLs only cycle 6
#data.lumi <- data.lumi[,c(-3,-10,-18,-23,-34,-44)]

## LCLs only cycle 0
#data.lumi <- data.lumi[,c(-4,-8,-17,-27,-37,-40)]

## My LCLs only cycle 0
data.lumi <- data.lumi[,c(5,6,10,14,19,20,21,22,24,26,28,29,31,33,34,35,42)]

## My LCLs only cycle 6
#data.lumi <- data.lumi[,c(4,10,11,16,17,18,20,23,24,29,30,31,33,37,39,41,43)]

### NORMALIZATION ###
# Convert raw Illumina probe intensities to expression values
# Corrects background, log2 stabiliizes variance, and quantile normalize
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
#data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)

#Take only the probes that have a detection p-value<.05 in at least one individual
hist(data.norm.all@assayData$detection)
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05)
detect.ind.all <- which(detect_quant.all > 0)
#With this threshold 28,663 probes out of 47,315 are expressed
expr_quant.all <- data.norm.all@assayData$exprs[detect.ind.all,]

###Find the column that is lumi_ID in feature data usually column 1
head(data.norm.all@featureData[[5]])
##[1] "ILMN_1762337" "ILMN_2055271" "ILMN_1736007" "ILMN_2383229" "ILMN_1806310" "ILMN_1779670"....
###Convert expr_quant rownames to these lumi IDs???
rownames(expr_quant.all)=data.norm.all@featureData[[5]][detect.ind.all]

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

#Label your columns by sample name
samplenames = read.delim('arraynames_cycle0_ST.txt', header=FALSE)
#LCL_names <- samplenames[c(-3,-10,-18,-23,-34,-44),1] #cycle 6
#LCL_names <- samplenames[c(-4,-8,-17,-27,-37,-40),1] #cycle 0
LCL_names <- samplenames[c(5,6,10,14,19,20,21,22,24,26,28,29,31,33,34,35,42),1] #Cycle 0 my LCLs
#LCL_names <- samplenames[c(4,10,11,16,17,18,20,23,24,29,30,31,33,37,39,41,43),3] #cycle 6 my LCLs


colnames(expr_quant.all) = LCL_names
#indiv <- samplenames[c(5,6,10,14,19,20,21,22,24,26,28,29,31,33,34,35,42),4]
indiv <- LCL_names
indiv.fs <- as.factor(indiv)
#Define covariates 
covars  <- indiv.fs

#Converted categorical covariates to a factor so they are levels

#Subtract out mean
expr_nice <- t(apply(expr_quant.all,1,function(x){x-mean(x)}))
expr_stem <- expr_nice
expr_PC <- expr_nice

#Don't subract out mean
expr_stem <- expr_quant.all

## Variance within vs among indvls
expr_LCL <- expr_stem

var_1 <- apply(expr_LCL[,c(1,16,9)],1,var)
var_2 <- apply(expr_LCL[,c(2,11,5)],1,var)
var_3 <- apply(expr_LCL[,c(13,17,8)],1,var)
var_4 <- apply(expr_LCL[,c(14,10,12)],1,var)
var_5 <- apply(expr_LCL[,c(3,7)],1,var)
var_6 <- apply(expr_LCL[,c(15,4,6)],1,var)

var_win_lcl0 <- cbind(var_1,var_2,var_3,var_4,var_5,var_6)
var_lcl0 <- apply(expr_LCL,1,var)


#Generate Heat map before regressing anything
cor <- cor(expr_stem,method="pearson", use="complete.obs")
heatmap.2(cor, main="Probe expression correlation", key=T, revC=T, density.info="none", trace="none", margins=c(5,9))

dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)
color <- indiv.fs
plot(hc, main = "Cluster Dendrogram: New LCLS", cex.main = 1.5 , col = "#487AA1", col.main = "#45ADA8", col.lab = "lightcoral", col.axis = "#F38630", lwd = 3, lty = 1, sub = "", hang = -1, axes = FALSE)

axis(side = 2, at = seq(0.0, 1.4, .2), col = "#F38630", labels = FALSE, lwd = 2)
# add text in margin
mtext(seq(0, 1.4, .2), side = 2, at = seq(0, 1.4, .2), line = 1, col = "#A38630", las = 2)

#PCA 
sum.PC <- prcomp(na.omit(expr_PC))
sumsum <- summary(sum.PC)

library(TeachingDemos)
op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
cor <- cor(abatch.residual.int, method="pearson", use="complete.obs")
title.PC = "PCA of Gene Expression in New LCLs"
sum.PC <- prcomp(na.omit(expr_stem))
sumsum <- summary(sum.PC)
#prints out plots in c(#rows, #columns)
color = indiv.fs
par(mfrow = c(1,1),oma=c(0,0,2,0)) 
plot(c(1:17),sum.PC$rotation[,1],cex=1.5,col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:17),sum.PC$rotation[,1], LCL_names, cex = 0.55, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC$rotation[,1], sum.PC$rotation[,i],labels=LCL_names, cex = 0.6, pos=3)   
}  
title.PC = "PCA of Gene Expression in LCLs"

plot(sum.PC$rotation[,2], sum.PC$rotation[,3],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 2 -", (sumsum$importance[2,2]*100),"% of variance", sep=" "), ylab=paste("PC 3 -",(sumsum$importance[2,3]*100),"% of variance", sep=" "))
text(sum.PC$rotation[,2], sum.PC$rotation[,3],labels=samplenames[(stems),2], cex = 0.8, pos=3) 
#Relationship between PCs and covariates
# this function takes your pca object, the list of covars, and the number of PCS you want to regress, the result is a matrix, N rows is the # of covariates, N columns is 2* the number of PCS
# each PC has 2 columns in the resulting matrix, column 1 is the pvalue, column 2 is the model adj R2, this is the % of variation in PC1 explained by your covariate.

npcs = 4
sum.PC <- prcomp(na.omit(expr_stem), scale.=TRUE)
results<-c()
for (i in 1:npcs)
{
  s = summary(lm(sum.PC$rotation[,i]~indiv.fs));
  results<-c(results,pf(s$fstatistic[[1]],
                        s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
             s$adj.r.squared)
}
resultsM<-matrix(nrow = 1, ncol = 2*npcs, data =
                   results, byrow = TRUE)
rownames(resultsM) = "indiv"
colnames(resultsM) = c("PC1 p value","PC1 R2","PC2 p value","PC2 R2","PC3 p value","PC3 R2","PC4 p value","PC4 R2")
resultsM
results_minal <- resultsM

#Now do the analysis by gene
## Finding the unique gene names matching probes to gene names using Darren's good probe list
gene_names=c()
for(i in 1:dim(expr_quant.all)[1]){
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_quant.all)[i],8]))
}

expr_genes <- expr_stem
rownames(expr_genes)=data.norm.all@featureData[[2]][detect.ind.all]
Unique_genes = unique(rownames(expr_genes))
length(Unique_genes)


## This loop will give the most 3' value for multiple probes within the same gene. In the end you get a simple file with all genes that are expressed with the corresponding mean intensity expression levels across its different probes.
# ncol is number of individuals
expr_gene = matrix(NA, ncol=17, nrow=length(Unique_genes))
i=0
for(gene in Unique_genes){
  
  i = i+1
  
  currRows = which(Unique_genes == gene)
  if(length(currRows)>1){
    if(expr_genes[currRows[1],6]=="+"){
      keepRow = currRows[which.max(expr_genes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(expr_genes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  expr_gene[i,] = expr_stem[keepRow,]
  
} 
dim(expr_gene)
rownames(expr_gene) = Unique_genes
colnames(expr_gene) = colnames(expr_stem)

#PCA for gene data
sum.PC <- prcomp(na.omit(expr_gene))
sumsum_gene <- summary(sum.PC)

#PCA by covariate
npcs = 4
sum.PC <- prcomp(na.omit(expr_gene))
results<-c()
for (f in covars) {
  for (i in 1:npcs)
  {
    s = summary(lm(sum.PC$rotation[,i]~f));
    results<-c(results,pf(s$fstatistic[[1]],
                          s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
               s$adj.r.squared)
  }
}
resultsM_gene<-matrix(nrow = length(covars), ncol = 2*npcs, data =
                        results, byrow = TRUE)
rownames(resultsM_gene) = covars_names
colnames(resultsM_gene) = c("PC1 p value","PC1 R2","PC2 p value","PC2 R2","PC3 p value","PC3 R2","PC4 p value","PC4 R2")
resultsM_gene


## Order genes by variance
var<-apply(expr_gene,1,var)
mean(var)
expr_plusvar <- cbind(expr_gene,var)

## Sort genes by variance
expr_sorted <- sort(expr_plusvar[,18])
expr_sorted <- expr_plusvar[order(expr_plusvar[,18], decreasing=TRUE),]
dim(expr_sorted)

## Subset most highly variable genes then remove variance column
expr_topvar <- expr_sorted[c(1:2000),]
expr_topvar <- expr_topvar[,-18]
dim(expr_topvar)

colnames(expr_topvar) = LCL_names

#cor.q <- cor(expr_topvar,method="pearson", use="complete.obs")
#heatmap.2(cor.q, main="Top 2500, stem, gene, reg array", key=T, revC=T, density.info="none", trace="none")

## Variance within vs among indvls

var_1 <- apply(expr_LCL[,c(1,6,12)],1,var)
var_2 <- apply(expr_LCL[,c(2,7,13)],1,var)
var_3 <- apply(expr_LCL[,c(3,8,14)],1,var)
var_4 <- apply(expr_LCL[,c(4,9,15)],1,var)
var_5 <- apply(expr_LCL[,c(5,10,16)],1,var)
var_6 <- apply(expr_LCL[,c(6,11,17)],1,var)

m1l <- mean(var_1)
m2l <- mean(var_2)
m3l <- mean(var_3)
m4l <- mean(var_4)
m5l <- mean(var_5)
m6l <- mean(var_6)

mean_l <- c(m1l,m2l,m3l,m4l,m5l,m6l)
mean(mean_l)
var_mean_l <- var(mean_l)
var_mean_l

var_win_lcl <- cbind(var_1,var_2,var_3,var_4,var_5,var_6)
var_lcl <- apply(expr_gene,1,var)

#Calculate ratio of variance btw and within individuals for every gene
var_ratio_new <- vector()
for(i in 1:length(var_lcl)) {
  var_ratio_new <- c(var_ratio_new, (mean(var_win_lcl[i,]))/(var_lcl[i]))  
}
length(var_ratio_new)


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
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05) #How many individuals have detection p<.05
detect.ind.all <- which(detect_quant.all > 1)
hist(detect_quant.all)

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
samplenames = read.delim('YGilad-ST sample names switched 8-28.txt', header=TRUE)
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
#expr_stem <- expr_nice
expr_stem <- expr_nice
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
  #abatch.residual.int.g[i,] = resid(model) + model$coefficients[1]
  abatch.residual.int.g[i,] = resid(model)
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
sum.PC <- prcomp(na.omit(abatch_lcl))
sumsum_gene_corrected <- summary(sum.PC)

op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
title.PC = "PCA of Gene Expression in lcls"
sum.PC <- prcomp(na.omit(abatch_lcl), scale=TRUE)
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
cor <- cor(abatch_lcl, method="pearson")
dis  <- 1-cor
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

#cor_wmeans_l <- sapply(cor_within_l, mean)

# across individual correlation coefficients
cor_plus <- cbind(cor, indiv.fs)
#indiv_plus <- c(indiv.fs, 0)

cor_plus <- rbind(cor_plus, indiv.fs)
cor_plus

cor_lcl <- cor_plus
cor_all_lcl <- c()
for (i in 1:(nrow(cor_lcl)-1)) {
  for (j in 1:(ncol(cor_lcl)-1)) {
    #print(i)
    if(cor_lcl[i,18]!= cor_lcl[18,j]) {
      cor_all_lcl <- c(cor_all_lcl, cor_lcl[i,j])
    }
    #else (i==j) {print(i)}
  }
}
cor_all_lcl
boxplot(cor_all_lcl)

# Euclidean distance within and between indvl of PC projections 1 & 2
dist12_ind <- c()
PC12_L <- cbind(sum.PC$rotation[,1], sum.PC$rotation[,2])
#Dist by individual
ind1 <- PC12_L[c(1,6,12),]
ind2 <- PC12_L[c(2,7,13),]
ind3 <- PC12_L[c(3,8,14),]
ind4 <- PC12_L[c(4,9,15),]
ind5 <- PC12_L[c(10,16),]
ind6 <- PC12_L[c(5,11,17)]

md1 <- mean(dist(ind1))
md2 <- mean(dist(ind2))
md3 <- mean(dist(ind3))
md4 <- mean(dist(ind4))
md5 <- mean(dist(ind5))
md6 <- mean(dist(ind6))
mdL <- c(md1,md2,md3,md4,md5,md6)
mdL_mean <- mean(mdL)

#Dist all

Edist_lcl <- c(mdL_mean, mean(dist(PC12_L)))
names(Edist_lcl)<- c("ED12 within LCL", "ED12 between LCL")
Edist_lcl

t.test(mdL, dist(PC12_L))

# STEMS

setwd("~/Lab/Variation Recovery")

##Create object with name of data file:
data = c('YGilad-ST-May18-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

plot(data.lumi, what='boxplot')

##LCL only
data.lumi = data.lumi[,c(2,4,6,8,10,14,16,18,20,22,24,26,28,30,32,34,36)] ## CHANGE TO STEM
stems = c(2,4,6,8,10,14,16,18,20,22,24,26,28,30,32,34,36)

### NORMALIZATION ###
# Convert raw Illumina probe intensities to expression values
# Corrects background, log2 stabiliizes variance, and quantile normalize
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)

#Take only the probes that have a detection p-value<.05 in at least one individual
hist(data.norm.all@assayData$detection)
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05)
detect.ind.all <- which(detect_quant.all > 1)

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
samplenames = read.delim('YGilad-ST sample names switched 8-28.txt', header=TRUE)
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
RMSD = samplenames[,10]

#Converted categorical covariates to a factor so they are levels. This will allow you to do math on them.
indiv.f = as.factor(indiv)
type.f = as.factor(type)
ID.f=as.factor(ID)
repr_batch.f = as.factor(repr_batch)
array_batch.f = as.factor(array_batch)
gender.f=as.factor(gender)
extr_batch.f = as.factor(extr_batch)
extr_date.f = as.factor(extr_date)
rmsd.f <- as.factor(RMSD)

#Subset stem cells and put your covariates in a list
indiv.fs = indiv.f[stems]
type.fs <- type.f[stems]
ID.fs <- ID.f[stems]
repr_batch.fs <- repr_batch.f[stems]
array_batch.fs <- array_batch.f[stems]
gender.fs <- gender.f[stems]
extr_batch.fs <- extr_batch.f[stems]
extr_date.fs <- extr_date.f[stems]
rmsd.fs <- rmsd.f[stems]
covars<-list(array_batch.fs,indiv.fs,gender.fs,extr_date.fs,extr_batch.fs) #leave out repr batch for lcls, will bug later
covars_names <- factor(c("array_batch.fs","indiv.fs","gender.fs","extr_date.fs","extr_batch.fs"))

#Subtract out mean. Make sure this is what you want to do.
expr_nice <- t(apply(expr_quant.all,1,function(x){x-mean(x)}))
#expr_stem <- expr_nice
expr_stem <- expr_nice
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
  #abatch.residual.int.g[i,] = resid(model) + model$coefficients[1]
  abatch.residual.int.g[i,] = resid(model)
}
cor.abatch.int.g <- cor(abatch.residual.int.g,method="pearson", use="complete.obs")
heatmap.2(cor.abatch.int.g, Rowv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))),
          Colv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))), margins=c(5,9),key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")

abatch_stem <- abatch.residual.int.g ## CHANGE TO STEM

abatch_small_stem <- abatch_stem[1:14828,]
dim(abatch_small_stem)

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
PC_table_stem <- resultsM_gene_corrected ##CHANGE TO STEM

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
cor <- cor(abatch_stem, method="pearson")
#cor <- cor(abatch_small_stem)
dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)
color <- indiv.fs
plot(hc, main = "Cluster Dendrogram: iPSCs", cex.main = 1.5 , col = "#487AA1", col.main = "#45ADA8", col.lab = "lightcoral", col.axis = "#F38630", lwd = 3, lty = 1, sub = "", hang = -1, axes = FALSE)
axis(side = 2, at = seq(0.0, 1.4, .2), col = "#F38630", labels = FALSE, lwd = 2)
# add text in margin
mtext(seq(0, 1.4, .2), side = 2, at = seq(0, 1.4, .2), line = 1, col = "#A38630", las = 2)

hclust(reorder(distance))

#Within individual pearson correlation coefficient
## CHANGE TO STEM
cor_2s <- c(cor[1,6],cor[1,12],cor[6,12])
cor_5s <- c(cor[2,7],cor[2,13],cor[7,13])
cor_6s <- c(cor[4,8],cor[8,14],cor[4,14]) ## CHANGE TO STEM
cor_9s <- c(cor[3,9],cor[3,15],cor[9,15]) ## CHANGE TO STEM
cor_10s <- c(cor[10,16])
cor_14s <- c(cor[5,11],cor[5,17],cor[11,17])

cor_within_s <- cbind(cor_2s,cor_5s,cor_6s,cor_9s,cor_10s,cor_14s)

#cor_wmeans_s <- sapply(cor_within_s, mean)

# across individual correlation coefficients
cor_plus <- cbind(cor, indiv.fs)
#indiv_plus <- c(indiv.fs, 0)

cor_plus <- rbind(cor_plus, indiv.fs)
cor_plus

cor_stem <- cor_plus
cor_all_stem <- c()
for (i in 1:(nrow(cor_stem)-1)) {
  for (j in 1:(ncol(cor_stem)-1)) {
    #print(i)
    if(cor_stem[i,18]!= cor_stem[18,j]) {
      cor_all_stem <- c(cor_all_stem, cor_stem[i,j])
    }
    #else (i==j) {print(i)}
  }
}
cor_all_stem
boxplot(cor_all_stem)

# Euclidean distance within and between indvl of PC projections 1 & 2
dist12_ind <- c()
PC12 <- cbind(sum.PC$rotation[,1], sum.PC$rotation[,2])
#Dist by individual
ind1 <- PC12[c(1,6,12),]
ind2 <- PC12[c(2,7,13),]
ind3 <- PC12[c(4,8,14),] ## CHANGE TO STEM
ind4 <- PC12[c(3,9,15),] ## CHANGE TO STEM
ind5 <- PC12[c(10,16),]
ind6 <- PC12[c(5,11,17)]

md1 <- mean(dist(ind1))
md2 <- mean(dist(ind2))
md3 <- mean(dist(ind3))
md4 <- mean(dist(ind4))
md5 <- mean(dist(ind5))
md6 <- mean(dist(ind6))
mdS <- c(md1,md2,md3,md4,md5,md6) ## CHANGE TO STEM
mdS_mean <- mean(mdS) ## CHANGE TO STEM

#Dist all
dist12S <- dist(cbind(sum.PC$rotation[,1], sum.PC$rotation[,2])) ## CHANGE TO STEM
mean(dist12S) ## CHANGE TO STEM
# Euclidean distance within and between indvl of PC projections 1 & 2
dist12_ind <- c()
PC12 <- cbind(sum.PC$rotation[,1], sum.PC$rotation[,2])
#Dist by individual
ind1 <- PC12[c(1,6,12),]
ind2 <- PC12[c(2,7,13),]
ind3 <- PC12[c(4,8,14),] ## CHANGE TO STEM
ind4 <- PC12[c(3,9,15),] ## CHANGE TO STEM
ind5 <- PC12[c(10,16),]
ind6 <- PC12[c(5,11,17)]

md1 <- mean(dist(ind1))
md2 <- mean(dist(ind2))
md3 <- mean(dist(ind3))
md4 <- mean(dist(ind4))
md5 <- mean(dist(ind5))
md6 <- mean(dist(ind6))
mdS <- c(md1,md2,md3,md4,md5,md6) ## CHANGE TO STEM
mdS_mean <- mean(mdS) ## CHANGE TO STEM

#Dist all
ED_t_stem <- t.test(mdS, dist(PC12))
Edist_stem <- c(mdS_mean, mean(dist(PC12)), ED_t_stem$p.value)
names(Edist_stem)<- c("ED12 within", "ED12 between", "p-value")
Edist_stem

ED_t_lcl <- t.test(mdL, dist(PC12_L))
Edist_lcl <- c(mdL_mean, mean(dist(PC12_L)), ED_t_lcl$p.value)
names(Edist_lcl)<- c("ED12 within", "ED12 between", "p-value")
Edist_lcl

Edist <- cbind(Edist_stem, Edist_lcl)
colnames(Edist) <- c("iPSCs", "LCLs")
Edist_t <- t(Edist)
Edist_t

###########################################################################################
# All together now
cor_wmeanl <- sapply(cor_within_l, mean)
cor_wmeans <- sapply(cor_within_s, mean) ## CHANGE TO STEM
cor_both <- cbind(cor_within_s, cor_within_l) ## CHANGE TO STEM
cor_bmeans <- cbind(cor_wmeans, cor_wmeansl)
boxplot(cor_bmeans)
boxplot(cor_both)
boxplot(cor_bmeans, main="Within Individual Pearson correlation coefficients for Stem cells and LCLs")
cor_all_both <- cbind(cor_all_lcl, cor_all_stem)
boxplot(cor_all_both)

#max <- max(nrow(cor_bmeans), nrow(cor_all_both))

#length(cor_within_s) <- max
#length(cor_within_l) <- max

cor_total <- cbind(cor_wmeans, cor_wmeanl, cor_all_stem, cor_all_lcl)
cor_total_vec <- c(cor_wmeans, cor_wmeanl, cor_all_stem, cor_all_lcl)

boxplot(cor_total)

type <- c(rep("stem", times=length(cor_wmeans)), rep("lcl", times=length(cor_wmeanl)), rep("stem", times=length(cor_all_stem)), rep("lcl", times=length(cor_all_lcl)))
group <- c(rep(as.factor("Within Individuals"), times=(length(cor_wmeans)+length(cor_wmeanl))), rep(as.factor("Between Individuals"), times=(length(cor_all_stem)+length(cor_all_lcl))))
groupn <- c(rep("1", times=(length(cor_wmeans)+length(cor_wmeanl))), rep("2", times=length(cor_all_stem)+length(cor_all_lcl)))

cor_df <- data.frame(cor_total_vec, type, group, groupn)
#cor_df$group = factor(cor_df$group,levels=c("Within Individuals", "Between Individuals"))
t_cor_w <- t.test(cor_within_s, cor_within_l)
pv_w <- signif(t_cor_w$p.value,3)
pv_w
t_cor_b <- t.test(cor_all_lcl, cor_all_stem)
pv_b <- signif(t_cor_b$p.value, 3)
pv_b

ggplot(cor_df, aes(x=groupn, y=cor_total_vec)) + geom_boxplot(aes(fill=type), xlab=FALSE) + theme(text = element_text(size=18)) + annotate(geom = "text", label=paste("**p =", pv_w), x=1, y=0.1) + annotate(geom = "text", label=paste("**p =", pv_b), x=2, y=-.45) + scale_x_discrete(labels=c("Within Individuals", "Between Individuals")) + theme(axis.title.x=element_blank(), panel.background=element_rect(fill='white')) + ylab("Pearson Correlation") + theme(axis.title.y = element_text(vjust=1.0)) + #stat_boxplot(geom ='errorbar', aes(x=group))  
mean(cor_wmeans)
mean(cor_wmeansl)
median(cor_wmeans)
median(cor_wmeansl)

# Euclidean distance within and between indvl of PC projections 1 & 2
t.test(dist12S, mdS)
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
  var_ratio_lcl <- c(var_ratio_lcl, (var_lcl[i]/(mean(var_win_lcl[i,]))))  
}
length(var_ratio_lcl)
# for(i in 1:length(var_lcl)) {
#   var_ratio_lcl <- c(var_ratio_lcl, ((mean(var_win_lcl[i,]))/var_lcl[i]))  
# }
# length(var_ratio_lcl)

#Do same for stem cells REMEMBER 9 AND 6 ARE SWITCHED SO THESE LABELS AREN'T THE SAME AS LCLS
expr_stem <- abatch_stem #if you've stored it this way

var_1s <- apply(expr_stem[,c(1,6,12)],1,var)
var_2s <- apply(expr_stem[,c(2,7,13)],1,var)
var_3s <- apply(expr_stem[,c(4,8,14)],1,var)
var_4s <- apply(expr_stem[,c(3,9,15)],1,var)
var_5s <- apply(expr_stem[,c(10,16)],1,var)
var_6s <- apply(expr_stem[,c(5,11,17)],1,var)

var_win_stem <- cbind(var_1s,var_2s,var_3s,var_4s,var_5s,var_6s)
var_stem <- apply(expr_stem,1,var) #calculate variance for each gene across all individuals
#rownames(var_stem) <- rownames(expr_stem)

#Calculate ratio of variance btw and within individuals for every gene
var_ratio_stem <- vector()
for(i in 1:length(var_stem)) {
  var_ratio_stem <- c(var_ratio_stem, (var_stem[i]/(mean(var_win_stem[i,]))))
}
length(var_ratio_stem)

mean(var_stem)
mean(var_lcl)

# for(i in 1:length(var_stem)) {
#   var_ratio_stem <- c(var_ratio_stem, ((mean(var_win_stem[i,]))/var_stem[i]))
# }
# length(var_ratio_stem)

## Calculate variance across stem cells and variance across LCLS
stem_var_rat <- var_ratio_stem
mean(stem_var_rat)
lcl_var_rat <- var_ratio_lcl
mean(lcl_var_rat)
t_var_rat <- t.test(stem_var_rat,lcl_var_rat)
pv_var_rat <- as.character(signif(t_var_rat$p.value, 5)) #round() will give you 0
pv_var_rat

t_var <- t.test(var_stem, var_lcl)
t_var
pv_var <- as.character(signif(t_var$p.value, 5))

## Rank genes

#Top variable genes
stems_list <- data.frame(var_ratio_stem, names(var_ratio_stem))
lcl_list <- data.frame(var_ratio_lcl, names(var_ratio_lcl))

# For Courtney
stems_list_both <- data.frame(var_stem, var_ratio_stem)
write.table(stems_list_both, "C:/Users/a a/Documents/Lab/Variation Recovery/Variance iPSCs.txt", sep="\t")

#Order dataframe
ordereds<- stems_list[order(stems_list[,1]),]
orderedl<- lcl_list[order(lcl_list[,1]),]  

#Reverse order
rordereds<- stems_list[order(stems_list[,1], decreasing=TRUE),]
rorderedl<- lcl_list[order(lcl_list[,1], decreasing=TRUE),]  

#Subset top and bottom
tops <- ordereds[c(1:2341),]
topl <- orderedl[c(1:2341),]
bottom_stem <- rordereds[c(1:5177),]
bottom_lcl <- rorderedl[c(1:2341),]
bottom_stemnames <- row.names(bottom_stem)
bottom_lclnames <- row.names(bottom_lcl)

#Overlap of top and bottom between cell types
overlap <- intersect(tops[,2], topl[,2])
length(overlap)
roverlap <- intersect(bottom_stem[,2], bottom_lcl[,2])
length(roverlap)

#write(bottom_lclnames, "C:/Users/a a/Documents/Lab/Variation Recovery/high var lcl.txt", sep="\t")
#write(bottom_stemnames, "C:/Users/a a/Documents/Lab/Variation Recovery/high var stem.txt", sep="\t")
write(roverlap, "C:/Users/a a/Documents/Lab/Variation Recovery/roverlap.txt", sep="\t")

# Order by variance ratio. Plot expression by individual CHANGE INDIVIDUAL IDENTIFIERS!!!!

color=indiv.fs

stem_bind <- data.frame(expr_stem, var_ratio_stem)
ordered_stem_bind <- stem_bind[order(stem_bind[,18], decreasing=TRUE),]
stem_ordered <- ordered_stem_bind[,-18]

lcl_bind <- data.frame(expr_LCL, var_ratio_lcl)
ordered_lcl_bind <- lcl_bind[order(lcl_bind[,18], decreasing=TRUE),]
lcl_ordered <- ordered_lcl_bind[,-18]

for(i in 1:2) {
  #plot(c(1:17), ordered_stem_bind[i,-18], col=color, cex=2, pch=20, main=rownames(stem_bind[i,]), xlab=paste("Individual"), ylab=paste("Relative Gene Expression")) 
  vectori <-  unlist(stem_ordered[i,])
  boxplot(vectori~indiv.fs, stem_ordered, names=(c("ind 1", "ind 2", "ind 3", "ind 4", "ind 5", "ind 6", "")), ylab="Gene Expression", main= rownames(stem_ordered[i,]))
}

for(i in 1:2) {
  #plot(c(1:17), ordered_stem_bind[i,-18], col=color, cex=2, pch=20, main=rownames(stem_bind[i,]), xlab=paste("Individual"), ylab=paste("Relative Gene Expression")) 
  vectori <-  unlist(lcl_ordered[i,])
  boxplot(vectori~indiv.fs, lcl_ordered, names=(c("ind 1", "ind 2", "ind 3", "ind 4", "ind 5", "ind 6", "")), ylab="Gene Expression", main= rownames(lcl_ordered[i,]))
}

#pvalues for each gene association with factor individual

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

plot(stem_var_rat, pvals_stems, xlim=c(0,5))
plot(lcl_var_rat, pvals_lcls)

mean(pvals_stems)
mean(pvals_lcls)

#count significant genes
sig_stems <- length(pvals_stems[pvals_stems<.05])
sig_lcls <- length(pvals_lcls[pvals_lcls<.05])

cutoff_stems <- ordereds[(nrow(ordereds)-sig_stems),]
cutoff_lcls <- orderedl[(nrow(ordereds)-sig_lcls),]
cutoff_stems
cutoff_lcls
cutoff_avg <- mean(c(cutoff_stems[,1], cutoff_lcls[,1]))
cutoff_avg

# Density plots of variance across sample type
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = gg_color_hue(2)

#Ratio
#pretty <- c("palevioletred2", "royalblue2")
var_all <- data.frame(var=c(stem_var_rat, lcl_var_rat), type = rep(c("iPSC", "LCL"), times=c(length(stem_var_rat),length(lcl_var_rat))), Cell_type = rep(c("2", "1"), times=c(length(stem_var_rat),length(lcl_var_rat))))
#var_all <- data.frame(var=c(lcl_var_rat, stem_var_rat), type = rep(c("LCL", "iPSC"), times=c(length(lcl_var_rat),length(stem_var_rat))))
ggplot(var_all, aes(x=var, fill=Cell_type)) + geom_density(alpha=0.5) +xlim(-.25,8.0)+xlab("Variance Between/Variance Within") + geom_vline(xintercept=cutoff_avg, linetype="dotted") + theme(legend.position=c(.75,.75), panel.background=element_rect(fill='white')) + theme(text = element_text(size=18)) +annotate(geom = "text", label=paste("p-value = ", pv_var_rat), x=6, y=.9) +ylab("Number of Genes") + scale_fill_manual(values=cols, labels=c("LCLs", "iPSCs"))
# 
# #incl cycle 0
#var_all <- data.frame(var=c(var_ratio0,lcl_var_rat, stem_var_rat), type = rep(c("LCL0", "LCL7", "Stem"), times=c(length(var_lcl0),length(lcl_var_rat),length(stem_var_rat))))
#ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.25,2.5)+xlab("Variance Within/Variance Between") + ggtitle("Gene Expression: Variance within vs between individuals") + theme(legend.position=c(.75,.75)) +annotate(geom = "text", label="**p-value = 2.2e-16", x=1.9, y=.9)+theme(text = element_text(size=23)) 

#Absolute

#pretty <- rev(pretty)
# var_all <- data.frame(var=c(var_lcl, var_stem), type = rep(c("LCL", "Stem"), times=c(length(lcl_var_rat),length(stem_var_rat))))
var_all <- data.frame(var=c(var_stem, var_lcl), type = rep(c("iPSC", "LCL"), times=c(length(stem_var_rat),length(lcl_var_rat))), Cell_type = rep(c("1", "2"), times=c(length(stem_var_rat),length(lcl_var_rat))))
ggplot(var_all, aes(x=var, fill=Cell_type)) + annotate(geom = "text", label=paste("p-value = ", pv_var), x=.075, y=65) +geom_density(alpha=0.5) +xlim(-.01,.1)+xlab("Variance")  + theme(legend.position=c(.75,.75), panel.background=element_rect(fill='white')) +theme(text = element_text(size=18))+ylab("Number of Genes") + scale_fill_manual(values=rev(cols), labels=c("iPSCs", "LCLs"))

#+ annotate(geom = "text", label="LCL mean = 0.923", x=2.5, y=.9) +annotate(geom = "text", label="Stem mean = 3.002", x=2.5, y=.8) +theme(text = element_text(size=23)) 
# 
# #incl cycle 0
# var_all <- data.frame(var=c(var_lcl0, var_lcl, var_stem), type = rep(c("LCL0","LCL7", "Stem"), times=c(length(var_lcl0),length(lcl_var_rat),length(stem_var_rat))))
# ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23)) #+ annotate(geom = "text", label="LCL mean = 0.923", x=2.5, y=.9) +annotate(geom = "text", label="Stem mean = 3.002", x=2.5, y=.8) +theme(text = element_text(size=23)) 


## eQTL enrichment

#read in gtex eQTL data
#eQTL_lcl <- read.table(c("eQTL lcl overlap results top 1500.txt"), fill=TRUE)
#eQTL_stem <- read.table(c("eQTL stem overlap results top 1500.txt"),  fill=TRUE)

eQTL_LCL_L <- read.table(c("8-28 sig LCL-LCL eQTLs.tab"), fill=TRUE)
unique_LL <- unique(eQTL_LCL_L[,9])
length(unique_LL)

eQTL_LCL_S <- read.table(c("8-28 sig LCL-Stem eQTLs.tab"), fill=TRUE)
unique_LS <- unique(eQTL_LCL_S[,9])
length(unique_LS)

eQTL_liv_s <- read.table(c("8-28 sig liver stem eQTLs.tab"), fill=TRUE)
unique_liv_s <- unique(eQTL_liv_s[,9])
length(unique_liv_s)

eQTL_liv_l <- read.table(c("8-28 sig liver  lcl eQTLs.tab"), fill=TRUE)
unique_liv_l <- unique(eQTL_liv_l[,9])
length(unique_liv_l)

eQTL_cereb_l <- read.table(c("8-28 sig cereb lcl eQTLs.tab"), fill=TRUE)
unique_cereb_l <- unique(eQTL_cereb_l[,9])
length(unique_cereb_l)

eQTL_cereb_s <- read.table(c("8-28 sig cereb stem eQTLs.tab"), fill=TRUE)
unique_cereb_s <- unique(eQTL_cereb_s[,9])
length(unique_cereb_s)

eQTL_cortex_s <- read.table(c("8-28 sig cortex stem eQTLs.tab"), fill=TRUE)
unique_cortex_s <- unique(eQTL_cortex_s[,9])
length(unique_cortex_s)

eQTL_cortex_l <- read.table(c("8-28 sig cortex lcl eQTLs.tab"), fill=TRUE)
unique_cortex_l <- unique(eQTL_cortex_l[,9])
length(unique_cortex_l)

eQTL_tcortex_l <- read.table(c("8-28 sig tcortex lcl eQTLs.tab"), fill=TRUE)
unique_tcortex_l <- unique(eQTL_tcortex_l[,9])
length(unique_tcortex_l)

eQTL_tcortex_s <- read.table(c("8-28 sig tcortex stem eQTLs.tab"), fill=TRUE)
unique_tcortex_s <- unique(eQTL_tcortex_s[,9])
length(unique_tcortex_s)

eQTL_pons_l <- read.table(c("8-28 sig pons lcl eQTLs.tab"), fill=TRUE)
unique_pons_l <- unique(eQTL_pons_l[,9])
length(unique_pons_l)

eQTL_pons_s <- read.table(c("8-28 sig pons stem eQTLs.tab"), fill=TRUE)
unique_pons_s <- unique(eQTL_pons_s[,9])
length(unique_pons_s)

eQTL_lcl2_s <- read.table(c("8-28 sig lcl stem 2 eQTLs.tab"), fill=TRUE)
unique_lcl2_s <- unique(eQTL_lcl2_s[,9])
length(unique_lcl2_s)

eQTL_lcl2_l <- read.table(c("8-28 sig lcl lcl 2 eQTLs.tab"), fill=TRUE)
unique_lcl2_l <- unique(eQTL_lcl2_l[,9])
length(unique_lcl2_l)

#eQTLs in overlap
eQTL_roverlap_lcl <- read.table(c("8-31 overlap LCL eQTLs.tab"), fill=TRUE)
unique_overlap_lcl <- unique(eQTL_roverlap_lcl[,9])
length(unique_overlap_lcl)

eQTL_roverlap_liv <- read.table(c("8-31 overlap Liver eQTLs.tab"), fill=TRUE)
unique_overlap_liv <- unique(eQTL_roverlap_liv[,9])
length(unique_overlap_liv)

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
varat_leQTL.s <- var_ratio_stem[stem_leQTL]
nvarat_leQTL.s <- var_ratio_stem[!(names(var_ratio_stem) %in% rownames(leQTL.s))]
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

#Variance and Variance ratio in eQTLs across all tissues

eQTL_alltis <- read.table(c("eQTLs all tissues.tab"), fill=TRUE)
eQTL_alltis <- unique(eQTL_alltis[,9])

length(eQTL_alltis)

alltis_stems <- ordereds[rownames(ordereds) %in% eQTL_alltis,]
dim(alltis_stems)
mean(alltis_stems[,1])
#intersect(alltis_stems[,2], eQTL_alltis)

alltis_lcls <- orderedl[rownames(orderedl) %in% eQTL_alltis,]
dim(alltis_lcls)
mean(alltis_lcls[,1])

all_stemvar <- var_stem[names(var_stem) %in% eQTL_alltis]
length(all_stemvar)
mean(all_stemvar)

all_lclvar <- var_lcl[names(var_lcl) %in% eQTL_alltis]
length(all_lclvar)
mean(all_lclvar)

## LCL eQTLs only
eQTL_lcls <- read.table(c("eQTLs lcls.tab"), fill=TRUE)
eQTL_lcls <- unique(eQTL_lcls[,9])

length(eQTL_lcls)

lcls_stems <- ordereds[rownames(ordereds) %in% eQTL_lcls,]
dim(lcls_stems)
mean_elcl_stemsrat <-lcls_stems[,1]
#intersect(lcls_stems[,2], eQTL_lcls)

lcls_lcls <- orderedl[rownames(orderedl) %in% eQTL_lcls,]
dim(lcls_lcls)
mean_elcl_lclrat <- lcls_lcls[,1]

lcl_stemvar <- var_stem[names(var_stem) %in% eQTL_lcls]
length(lcl_stemvar)
mean_elcl_stemvar <- lcl_stemvar

lcl_lclvar <- var_lcl[names(var_lcl) %in% eQTL_lcls]
length(lcl_lclvar)
mean_elcl_lclvar <- (lcl_lclvar)

bothvar_lcl_eQTLs <- c(mean_elcl_stemsrat, mean_elcl_lclrat, mean_elcl_stemvar, mean_elcl_lclvar)
var_lcl_eQTLs <- c(mean_elcl_stemvar, mean_elcl_lclvar, mean(var_stem), mean(var_lcl))
names_var <- c(rep("iPSC", times=length(lcls_stems)),  rep("LCL", times=length(lcls_lcls)), rep("iPSC", times=length(var_stem)), rep("LCL", times=length(var_lcl)))

type_var <- c("eQTL var", "eQTL var", "Mean var", "Mean var")
type_var <- c((rep("eQTL var",  times=(length(lcls_stems)+length(lcls_lcls)))), (rep("mean var", times=(length(var_stem)+length(var_lcl))))) 

df_var <- data.frame(var_lcl_eQTLs, names_var, type_var)
df_var

df_var_se <- summarySE(df_var, measurevar="len", groupvars=c("supp","dose"))

ggplot(df_var, aes(type_var,var_lcl_eQTLs), fill=names_var) +geom_bar(aes(fill=names_var), stat="identity", position=position_dodge())

# For GSEA
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

# Calculate proportion variation explained by individual of origin

indiv = samplenames[,2]
indiv.f = as.factor(indiv)

#Stems
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
exp_s <- mean(resultsM_s[,2])
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
exp_l <- mean(resultsM[,2])
mean(resultsM[,1])
boxplot(resultsM[,1], resultsM_s[,1])

boxplot(resultsM[,2], resultsM_s[,2])
barplot(resultsM[,2], resultsM_s[,2])

plot(resultsM[,1], resultsM_s[,1], what="density")
hist(resultsM[,2])

exp_PC1s <- PC_table_stem[2,2]
exp_PC1l <- PC_table_lcl[2,2]
explained_avg <- c(exp_l, exp_s)
df_avg <- data.frame(explained_avg)
explained_PC <- c(exp_PC1s, exp_PC1l)
explained <- c(exp_l, exp_s, exp_PC1l, exp_PC1s)

#type <- c("iPSC","LCL", "iPSC", "LCL")
type <- c("LCL", "iPSC","LCL", "iPSC")
cat <- c("mean explained", "mean explained","PC1", "PC1")
means_all <- cbind(explained, type)
means_all
df_ma <- data.frame(explained)
df_ma <- cbind (df_ma, type2, cat)
row.names(df_ma) <- cat2
df_ma

df_ma_me <- df_ma[c(1,2),]
df_ma_me

ggplot(df_ma_me, fill=type) + geom_bar(aes(row.names(df_ma_me), explained), stat="identity")
ggplot(df_ma, aes(x=cat, y=explained), fill=type) + geom_bar(aes(fill=type), stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(color="black"), axis.title.x=element_blank(), text=element_text(size=18)) + theme(panel.background=element_rect(fill='white'))#+ geom_errorbar(aes(ymin=val-se, ymax=val+se))
ggplot(df_avg, aes(y=explained_avg), fill=type) + geom_bar(aes(fill=type), stat="identity", position=position_dodge()) #+ geom_errorbar(aes(ymin=val-se, ymax=val+se))


ggplot(df_ma, aes(x=cat, y=explained), fill=type) + geom_bar(aes(fill=type), stat="identity")

both <- data.frame(var=c(resultsM[,2], resultsM_s[,2]), type = rep(c("LCL", "Stem"), times=c(nrow(resultsM),nrow(resultsM_s))))
ggplot(both, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.25,2.5)+xlab("p-Value") + ggtitle("Gene expression correlation with covariate individual") + theme(legend.position=c(.75,.75)) + theme(text = element_text(size=23)) 


#Add PCs?
PC_table_stem
PC_table_lcl


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
PC_table_stem <- resultsM_gene_corrected ##CHANGE TO STEM

library(plyr)
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## LCL eQTLs only
eQTL_lcls <- read.table(c("eQTLs lcls.tab"), fill=TRUE)
eQTL_lcls <- unique(eQTL_lcls[,9])

length(eQTL_lcls)

lcls_stems <- ordereds[rownames(ordereds) %in% eQTL_lcls,]
dim(lcls_stems)
mean_elcl_stemsrat <-lcls_stems[,1]
#intersect(lcls_stems[,2], eQTL_lcls)

lcls_lcls <- orderedl[rownames(orderedl) %in% eQTL_lcls,]
dim(lcls_lcls)
mean_elcl_lclrat <- lcls_lcls[,1]

lcl_stemvar <- var_stem[names(var_stem) %in% eQTL_lcls]
length(lcl_stemvar)
elcl_stemvar <- lcl_stemvar

lcl_lclvar <- var_lcl[names(var_lcl) %in% eQTL_lcls]
length(lcl_lclvar)
elcl_lclvar <- lcl_lclvar

bothvar_lcl_eQTLs <- c(mean_elcl_stemsrat, mean_elcl_lclrat, mean_elcl_stemvar, mean_elcl_lclvar)
var_lcl_eQTLs <- c(elcl_stemvar, elcl_lclvar, var_stem, var_lcl)
names_var <- c(rep("iPSC", times=length(lcl_stemvar)),  rep("LCL", times=length(lcl_lclvar)), rep("iPSC", times=length(var_stem)), rep("LCL", times=length(var_lcl)))

type_var <- c("eQTL var", "eQTL var", "Mean var", "Mean var")
type_var <- c((rep("eQTL var",  times=(length(lcl_stemvar)+length(lcl_lclvar)))), (rep("mean var", times=(length(var_stem)+length(var_lcl))))) 

df_var <- data.frame(var_lcl_eQTLs, names_var, type_var, fill=TRUE)
df_var

df_var_se <- summarySE(df_var, measurevar="var_lcl_eQTLs", groupvars=c("type_var","names_var"))
df_var_se
dodge <- position_dodge(width=0.9)
ggplot(df_var_se, aes(type_var,var_lcl_eQTLs, fill=names_var)) +geom_bar(stat="identity", position=dodge) + geom_errorbar(aes(ymin=var_lcl_eQTLs-se, ymax=var_lcl_eQTLs+se), width=0.2, position=dodge)


t.test(elcl_stemvar, elcl_lclvar)
t.test(elcl_stemvar, var_stem)

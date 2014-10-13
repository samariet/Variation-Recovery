# 
# # updateR()
# # Instead of updating packageess, reinstall.
#  source("http://bioconductor.org/biocLite.R")
#  biocLite()  
# # 
#  update.packages(ask = FALSE, dependencies = c('Suggests'))

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
library(sva)
library(qpcR)
library(limma)
library(ClassDiscovery)
library(stats)

############################################################################################
#Input Data and low level analysis: probe filtering, normalization.

setwd("~/Lab/Variation Recovery")

##Create object with name of data file:
data = c('YGilad-ST-May18-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes. Mapping arguments are Null because Probe list is annotaed:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

plot(data.lumi, what='boxplot')
plot(data.lumi, what='density')

#Keep only stem cells and iPSCs.(No hearts)
data.lumi <- data.lumi[,c(1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]

#Define covariates for both cell types, independently and together. 
both = c(1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)
stems = c(2,4,6,8,10,14,16,18,20,22,24,26,28,30,32,34,36)
lcls = c(1,3,5,7,9,13,15,17,19,21,23,25,27,29,31,33,35)


#Covariates indexes from file. 
samplenames = read.delim('YGilad-ST sample names switched 8-28.txt', header=TRUE)

#First both togethre
name =samplenames[,1]
indiv = samplenames[,2]
type = samplenames[,3]
ID = samplenames[,4]
repr_batch = samplenames[,5]
array_batch = samplenames[,6]
gender = samplenames[,7]
extr_batch = samplenames[,8]
extr_date = samplenames[,9]
RMSD = samplenames[,12]

#Converted categorical covariates to a factor so they are levels. This will allow you to do math on them.
name.f = as.factor(name)
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
name.fb =name.f[both]
indiv.fb = indiv.f[both]
type.fb <- type.f[both]
ID.fb <- ID.f[both]
repr_batch.fb <- repr_batch.f[both]
array_batch.fb <- array_batch.f[both]
gender.fb <- gender.f[both]
extr_batch.fb <- extr_batch.f[both]
extr_date.fb <- extr_date.f[both]
rmsd.fb <- rmsd.f[both]
covars<-list(array_batch.fb,indiv.fb,gender.fb,extr_date.fb,extr_batch.fb) #leave out repr batch for lcls, will bug later
covars_names <- factor(c("array_batch.fb","indiv.fb","gender.fb","extr_date.fb","extr_batch.fb"))

#Subset stem cells and put your covariates in a list
name.fs =name.f[stems]
indiv.fs = indiv.f[stems]
type.fs <- type.f[stems]
ID.fs <- ID.f[stems]
repr_batch.fs <- repr_batch.f[stems]
array_batch.fs <- array_batch.f[stems]
gender.fs <- gender.f[stems]
extr_batch.fs <- extr_batch.f[stems]
extr_date.fs <- extr_date.f[stems]
rmsd.fs <- rmsd.f[stems]
covars.s<-list(array_batch.fs,indiv.fs,gender.fs,extr_date.fs,extr_batch.fs) #leave out repr batch for lcls, will bug later
covars_names.s <- factor(c("array_batch.fs","indiv.fs","gender.fs","extr_date.fs","extr_batch.fs"))

#Subset LCLs and put your covariates in a list
name.fl =name.f[lcls]
indiv.fl = indiv.f[lcls]
type.fl <- type.f[lcls]
ID.fl <- ID.f[lcls]
repr_batch.fl <- repr_batch.f[lcls]
array_batch.fl <- array_batch.f[lcls]
gender.fl <- gender.f[lcls]
extr_batch.fl <- extr_batch.f[lcls]
extr_date.fl <- extr_date.f[lcls]
rmsd.fl <- rmsd.f[lcls]
covars.l<-list(array_batch.fl,indiv.fl,gender.fl,extr_date.fl,extr_batch.fl) #leave out repr batch for lcls, will bug later
covars_names.l <- factor(c("array_batch.fl","indiv.fl","gender.fl","extr_date.fl","extr_batch.fl"))

#Take only the probes that have a detection p-value<.05 in at least one individual
detect_quant.all= rowSums(data.lumi@assayData$detection<0.05) #47,311
detect.ind.all <- which(detect_quant.all > 1) #31,945
data.lumi <- data.lumi[detect.ind.all,]

###Find the column that is lumi_ID in feature data usually column 1
head(data.lumi@featureData[[5]]) ## Should read: [1] "ILMN_1762337" "ILMN_2055271" "ILMN_1736007" "ILMN_2383229" "ILMN_1806310" "ILMN_1779670"....

###Subset expression by probes with no CEU Hapmap SNPs.
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
probes = as.character(goodprobes$probeID) ## Convert from factor to character
data.lumi.clean = data.lumi[data.lumi@featureData[[5]] %in% probes, ] #featureData[[5]] should be probe names... Goes from 31,945 to 20,224


### NORMALIZATION ###
# Convert raw Illumina probe intensities to expression values (Both iPSCs and LCLs, no hearts)
# Corrects background, log2 stabiliizes variance, and quantile normalize
data.norm.all <- lumiExpresso(data.lumi.clean, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
expr_quant.all <- data.norm.all@assayData$exprs
dim(expr_quant.all)

###Convert expr_quant rownames to these lumi IDs
rownames(expr_quant.all)=data.norm.all@featureData[[5]]
#Label your columns by sample name
samplenames = read.delim('YGilad-ST sample names switched 8-28.txt', header=TRUE)
colnames(expr_quant.all) = samplenames[(both),1]

#Subtract out mean. Make sure this is what you want to do.
#expr_stem <- rbind(expr_quant.all, type.fb)
expr_stem <- expr_quant.all
dim(expr_stem)
expr_stem_s <- expr_stem[,c(type.fb=="i")]
dim(expr_stem_s)
expr_stem_l <- expr_stem[,c(type.fb=="L")]
dim(expr_stem_l)

####################################################################################################################
## This next section subsets your data to include only 1 probe per gene.The 3' most probe. 
#expr_stem <- expr_stem[-nrow(expr_stem),]
expr_genes <- expr_stem

# Convert probe IDs to gene names using Darren's file (has HGNC & ensembl)
gene_names=c()
for(i in 1:dim(expr_stem)[1]){ 
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_stem)[i],8])) #creates a list of gene names the same length as expression data
}
rownames(expr_genes)=gene_names
Unique_genes = unique(rownames(expr_genes))
length(Unique_genes) #this will be the length of your final data set. One probe per gene.

## This loop will give the most 3' value for multiple probes within the same gene. 
expr_gene = matrix(NA, ncol=length(both), nrow=length(Unique_genes))
i=0
d=0
e=0
f=0
for(gene in unique(gene_names)){
  i = i+1
  
  currRows = which(gene_names == gene)
  if(length(currRows)>1){
    if(goodprobes[currRows[1],6]=="+"){
      keepRow = currRows[which.max(goodprobes[currRows,2])]
      #print(keepRow)
      d=d+1
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
      #print(keepRow)
      e=e+1
    }
  }
  else{
    keepRow=currRows[1]
    f=f+1
  }
  expr_gene[i,] = expr_stem[keepRow,]
  
} 
d #1952 + strand multiple probes
e #1823 - strand multiple probes
f #11531= single probe

dim(expr_gene) #From 20224 probes to 15306 genes
rownames(expr_gene) = Unique_genes
colnames(expr_gene) = colnames(expr_stem)


#Regress out arrray
abatch_all <- ComBat(expr_gene, array_batch.fb, mod=NULL)
cor.abatch.int.g <- cor(abatch_all,method="pearson", use="complete.obs")
heatmap.2(cor.abatch.int.g, Rowv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))),
          Colv=as.dendrogram(hclust(as.dist(1-cor.abatch.int.g))), margins=c(5,9),key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")


#Making dendrograms
cor <- cor(abatch_all, method="pearson")
dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)
hc <- as.dendrogram(hc)
color <- indiv.fb
plot(hc, nodePar = list(col =c("red", "green")))
plot(hc, main = "Cluster Dendrogram: Both Cell types", cex.main = 1.5 , nodePar = list(pch = c(1,NA), cex = 0.8, col=indiv.fb), col = "#487AA1", col.main = "#45ADA8", col.lab = "lightcoral", col.axis = "#F38630", lwd = 3, lty = 1, sub = "", hang = -1, axes = FALSE)
axis(side = 2, at = seq(0.0, 1.4, .2), col = "#F38630", labels = FALSE, lwd = 2)
# add text in margin
mtext(seq(0, 1.4, .2), side = 2, at = seq(0, 1.4, .2), line = 1, col = "#A38630", las = 2)

# vector of colors labelColors = c('red', 'blue', 'darkgreen', 'darkgrey',
# 'purple')
labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951", "red", "bllue")
# cut dendrogram in 4 clusters
clusMember = cutree(hc, 6)
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
# using dendrapply
clusDendro = dendrapply(hc, colLab)
# make plot
plot(clusDendro, main = "Cool Dendrogram", type = "triangle")


#Now prep for separate analysis in iPSCs and LCLs by splitting the dataset.
abatch_bind <- rbind(abatch_all, type.fb)
dim(abatch_bind)
abatch_stem <- abatch_bind[,c(type.fb=="i")]
abatch_stem <- abatch_stem[-nrow(abatch_stem),]
dim(abatch_stem)
abatch_lcl<- abatch_bind[,c(type.fb=="L")]
abatch_lcl <- abatch_lcl[-nrow(abatch_lcl),]
dim(abatch_lcl)


 ##############################################################################################################
# Stems Analysis

#Heatmap
cor.stem <- cor(abatch_stem,method="pearson", use="complete.obs")
heatmap.2(cor.stem, Rowv=as.dendrogram(hclust(as.dist(1-cor.stem))),
          Colv=as.dendrogram(hclust(as.dist(1-cor.stem))), margins=c(5,9),key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")

#Dendrogram

cor <- cor(abatch_stem, method="pearson")
dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)
#hc <- as.dendrogram(hc)
#hccol <- rep("red", times=17)
#indiv_vec <- as.vector(indiv.fs)
hc_col <- c("red", "black", "forestgreen", "darkorchid2", "navy", "red", "black", "darkorchid2","forestgreen", "darkorange1", "navy", "red", "black","darkorchid2", "forestgreen", "darkorange1", "navy")
#hc_col <- c("red", "black")
hc_col
hc_names <- as.vector(name.fs)
plotColoredClusters(hd=hc, labs=hc_names, cols=hc_col, lwd=2, font.lab=2, font.axis=2, font.sub=2, cex=1.2, sub = "", xlab="") #lty=3)
hc_order <- c(1,3,2,4,6,5,8,7,9,10,11,12,13,14,15,16,17)
hc_order <- c(1:17)
plot(reorder(hc, hc_order))
plot(hc)
color <- indiv.fs
plot(hc, main = "Cluster Dendrogram: Stems", cex.main = 1.5 , col = "#487AA1", col.main = "#45ADA8", col.lab = "lightcoral", col.axis = "#F38630", lwd = 3, lty = 1, sub = "", hang = -1) #, axes = FALSE)
axis(side = 2, at = seq(0.0, 1.4, .2), col = "#F38630", labels = FALSE, lwd = 2)
# add text in margin
mtext(seq(0, 1.4, .2), side = 2, at = seq(0, 1.4, .2), line = 1, col = "#A38630", las = 2)

#Relationship between PCs and covariates for regressed data
npcs = 4
sum.PC <- prcomp(na.omit(abatch_stem), scale=TRUE)
results<-c()
for (f in covars.s) {
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
rownames(resultsM_gene_corrected) = covars_names.s
colnames(resultsM_gene_corrected) = c("PC1 p value","PC1 R2","PC2 p value","PC2 R2","PC3 p value","PC3 R2","PC4 p value","PC4 R2")
resultsM_gene_corrected
PC_table_stem <- resultsM_gene_corrected


#PCA plots for regressed data
sum.PC <- prcomp(na.omit(abatch_stem), scale=TRUE)
sumsum <- summary(sum.PC)
op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
title.PC = "PCA of Gene Expression in Stemss"

#prints out plots in c(#rows, #columns)
color = indiv.fs
par(mfrow = c(1,1),oma=c(0,0,2,0)) 
plot(c(1:ncol(abatch_stem)),sum.PC$rotation[,1],cex=1.5,col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:ncol(abatch_stem)),sum.PC$rotation[,1], indiv.fs, cex = 0.55, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC$rotation[,1], sum.PC$rotation[,i],labels=indiv.fs, cex = 0.8, pos=3)   
}  

#Within individual pearson correlation coefficient
cor_2s <- c(cor[1,6],cor[1,12],cor[6,12])
cor_5s <- c(cor[2,7],cor[2,13],cor[7,13])
cor_6s <- c(cor[4,8],cor[8,14],cor[4,14]) ## CHANGE TO STEM
cor_9s <- c(cor[3,9],cor[3,15],cor[9,15]) ## CHANGE TO STEM
cor_10s <- c(cor[10,16])
cor_14s <- c(cor[5,11],cor[5,17],cor[11,17])

cor_within_s <- qpcR:::cbind.na(cor_2s,cor_5s,cor_6s,cor_9s,cor_10s,cor_14s)

cor_wmeans_s <- apply(cor_within_s, 2, mean)

# across individual correlation coefficients
cor_plus <- cbind(cor, indiv.fs)
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
PC12 <- cbind(sum.PC$rotation[,1], sum.PC$rotation[,2])

# Within Individual
dist12_ind <- c()
ind1 <- PC12[c(1,6,12),]
ind2 <- PC12[c(2,7,13),]
ind3 <- PC12[c(4,8,14),]
ind4 <- PC12[c(3,9,15),]
ind5 <- PC12[c(10,16),]
ind6 <- PC12[c(5,11,17),]

md1 <- mean(dist(ind1))
md2 <- mean(dist(ind2))
md3 <- mean(dist(ind3))
md4 <- mean(dist(ind4))
md5 <- mean(dist(ind5))
md6 <- mean(dist(ind6))
mdS <- c(md1,md2,md3,md4,md5,md6)
mdS_mean <- mean(mdS)

# Across Individual

dist_PC12_S <- as.matrix(dist(PC12))
dist_PC12_S

dist_plus <- cbind(dist_PC12_S, indiv.fs)

dist_plus <- rbind(dist_plus, indiv.fs)
dist_plus

dist_stem <- dist_plus
dist_all_stem <- c()
for (i in 1:(nrow(dist_stem)-1)) {
  for (j in 1:(ncol(dist_stem)-1)) {
    #print(i)
    if(dist_stem[i,18]!= dist_stem[18,j]) {
      dist_all_stem <- c(dist_all_stem, dist_stem[i,j])
    }
    else {print(i)}
  }
}
dist_all_stem
mean_dist_all_stem <- mean(dist_all_stem)


##LCLs


cor.lcl <- cor(abatch_lcl,method="pearson", use="complete.obs")
heatmap.2(cor.lcl, Rowv=as.dendrogram(hclust(as.dist(1-cor.lcl))),
          Colv=as.dendrogram(hclust(as.dist(1-cor.lcl))), margins=c(5,9),key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")

cor <- cor(abatch_lcl, method="pearson")
dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)

cor <- cor(abatch_lcl, method="pearson")
dis  <- 1-cor
distance <- as.dist(dis)
hc <- hclust(distance)
hc <- as.dendrogram(hc)
#hc <- as.dendrogram(hc)
#hccol <- rep("red", times=17)
#indiv_vec <- as.vector(indiv.fs)
hc_col <- c("red", "black","darkorchid2", "forestgreen",  "navy", "red", "black", "darkorchid2","forestgreen", "darkorange1", "navy", "red", "black","darkorchid2", "forestgreen", "darkorange1", "navy")
#hc_col <- c("red", "black")
hc_col
hc_names <- as.vector(name.fl)
plotColoredClusters(hc, labs=hc_names, cols=hc_col, lwd=2, font.lab=2, font.axis=2, font.sub=2, cex=1.2, sub = "", xlab="") #lty=3)
 par(xaxs="r",yaxs="r")
#Relationship between PCs and covariates for regressed data
hc_order <- c(1,3,2,4,6,5,8,7,9,10,11,12,13,14,15,16,17)
hc_order <- c(1:17)
hc_order <- c(1:7)
plot(reorder(as.dendrogram(hc), hc_order))
reorder.dendrogram(hc, hc_order)
rhc<- reorder(hc,hc_order)
plot(hc)
plot(rhc)

npcs = 4
sum.PC <- prcomp(na.omit(abatch_lcl), scale=TRUE)
results<-c()
for (f in covars.s) {
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
rownames(resultsM_gene_corrected) = covars_names.l
colnames(resultsM_gene_corrected) = c("PC1 p value","PC1 R2","PC2 p value","PC2 R2","PC3 p value","PC3 R2","PC4 p value","PC4 R2")
resultsM_gene_corrected
PC_table_lcl <- resultsM_gene_corrected


#PCA plots for regressed data
sum.PC <- prcomp(na.omit(abatch_lcl), scale=TRUE)
sumsum <- summary(sum.PC)
op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
title.PC = "PCA of Gene Expression in lcls"

#prints out plots in c(#rows, #columns)
color = indiv.fl
par(mfrow = c(1,1),oma=c(0,0,2,0)) 
plot(c(1:ncol(abatch_stem)),sum.PC$rotation[,1],cex=1.5,col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:ncol(abatch_stem)),sum.PC$rotation[,1], indiv.fl, cex = 0.55, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC$rotation[,1], sum.PC$rotation[,i],labels=indiv.fl, cex = 0.8, pos=3)   
}  

# plot(sum.PC$rotation[,1], sum.PC$rotation[,2],cex=1.5, col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,2]*100),"% of variance", sep=" "))
# text(sum.PC$rotation[,1], sum.PC$rotation[,2],labels=indiv.fl, cex = 0.8, pos=3) 
# text(.245, 0.4, labels=c("Euclidean Distance of Clusters:\n p=\n",ED_t_lclx$p.value))

#Within individual pearson correlation coefficient
cor_2l <- c(cor[1,6],cor[1,12],cor[6,12])
cor_5l <- c(cor[2,7],cor[2,13],cor[7,13])
cor_6l <- c(cor[3,8],cor[8,14],cor[3,14])
cor_9l <- c(cor[4,9],cor[4,15],cor[9,15])
cor_10l <- c(cor[10,16])
cor_14l <- c(cor[5,11],cor[5,17],cor[11,17])

cor_within_l <- qpcR:::cbind.na(cor_2l,cor_5l,cor_6l,cor_9l,cor_10l,cor_14l)

cor_wmeans_l <- apply(cor_within_l, 2, mean)

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
ind6 <- PC12_L[c(5,11,17),]

md1 <- mean(dist(ind1))
md2 <- mean(dist(ind2))
md3 <- mean(dist(ind3))
md4 <- mean(dist(ind4))
md5 <- mean(dist(ind5))
md6 <- mean(dist(ind6))
mdL <- c(md1,md2,md3,md4,md5,md6)
mdL_mean <- mean(mdL)

#Dist across individuals

dist_PC12_L <- as.matrix(dist(PC12_L))
dist_PC12_L

dist_plus_L <- cbind(dist_PC12_L, indiv.fl)

dist_plus_L <- rbind(dist_plus_L, indiv.fl)
dist_plus_L

dist_lcl <- dist_plus_L
dist_all_lcl <- c()
for (i in 1:(nrow(dist_lcl)-1)) {
  for (j in 1:(ncol(dist_lcl)-1)) {
    #print(i)
    if(dist_lcl[i,18]!= dist_lcl[18,j]) {
      dist_all_lcl <- c(dist_all_lcl, dist_lcl[i,j])
    }
    else {print(i)}
  }
}
dist_all_lcl
mean_dist_all_lcl <- mean(dist_all_lcl)




######################################################################################
#Differential expression with limma

design <- model.matrix(~0+type.fb)

colnames(design) <- c("heart", "iPSC", "LCL")

design

fit <- lmFit(abatch_all, design)

contrast.matrix<-makeContrasts(iPSC-LCL,levels=design)

fit2<-contrasts.fit(fit,contrast.matrix)

fit2<-eBayes(fit2)

output <- topTable(fit2, number=20000, p.value=.01)

output
head(output)

adjust <- output[,5]
head(adjust)
length(adjust)
adjust

result <- decideTests(fit2)

result

######################################################################################
#Distance between all
ED_t_stem <- t.test(mdS, dist(PC12))
Edist_stem <- c(mdS_mean, mean(dist(PC12)), ED_t_stem$p.value)
names(Edist_stem)<- c("ED12 within", "ED12 between all", "p-value")
Edist_stem

ED_t_lcl <- t.test(mdL, dist(PC12_L))
Edist_lcl <- c(mdL_mean, mean(dist(PC12_L)), ED_t_lcl$p.value)
names(Edist_lcl)<- c("ED12 within", "ED12 between", "p-value")
Edist_lcl

#Distance between different clsuters
ED_t_stemx <- t.test(mdS, dist_all_stem)
Edist_stemx <- c(mdS_mean, mean_dist_all_stem, ED_t_stemx$p.value)
names(Edist_stemx)<- c("ED12 within", "ED12 between all", "p-value")
Edist_stemx

ED_t_lclx <- t.test(mdL, dist_all_lcl)
Edist_lclx <- c(mdL_mean, mean_dist_all_lcl, ED_t_lclx$p.value)
names(Edist_lclx)<- c("ED12 within", "ED12 between", "p-value")
Edist_lclx

#combine across
Edist <- cbind(Edist_stemx, Edist_lclx)
colnames(Edist) <- c("iPSCs", "LCLs")
Edist_t <- t(Edist)
Edist_t

#Correlations: mean for each individual and overall
cor_wmeanl <- apply(cor_within_l,2, mean)
cor_wmeans <- apply(cor_within_s,2, mean)
cor_both <- cbind(cor_within_s, cor_within_l)
cor_bmeans <- cbind(cor_wmeans, cor_wmeanl)
cor_all_both <- cbind(cor_all_lcl, cor_all_stem)
boxplot(cor_all_both)
boxplot(cor_both)
boxplot(cor_bmeans, main="Within Individual Pearson correlation coefficients for Stem cells and LCLs")

cor_total <- cbind(cor_wmeans, cor_wmeanl, cor_all_stem, cor_all_lcl)
cor_total_vec <- c(cor_wmeans, cor_wmeanl, cor_all_stem, cor_all_lcl)
cor_within_l.long <- c(as.vector(cor_within_l), rep(NA, length(cor_all_lcl)-length(cor_within_l)))
cor_within_s.long <- c(as.vector(cor_within_s), rep(NA, length(cor_all_stem)-length(cor_within_s)))
cor_tots <- cbind(cor_within_l.long, cor_within_s.long, cor_all_lcl, cor_all_stem)

boxplot(cor_tots) #between all samples
boxplot(cor_total) #mean for each group

mean(cor_wmeans)
mean(cor_wmeanl)
median(cor_wmeans)
median(cor_wmeanl)

# Get P values
t.test(cor_within_l, cor_all_lcl)
t.test(cor_within_s, cor_all_stem)
t.test(cor_wmeanl, cor_all_lcl)
t.test(cor_wmeans, cor_all_stem)

t_cor_w <- t.test(as.vector(cor_within_s), as.vector(cor_within_l))
t_cor_w
pv_w <- signif(t_cor_w$p.value,3)
pv_w
t_cor_wmean <- t.test(cor_wmeanl, cor_wmeans)
pv_wm <- signif(t_cor_wmean$p.value, 3)
pv_wm
t_cor_b <- t.test(cor_all_lcl, cor_all_stem)
pv_b <- signif(t_cor_b$p.value, 3)
pv_b

#Make nice boxplot
type <- c(rep("stem", times=length(cor_wmeans)), rep("lcl", times=length(cor_wmeanl)), rep("stem", times=length(cor_all_stem)), rep("lcl", times=length(cor_all_lcl)))
group <- c(rep(as.factor("Within Individuals"), times=(length(cor_wmeans)+length(cor_wmeanl))), rep(as.factor("Between Individuals"), times=(length(cor_all_stem)+length(cor_all_lcl))))
groupn <- c(rep("1", times=(length(cor_wmeans)+length(cor_wmeanl))), rep("2", times=length(cor_all_stem)+length(cor_all_lcl)))

cor_df <- data.frame(cor_total_vec, type, group, groupn)

ggplot(cor_df, aes(x=groupn, y=cor_total_vec)) + geom_boxplot(aes(fill=type), xlab=FALSE) + theme(text = element_text(size=18)) + annotate(geom = "text", label=paste("**p =", pv_w), x=1, y=0.98) + annotate(geom = "text", label=paste("**p =", pv_b), x=2, y=.97) + scale_x_discrete(labels=c("Within Individuals", "Between Individuals")) + theme(axis.title.x=element_blank(), panel.background=element_rect(fill='white')) + ylab("Pearson Correlation") + theme(axis.title.y = element_text(vjust=1.0)) #stat_boxplot(geom ='errorbar', aes(x=group))  

## Coefficient of variatatiionn
#Between samples coefficient of variation
cv <- function(x) (sd(x)/mean(x))
CV_S <- apply(abatch_stem, 1, cv)
CV_L <- apply(abatch_lcl, 1, cv)

mean(CV_S)
mean(CV_L)
t.test(CV_S, CV_L)

## Calculate CV between individuals ##
# Pick a random line from each individual
df_astem <- data.frame((abatch_stem))
df_s1 <- df_astem[,(indiv.fs==2)]
head(df_s1)
df_s2 <- df_astem[,(indiv.fs==5)]
df_s3 <- df_astem[,(indiv.fs==6)]
df_s4 <- df_astem[,(indiv.fs==9)]
df_s5 <- df_astem[,(indiv.fs==10)]
df_s6 <- df_astem[,(indiv.fs==14)]


df_alcl <- data.frame((abatch_lcl))
df_l1 <- df_alcl[,(indiv.fl==2)]
df_l2 <- df_alcl[,(indiv.fl==5)]
df_l3 <- df_alcl[,(indiv.fl==6)]
df_l4 <- df_alcl[,(indiv.fl==9)]
df_l5 <- df_alcl[,(indiv.fl==10)]
df_l6 <- df_alcl[,(indiv.fl==14)]

head(df_l5)

random_stem <- cbind(sample(df_s1,1),sample(df_s2,1), sample(df_s3,1), sample(df_s4,1),sample(df_s5,1),sample(df_s6,1))
length(random_stem)
names(random_stem)
dim(random_stem)

#random_stem_m <- as.matrix(unlist(random_stem), nrow=nrow(abatch_stem), byrow=TRUE)
head(random_stem)
dim(random_stem)

random_lcl <- cbind(sample(df_l1,1), sample(df_l2,1), sample(df_l3,1), sample(df_l4,1), sample(df_l5,1), sample(df_l6,1))
length(random_lcl)
names(random_lcl)
head(random_stem)

random_stem_cv <- apply(random_stem, 1, cv)
length(random_stem_cv)
high_CVSR <- random_stem_cv[which(random_stem_cv>.025)]
length(high_CVSR)
head(high_CVSR)

random_lcl_cv <- apply(random_lcl, 1, cv)
length(random_lcl_cv)
high_CVLR <- random_lcl_cv[which(random_lcl_cv>.025)]
length(high_CVLR)
head(high_CVLR)

#Density plots of CVs
cv_all <- data.frame(coefvar=c(CV_L, CV_S), type = rep(c("LCL", "Stem"), times=c(length(CV_L),length(CV_S))))
ggplot(cv_all, aes(x=coefvar, fill=type)) + geom_density(alpha=0.5) +xlim(0,0.05)+xlab("Coefficients of Variation") + ggtitle("Gene Expression: Coefficients of Variation")# + theme(legend.position=c(.75,.75)) +annotate(geom = "text", label="**p-value = 2.2e-16", x=1.9, y=.9)+theme(text = element_text(size=23)) 

high_CVS <- CV_S[which(CV_S>.025)]
head(high_CVS)
length(high_CVS)

high_CVL <- CV_L[which(CV_L>.025)]
head(high_CVL)
length(high_CVL)

high_CVLS <- intersect(names(high_CVS), names(high_CVL))
head(high_CVLS)
length(high_CVLS)

#Change CV to Random CV
CV_L <- random_lcl_cv
CV_S <- random_stem_cv

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
mean_lcl <- apply(expr_LCL, 1, mean)
 plot(mean_lcl, var_lcl)

 #Calculate ratio of variance btw and within individuals for every gene
var_ratio_lcl <- vector()
for(i in 1:length(var_lcl)) {
  var_ratio_lcl <- c(var_ratio_lcl, (var_lcl[i]/(mean(var_win_lcl[i,]))))  
}
length(var_ratio_lcl)

var_l <- cbind(var_win_lcl, var_lcl)
boxplot(var_l, ylim=c(0,.1))
t.test(var_win_lcl, var_lcl)
mean(var_win_lcl)
mean(var_lcl)

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
mean_stem <- apply(expr_stem, 1, mean)
 var_stem <- log(var_stem)
plot(mean_stem, cv_stem)
 
 cv_stem <- var_stem/mean_stem
 
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

#Order dataframe
ordereds<- stems_list[order(stems_list[,1]),]
orderedl<- lcl_list[order(lcl_list[,1]),]  

#Reverse order
rordereds<- stems_list[order(stems_list[,1], decreasing=TRUE),]
rorderedl<- lcl_list[order(lcl_list[,1], decreasing=TRUE),]  

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
  s <- summary(lm(abatch_lcl[i,]~indiv.fl));
  pvals_lcls <- c(pvals_lcls,pf(s$fstatistic[[1]], s$fstatistic[[2]], s$fstatistic[[3]], lower.tail = FALSE))
}

hist(pvals_stems)
hist(pvals_lcls)

adjust_stems <- p.adjust(pvals_stems, method="BH")
length(which(adjust_stems<.05))

adjust_lcls <- p.adjust(pvals_lcls, method="BH")
length(which(adjust_lcls<.05))

hist(adjust_stems)
hist(adjust_lcls)

summary(pvals_lcls)

#plot(stem_var_rat, pvals_stems, xlim=c(0,5))
#plot(lcl_var_rat, pvals_lcls)

mean(pvals_stems)
mean(pvals_lcls)

#count significant genes
sig_stems <- length(adjust_stems[adjust_stems<.05])
sig_stems
sig_lcls <- length(adjust_lcls[adjust_lcls<.05])
sig_lcls

cutoff_stems <- rordereds[sig_stems,]
cutoff_lcls <- rorderedl[sig_lcls,]
cutoff_stems
cutoff_lcls
cutoff_avg <- mean(c(cutoff_stems[,1], cutoff_lcls[,1]))
cutoff_avg


#Subset top and bottom
tops <- ordereds[c(1:sig_stems),]
topl <- orderedl[c(1:sig_lcls),]
bottom_stem <- rordereds[c(1:sig_stems),]
bottom_lcl <- rorderedl[c(1:sig_lcls),]
bottom_stemnames <- row.names(bottom_stem)
bottom_lclnames <- row.names(bottom_lcl)

dim(tops)
dim(topl)

#Overlap of top and bottom between cell types
overlap <- intersect(tops[,2], topl[,2])
length(overlap)
roverlap <- intersect(bottom_stem[,2], bottom_lcl[,2])
length(roverlap)

write(bottom_lclnames, "C:/Users/a a/Documents/Lab/Variation Recovery/high var lcl.txt", sep="\t")
write(bottom_stemnames, "C:/Users/a a/Documents/Lab/Variation Recovery/high var stem.txt", sep="\t")
write(roverlap, "C:/Users/a a/Documents/Lab/Variation Recovery/roverlap.txt", sep="\t")

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
ggplot(var_all, aes(x=var, fill=Cell_type)) + geom_density(alpha=0.5) +xlim(-.25,8.0)+xlab("Variance Between/Variance Within") + geom_vline(xintercept=cutoff_avg, linetype="dotted") + theme(legend.position=c(.75,.75), panel.background=element_rect(fill='white')) + theme(text = element_text(size=18)) +annotate(geom = "text", label=paste("p-value = ", pv_var_rat), x=6, y=.9) + scale_fill_manual(values=cols, labels=c("LCLs", "iPSCs"))

#Absolute

#pretty <- rev(pretty)
# var_all <- data.frame(var=c(var_lcl, var_stem), type = rep(c("LCL", "Stem"), times=c(length(lcl_var_rat),length(stem_var_rat))))
var_all <- data.frame(var=c(var_stem, var_lcl), type = rep(c("iPSC", "LCL"), times=c(length(stem_var_rat),length(lcl_var_rat))), Cell_type = rep(c("1", "2"), times=c(length(stem_var_rat),length(lcl_var_rat))))
ggplot(var_all, aes(x=var, fill=Cell_type)) + geom_density(alpha=0.5) + annotate(geom = "text", label=paste("p-value = ", pv_var), x=.075, y=50) +geom_density(alpha=0.5) +xlim(-.01,.1)+xlab("Variance")  + theme(legend.position=c(.75,.75), panel.background=element_rect(fill='white')) +theme(text = element_text(size=18)) + scale_fill_manual(values=rev(cols), labels=c("iPSCs", "LCLs"))



#eQTL enrichment
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

lcls_lcls <- orderedl[rownames(orderedl) %in% eQTL_lcls,]
dim(lcls_lcls)

elcl_stemvar <- var_stem[names(var_stem) %in% eQTL_lcls]
cv_elcl_stem <- CV_S[names(CV_S) %in% eQTL_lcls]
length(elcl_stemvar)
length(cv_elcl_stem)
mean(cv_elcl_stem)

elcl_lclvar <- var_lcl[names(var_lcl) %in% eQTL_lcls]
cv_elcl_lcl <- CV_L[names(CV_L) %in% eQTL_lcls]
length(elcl_lclvar)
length(cv_elcl_lcl)
mean(cv_elcl_lcl)
mean(cv_elcl_stem)
mean(CV_L)
mean(CV_S)
mean(elcl_lclvar)
mean(elcl_stemvar)
mean(var_stem)
mean(var_lcl)
cv_nelcl_stem <- CV_S[!(names(CV_S) %in% eQTL_lcls)]
cv_nelcl_lcl <- CV_L[!(names(CV_L) %in% eQTL_lcls)]

#bothvar_lcl_eQTLs <- c(mean_elcl_stemsrat, mean_elcl_lclrat, mean_elcl_stemvar, mean_elcl_lclvar)
#bothvar_lcl_eQTLs


#Plot Variance
dodge <- position_dodge(width=0.9)
var_lcl_eQTLs <- c(elcl_stemvar, elcl_lclvar, var_stem, var_lcl)
length(var_lcl_eQTLs)
var_lcl_eQTLs <- sapply(var_lcl_eQTLs, log)
names_var <- c(rep("iPSC", times=length(elcl_stemvar)), rep("LCL", times=length(elcl_lclvar)), rep("iPSC", times=length(var_stem)), rep("LCL", times=length(var_lcl)))
length(names_var)
type_var <- c((rep("eQTL var",  times=(length(elcl_stemvar)+length(elcl_lclvar)))), (rep("mean var", times=(length(var_stem)+length(var_lcl))))) 
length(type_var)
df_var <- data.frame(var_lcl_eQTLs, names_var, type_var, fill=TRUE)
df_var_se <- summarySE(df_var, measurevar="var_lcl_eQTLs", groupvars=c("type_var","names_var"))
dodge <- position_dodge(width=0.9)
ggplot(df_var_se, aes(type_var,var_lcl_eQTLs, fill=names_var)) +geom_bar(stat="identity", position=dodge) + geom_errorbar(aes(ymin=var_lcl_eQTLs-se, ymax=var_lcl_eQTLs+se), width=0.2, position=dodge)
ggplot(df_var, aes(type_var, var_lcl_eQTLs, fill=names_var)) + geom_boxplot(aes(fill=names_var), position=dodge)

#Plot coefficient of variation
cv_lcleQTLs <- c(cv_elcl_stem, cv_elcl_lcl, CV_S, CV_L)
cv_log2 <- log2(cv_lcleQTLs)
df_cv <- data.frame(cv_log2, names_var, type_var)
df_cv_se <- summarySE(df_cv, measurevar="cv_lcleQTLs", groupvars=c("type_var", "names_var"))
ggplot(df_cv, aes(type_var,cv_log2, fill=names_var)) +geom_boxplot(aes(fill=names_var), position=dodge) #+ geom_errorbar(aes(ymin=cv_lcleQTLs-se, ymax=cv_lcleQTLs+se), width=0.2, position=dodge)

t.test(elcl_stemvar, elcl_lclvar)
t.test(elcl_stemvar, var_stem)
t.test(CV_S, CV_L)
t.test(cv_elcl_stem, CV_S)
t.test(cv_elcl_lcl, CV_L)
t.test(cv_elcl_stem, cv_elcl_lcl)
t.test(cv_elcl_stem, cv_nelcl_stem)
t.test(cv_elcl_lcl, cv_nelcl_lcl)
t.test(cv_nelcl_lcl, cv_nelcl_stem)

#Plot coefficient of variation with vs without eqtls
cv_nlcleQTLs <- c(cv_elcl_stem, cv_elcl_lcl, cv_nelcl_stem, cv_nelcl_lcl)
cv_nlog2 <- log2(cv_nlcleQTLs)
names_var <- c(rep("iPSC", times=length(cv_elcl_stem)),  rep("LCL", times=length(cv_elcl_lcl)), rep("iPSC", times=length(cv_nelcl_stem)), rep("LCL", times=length(cv_nelcl_lcl)))
type_var <- c((rep("eQTL var",  times=(length(cv_elcl_stem)+length(cv_elcl_lcl)))), (rep("no eQTL var", times=(length(cv_nelcl_stem)+length(cv_nelcl_lcl))))) 
df_ncv <- data.frame(cv_nlog2, names_var, type_var)
df_ncv_se <- summarySE(df_ncv, measurevar="cv_lcleQTLs", groupvars=c("type_var", "names_var"))
ggplot(df_ncv, aes(type_var,cv_nlog2, fill=names_var)) +geom_boxplot(aes(fill=names_var), position=dodge) #+ geom_errorbar(aes(ymin=cv_lcleQTLs-se, ymax=cv_lcleQTLs+se), width=0.2, position=dodge)


write.table(bottom_lclnames, "C:/Users/a a/Documents/Lab/Variation Recovery/ensembl_lcl_names.txt", sep="\t")
write.table(bottom_stemnames, "C:/Users/a a/Documents/Lab/Variation Recovery/ensembl_top_1500_lcl.txt", sep="\t")


#Now for all
liver_eQTLs <- read.table(c("liver eQTLs.tab"), fill=TRUE)
pons_eQTLs <- read.table(c("Brain Pons eQTLs.tab"), fill=TRUE)
fcortex_eQTLs <- read.table(c("Brain Frontal Cortex eQTLs.tab"), fill=TRUE)
cerebellum_eQTLs <- read.table(c("Brain Cerebellum eQTLs.tab"), fill=TRUE)
tcortex_eQTLs <- read.table(c("Brain Temporal Cortex eQTLs.tab"), fill=TRUE)
all_eQTLs <- read.table(c("non-lcl tissues eQTLs.tab"), fill=TRUE)

#dim(liver_eQTLs) <- read.table(c("Liver eQTLs.tab"), fill=TRUE)
all_eQTLs <- all_eQTLs
all_eQTLs <- unique(all_eQTLs[,9])
length(all_eQTLs)

eall_stemvar <- var_stem[names(var_stem) %in% all_eQTLs]
cv_eall_stem <- CV_S[names(CV_S) %in% all_eQTLs]
cv_neall_stem <- CV_S[!(names(CV_S) %in% all_eQTLs)]
length(eall_stemvar)
length(cv_eall_stem)
mean(cv_eall_stem)

eall_lclvar <- var_lcl[names(var_lcl) %in% all_eQTLs]
cv_eall_lcl <- CV_L[names(CV_L) %in% all_eQTLs]
cv_neall_lcl <- CV_L[!(names(CV_L) %in% all_eQTLs)]
length(eall_lclvar)
length(cv_eall_lcl)
mean(cv_eall_lcl)
mean(cv_eall_stem)
mean(CV_L)
mean(CV_S)

var_alleQTLs <- c(eall_stemvar, eall_lclvar, var_stem, var_lcl)

cv_alleQTLs <- c(cv_eall_stem, cv_eall_lcl, CV_S, CV_L)
cv_alleQTLs <- log2(cv_alleQTLs)
names_var_all <- c(rep("iPSC", times=length(cv_eall_stem)),  rep("LCL", times=length(cv_eall_lcl)), rep("iPSC", times=length(CV_S)), rep("LCL", times=length(CV_L)))
type_var_all <- c((rep("eQTL CV",  times=(length(cv_eall_stem)+length(cv_eall_lcl)))), (rep("mean CV", times=(length(CV_S)+length(CV_L))))) 

dodge <- position_dodge(width=0.9)
df_cv_all <- data.frame(cv_alleQTLs, names_var_all, type_var_all)
df_cv_all_se <- summarySE(df_cv_all, measurevar="cv_alleQTLs", groupvars=c("type_var_all", "names_var_all"))
ggplot(df_cv_all, aes(type_var_all, cv_alleQTLs, fill=names_var_all)) +geom_boxplot(aes(fill=names_var_all), position=dodge) #+ geom_errorbar(aes(ymin=cv_lcleQTLs-se, ymax=cv_lcleQTLs+se), width=0.2, position=dodge)

var_alleQTLs <- log(var_alleQTLs)
df_var_all <- data.frame(var_alleQTLs, names_var_all, type_var_all)
ggplot(df_var_all, aes(type_var_all, var_alleQTLs, fill=names_var_all)) + geom_boxplot(aes(fill=names_var_all), position=dodge)

t.test(elcl_stemvar, elcl_lclvar)
t.test(elcl_stemvar, var_stem)
t.test(CV_S, CV_L)
t.test(cv_eall_stem, cv_eall_lcl)
t.test(cv_elcl_stem, cv_elcl_lcl)
t.test(cv_eall_stem, CV_S)
t.test(cv_eall_lcl, CV_L)
t.test(cv_eall_lcl, cv_neall_lcl)
t.test(cv_eall_stem, cv_neall_stem)
t.test(cv_elcl_lcl, cv_nelcl_lcl)
t.test(cv_elcl_stem, cv_nelcl_stem)
#all
cv_alleQTLs <- c(cv_eall_stem, cv_elcl_stem, cv_eall_lcl, cv_elcl_lcl, CV_S, CV_L)
log_cv_alleQTLs <- log(cv_alleQTLs) #
#cv_alleQTLs <- c(eall_stemvar, elcl_stemvar, eall_lclvar, elcl_lclvar, var_stem, var_lcl)
names_var_all <- c(rep("iPSC", times=length(cv_eall_stem)+length(cv_elcl_stem)),  rep("LCL", times=length(cv_eall_lcl)+length(cv_elcl_lcl)), rep("iPSC", times=length(CV_S)), rep("LCL", times=length(CV_L)))
type_var_all <- c(rep("all eQTL CV",  times=(length(cv_eall_stem))), rep("LCL eQTL CV", times=length(cv_elcl_stem)), rep("all eQTL CV", times=length(cv_eall_lcl)), rep("LCL eQTL CV", times=length(cv_elcl_lcl)), (rep("mean CV", times=(length(CV_S)+length(CV_L)))))
df_cv_all <- data.frame(log_cv_alleQTLs, names_var_all, type_var_all)
ggplot(df_cv_all, aes(type_var_all, log_cv_alleQTLs, fill=names_var_all)) +geom_boxplot(aes(fill=names_var_all), position=dodge) #+ geom_errorbar(aes(ymin=cv_lcleQTLs-se, ymax=cv_lcleQTLs+se), width=0.2, position=dodge)


#Amount of variation explained by individual
results <- c()
for (i in 1:nrow(abatch_stem)) {
  s = summary(lm(abatch_stem[i,]~indiv.fs));
  results<-c(results,pf(s$fstatistic[[1]],
                        s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
             s$adj.r.squared)
}

resultsM_s <- matrix(ncol=2, data=results, byrow=TRUE)

results <- c()
for (i in 1:nrow(abatch_stem)) {
  s = summary(lm(abatch_stem[i,]~indiv.fs));
  results<-c(results,pf(s$fstatistic[[1]],
                        s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
             s$r.squared)
}
resultsM_s_noadj <- matrix(ncol=2, data=results, byrow=TRUE)
exp_s <- mean(resultsM_s[,2])
exp_s
mean(resultsM_s[,1])
boxplot(resultsM_s[,2])

a <- aov(abatch_stem[1,]~indiv.fs)
a

summary(a)
for (i in 1:nrow(abatch_stem)) {
  s=summary(aov(abatch_stem[i,])~indiv.fs)
  results <- c(s)
}

#LCLs

results <- c()
for (i in 1:nrow(abatch_lcl)) {
  s = summary(lm(abatch_lcl[i,]~indiv.fl));
  results<-c(results,pf(s$fstatistic[[1]],
                        s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
             s$adj.r.squared)
}

resultsM <- matrix(ncol=2, data=results, byrow=TRUE)

results <- c()
for (i in 1:nrow(abatch_lcl)) {
  s = summary(lm(abatch_lcl[i,]~indiv.fl));
  results<-c(results,pf(s$fstatistic[[1]],
                        s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
             s$r.squared)
}

resultsM_noadj <- matrix(ncol=2, data=results, byrow=TRUE)
exp_l <- mean(resultsM[,2])
exp_l
mean(resultsM[,2])
mean(resultsM_noadj[,2])
mean(resultsM_s[,2])
mean(resultsM_s_noadj[,2])

boxplot(resultsM[,1], resultsM_s[,1])

boxplot(resultsM[,2], resultsM_s[,2], resultsM_noadj, resultsM_s_noadj)
barplot(resultsM[,2], resultsM_s[,2])

t.test(resultsM[,2], resultsM_s[,2])
t.test(resultsM[,1], resultsM_s[,1])
t.test(resultsM_noadj, resultsM_s_noadj)


 plot(resultsM[,1], resultsM_s[,1], what="density")
hist(resultsM[,2])

exp_PC1s <- PC_table_stem[2,2]
exp_PC1l <- PC_table_lcl[2,2]
exp_PC2s <- PC_table_stem[2,4]
exp_PC2l <- PC_table_lcl[2,4]
exp_PC3s <- PC_table_stem[2,6]
exp_PC3l <- PC_table_lcl[2,6]
exp_PC4s <- PC_table_stem[2,8]
exp_PC4l <- PC_table_lcl[2,8]

print(c(exp_PC1s, exp_PC1l, exp_PC2s, exp_PC2l))
PC_table_stem
PC_table_lcl

explained_avg <- c(exp_l, exp_s)
df_avg <- data.frame(explained_avg)
explained_PC <- c(exp_PC1s, exp_PC1l, exp_PC2s, exp_PC2l, exp_PC3s, exp_PC3l, exp_PC4s, exp_PC4l)
type_PC <- c("iPSC", "LCL","iPSC", "LCL","iPSC", "LCL","iPSC", "LCL")
num_PC <- c("1", "1", "2", "2", "3", "3", "4", "4")
df_PC <- data.frame(explained_PC, type_PC, num_PC)
ggplot(df_PC, aes(num_PC, explained_PC, fill=type_PC))+geom_bar(stat="identity", position=dodge)
print(explained_PC)
explained <- c(exp_l, exp_s, exp_PC1l, exp_PC1s)
explained_gene <- c(resultsM_s[,2], resultsM[,2])

type <- c(rep("iPSC", times=nrow(resultsM_s)), rep("LCL", times=nrow(resultsM)))
df_exp_gene <- data.frame(explained_gene, type)
dim(df_exp_gene)
colnames(df_exp_gene) <- c("R2", "type")
df_exp_se <- summarySE(df_exp_gene, measurevar="R2", groupvars=c("type"))
df_exp_se
ggplot(df_exp_se, aes(type, R2, fill=type)) + geom_bar( stat="identity",position=dodge)  + geom_errorbar(aes(ymin=R2-se, ymax=R2+se), width=0.2, position=dodge)

df_var_se <- summarySE(df_var, measurevar="var_lcl_eQTLs", groupvars=c("type_var","names_var"))
dodge <- position_dodge(width=0.9)
ggplot(df_var_se, aes(type_var,var_lcl_eQTLs, fill=names_var)) +geom_bar(stat="identity", position=dodge) + geom_errorbar(aes(ymin=var_lcl_eQTLs-se, ymax=var_lcl_eQTLs+se), width=0.2, position=dodge)

#df_var_se <- summarySE(df_var, measurevar="var_lcl_eQTLs", groupvars=c("type_var","names_var"))
#type <- c("iPSC","LCL", "iPSC", "LCL")
type <- c("LCL", "iPSC","LCL", "iPSC")
cat <- c("mean explained", "mean explained","PC1", "PC1")
means_all <- cbind(explained, type)
means_all
df_ma <- data.frame(explained)
df_ma <- cbind (df_ma, type, cat)
row.names(df_ma) <- cat2
df_ma

df_ma_me <- df_ma[c(1,2),]
df_ma_me

ggplot(df_ma_me, fill=type) + geom_bar(aes(row.names(df_ma_me), explained), stat="identity")
ggplot(df_ma, aes(x=cat, y=explained), fill=type) + geom_bar(aes(fill=type), stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(color="black"), axis.title.x=element_blank(), text=element_text(size=18)) + theme(panel.background=element_rect(fill='white'))#+ geom_errorbar(aes(ymin=val-se, ymax=val+se))
ggplot(df_exp_se, aes(y=explained_avg), fill=type) + geom_bar(aes(fill=type), stat="identity", position=position_dodge()) #+ geom_errorbar(aes(ymin=val-se, ymax=val+se))


ggplot(df_ma, aes(x=cat, y=explained), fill=type) + geom_bar(aes(fill=type), stat="identity")

both <- data.frame(var=c(resultsM[,2], resultsM_s[,2]), type = rep(c("LCL", "Stem"), times=c(nrow(resultsM),nrow(resultsM_s))))
ggplot(both, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.25,2.5)+xlab("p-Value") + ggtitle("Gene expression correlation with covariate individual") + theme(legend.position=c(.75,.75)) + theme(text = element_text(size=23)) 


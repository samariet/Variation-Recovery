x <- rnorm(100)
y <- x + 2 + rnorm(100)
plot(x,y)
cor(x,y)
cor(x + 10,y + 10)
plot(x + 10,y + 10)
z <- rpois(100, 10)
plot(x + z, y + z)
cor(x + z, y + z)
range((x + z) - (y + z))
range(x - y)

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
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05) #How many individuals have detection p<.05
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

#relevel(array_batch.fs, "2")
##Regress out array batch and add the intercept back in
abatch.residual.int.g = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_stem))
rownames(abatch.residual.int.g) = rownames(expr_gene)
colnames(abatch.residual.int.g) = colnames(expr_gene)
for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ array_batch.fs, "3")
  abatch.residual.int.g[i,] = resid(model) # + mean(expr_gene[i,])
  #abatch.residual.int.g[i,] = resid(model)
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
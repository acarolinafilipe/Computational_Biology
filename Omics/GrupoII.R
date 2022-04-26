#Grupo II

#a)READ COVERAGE DISTRIBUTION

setwd("/Users/Carolina/Ambiente de Trabalho/BC/Lab 5 data/")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(edgeR)
library(DESeq2)
library(stringr)
library(tibble)
library(ggrepel)
library(GOstats)
library(biomaRt)

#import dataset

TCGA_BRCA_Gene_ReadCounts <- read.table(file = "TCGA_BRCA_Gene_ReadCounts.txt.gz", header=TRUE)

#Read Coverage Depth: soma das colunas, expressao genetica (nr total de reads para cada amostra)

genereads<- as.data.frame(colSums(TCGA_BRCA_Gene_ReadCounts[,-1]))
names(genereads)[names(genereads) == "colSums(TCGA_BRCA_Gene_ReadCounts[, -1])"] <- "Total_counts"

genereads$names <- rownames(genereads)

#Read Coverage Plot
hist(genereads$Total_counts, main= 'Read Coverage Distribution', xlab= 'Read Coverage Depth', ylab='Number of Samples', col="lightblue")

#LIBRARY COMPLEXITY

#first, get the genes and do cumulative sum of them

lib_matrix <- TCGA_BRCA_Gene_ReadCounts
colnames_lib<- colnames(lib_matrix[,-1])

#sorting columns by ascending order
lib_matrix_sorted <- apply(lib_matrix,2,sort,decreasing=T)

#normalize
norm_lib <- function(x) {return ((cumsum(x)/sum(as.numeric(x))))}

lib_matrix_n <- data.frame(apply(lib_matrix_sorted[,-1],2, norm_lib))
#add a column w/ the no. of genes
lib_matrix_n$num_genes <- seq.int(nrow(lib_matrix_n))

#turn all the columns into one
meltedlib <- melt(lib_matrix_n, id.vars = "num_genes")

#plot results
dev.new()
ggplot(meltedlib, aes(x = num_genes, y = value, color = variable)) + 
  geom_line() + 
  xlab(label="Number of Genes") + 
  ylab(label = "Cumulative Proportion of Reads") +
  ggtitle('Library Complexity') +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5))

#Find the sample which  70% of its estimated counts attributed to only 1000 genes
samples_i <- lib_matrix_n[1000,1:878]>0.70
samples_i <- which(samples_i)

samples_data <- data.frame(lib_matrix$Gene,lib_matrix[,samples_i+1])
samples_data <- samples_data[order(-lib_matrix[,samples_i+1]) ,]
samples_data[1:10,]

#b) NORMALIZATION

ReadCounts <- TCGA_BRCA_Gene_ReadCounts
ReadCounts[ReadCounts == 0] <- NA #To deal w/ logarithm, turn zeros into NA

#plot distributions prior to normalization
dev.new()
boxplot(ReadCounts[,2:101], log = "y", main = 'Distribution of Read Counts Before Normalization',
        xlab="Samples", ylab="Estimated Counts")

# Normalize data for each sample, as well a dealing w/bias

#renaming for simplicity and removing rows that are all zero
ReadCounts <- TCGA_BRCA_Gene_ReadCounts[rowSums(TCGA_BRCA_Gene_ReadCounts[, -1]) > 0, ]
rownames(ReadCounts)<- ReadCounts[,1]
ReadCounts<- ReadCounts[,-1]

dge <- DGEList(counts=data.matrix(ReadCounts))
dge <- calcNormFactors(dge)
v <- voom(dge,plot=TRUE)

#plot distributions post-normalization
dev.new()
boxplot(v$E[,1:100], main="Distribution of Read Counts After Normalization",
        xlab="Samples", ylab="Estimated Counts")



#c) PCA

#import annotations:
ClinicalInfo <- read.delim(file = 'TCGA_BRCA_ClinicalAnnotation.txt.gz', header=TRUE)

#Scree Plot
PCA_basic <- prcomp(t(v$E), center = TRUE)
dev.new()
pca.var.per <- PCA_basic$sdev^2 / sum (PCA_basic$sdev^2)
plot(pca.var.per[1:20], main="Scree Plot", xlab="Principal Component", ylab="Variance Explained", type="b" , ylim=c(0,1))

#cumultive sum of PCA
dev.new()
plot(cumsum(pca.var.per[1:800]), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained", main = "Cumulative Principal Components",
     type = "b")

#plot PC1 and PC2
pca.data <- data.frame(Sample=rownames(PCA_basic$x),
                       X=PCA_basic$x[,1],
                       Y=PCA_basic$x[,2])
dev.new()
ggplot(pca.data, aes(x=X, y=Y)) +
  geom_point() +
  xlab(paste("PC1 : ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 : ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle('PCA Analysis')

library(rgl)
pca3d.data <- data.frame(Sample=rownames(PCA_basic$x),
                         X=PCA_basic$x[,1],
                         Y=PCA_basic$x[,2],
                         Z=PCA_basic$x[,3])

dev.new()
plot3d(x=pca3d.data$X, y=pca3d.data$Y, z=pca3d.data$Z, 
       xlab = (paste("PC1 : ", round(pca.var.per[1], digits=4), "%", sep="")),
       ylab=(paste("PC2 : ", round(pca.var.per[2], digits = 4), "%", sep="")),
       zlab=(paste("PC3 : ", round(pca.var.per[3], digits = 4), "%", sep="")))


#import annotations:
ClinicalInfo <- read.delim(file = 'TCGA_BRCA_ClinicalAnnotation.txt.gz', header=TRUE,stringsAsFactors = FALSE)

#remove NA values, make them text
ClinicalInfo[is.na(ClinicalInfo)]<-"NA"
#make sure the no.of columns is the same in ReadCounts as no. of rows in ClinicalInfo

colnames_t <- str_sub(colnames_lib, start=1, end=-4) #disregards type of tumor
colnames_t<- gsub('\\.', '-', colnames_t)
colnames_index<-match(colnames_t,ClinicalInfo[,1])

new_clinical = ClinicalInfo[FALSE,]
i_set<- seq(1, length(colnames_t))
for (i in i_set) {
  new_clinical[i,]<-ClinicalInfo[colnames_index[i],]
}

#insert column in clinical annotations with tumor type
aa <- str_sub(colnames_lib, start=14)
for (i in i_set){
  if (aa[i]=="01") #tumor
  {aa[i]<- as.numeric(1) }
  else if (aa[i]=="06") #metasthasis
  {aa[i]<- as.numeric(2)}
  else {aa[i]<- as.numeric(0)} #healthy
}

#add type of tumor + gene coverage
new_clinical <- add_column(new_clinical, aa, .after = "Patient.ID")
new_clinical <- add_column(new_clinical, genereads$Total_counts, .after = "aa")


# get the name of the top 50 measurements (genes) that contribute most to pc1.
loading_scores <-PCA_basic$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_genes <- names(gene_score_ranked[1:20])

top_genes ## show the names of the top 10 genes

PCA_basic$rotation[top_genes,1] ## show the scores (and +/- sign)

#now, for the y axis
# get the name of the top 50 measurements (genes) that contribute most to pc2.
loading_scores <-PCA_basic$rotation[,2]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_genes <- names(gene_score_ranked[1:20])

top_genes ## show the names of the top 10 genes

PCA_basic$rotation[top_genes,2] ## show the scores (and +/- sign)

#Studying variation 
pca_data <- data.frame(PCA_basic$x)

#########################
#First study - Gene Coverage


label1 <-cut(genereads$Total_counts, c(0,6*10^7,1*10^8,1.4*10^8,1*10^9))

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Total No. of Reads")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Read Coverage")


dataplot1.3d <- data.frame(label1, X=pca_data[,1], Y=pca_data[,2], Z=pca_data[,3])
colnames(dataplot1.3d)<- c("label","x","y","z")

dev.new()
cols=c("red","blue","green","purple")
col = cols[as.numeric(label1)]
legend3d("topright", legend = levels(label1), pch=15,col=cols)
plot3d(x=dataplot1.3d$x, y=dataplot1.3d$y, z=dataplot1.3d$z, col=col,
       xlab = "PC1",
       ylab="PC2",
       zlab= "pc3") 

#This confirms normalization was well performed


#########################
#Second study - Gender

label1 <-new_clinical$Gender

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Gender")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Gender")

dataplot1.3d <- data.frame(label1, X=pca_data[,1], Y=pca_data[,2], Z=pca_data[,3])
colnames(dataplot1.3d)<- c("label","x","y","z")

mycolors <- c('royalblue1','pink')
dataplot1.3d$color <- mycolors[as.numeric(new_clinical$Gender)]

dev.new()
plot3d(x=dataplot1.3d$x, y=dataplot1.3d$y, z=dataplot1.3d$z, col= dataplot1.3d$color,
       xlab = "PC1",
       ylab="PC2",
       zlab= "pc3") #VER CORES!

#########################
#Third study - Age at diagnosis

label1 <- cut(new_clinical$Age.at.diagnosis..years., c(0,40,60,80,100))

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Age at Diagnosis")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Age at Diagnosis")

dataplot1.3d <- data.frame(label1, X=pca_data[,1], Y=pca_data[,2], Z=pca_data[,3])
colnames(dataplot1.3d)<- c("label","x","y","z")

mycolors <- c('royalblue1','yellow','orange','green','black')
dataplot1.3d$color <- mycolors[as.numeric(label1)]

dev.new()
plot3d(x=dataplot1.3d$x, y=dataplot1.3d$y, z=dataplot1.3d$z, col=dataplot1.3d$color,
       xlab = "PC1",
       ylab="PC2",
       zlab= "pc3") #VER CORES!

#########################
#Forth study - Ethnicity

label1 <-new_clinical$Ethnicity

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Ethnicity")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Ethnicity")

#########################
#Fifth study - Menopausal Status

label1 <-new_clinical$Menopausal.status


dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Menopausal Status")) +
  scale_color_manual(labels = c("Indeterminate", "NA","peri","post","pre"),values = c ("red","yellow","green","pink","blue"))+
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Menopausal Status")

#########################
#Sixth study - Anatomic Subdivision

label1 <-new_clinical$Anatomic.subdividion


dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Anatomic Subdivision")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Anatomic Subdivision")


#########################
#Seventh study - Tumor type

label1 <-new_clinical$aa


dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Type of Tumor")) +
  scale_color_manual(labels = c("Healthy","Primary", "Metasthasis"),values = c ("purple","pink","black"))+
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Type of Tumor")

#########################
#Eighth study - Surgical Procedure

label1 <-new_clinical$Surgical.procedure


dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Surgical Procedure")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Surgical Procedure")

#########################
#Nineth study - PAM50

label1 <-new_clinical$PAM50

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

dev.new()
ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("PAM50")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per PAM50")

#IId)

#building design matrix, with tumor type and age, considering intersept
design.matrix <- model.matrix(~aa*Age.at.diagnosis..years., data = new_clinical)

#remember that 0-healthy, 1-tumor, 2- metástase. First column is baseline (healthy), 
#2nd is tumor, etc.

dge2 <- DGEList(counts=data.matrix(ReadCounts))
dge2 <- calcNormFactors(dge2)
v2 <- voom(dge2,design.matrix,plot=TRUE)
#o lmfit faz-te o ajuste linear dos teus dados, 
#se lhe deres já as counts normalizadas pelo voom e a tua design matrix
linearfit = lmFit(v2$E,design.matrix) 
#o eBayes dá-te uma data de testes estatísticos para esse fit linear. 
#podes ver quais são se fizeres names(EBFit)
eBfit = eBayes(linearfit)

#TESTE##############################
summary(decideTests(eBfit))


################################
mypval=0.05
myfc=3

#tumorvsnormal
limma.res.pval <- topTable(eBfit,coef=2,sort.by = "t", number = Inf, p.val=mypval)
dim(limma.res.pval)

limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]
dim(limma.res.pval.FC)

#normalvstumorvsage
limmaage.res.pval <- topTable(eBfit,coef=5, sort.by = "t", number = Inf, p.val=mypval)
dim(limmaage.res.pval)

limmaage.res.pval.FC <- limmaage.res.pval[which(abs(limmaage.res.pval$logFC)>0.05),]
dim(limmaage.res.pval.FC)


#################################################

#o que nos interessa realmente é a função topTable, 
#que como o nome indica dá-te os top número que lá 
#puseres genes, de acordo com uma condição que eu quero testar

#por exemplo,  os genes diferencialmente expressos em tumores, 
#sabendo que a minha baseline é normal
HealthyTumor <- topTable(eBfit,coef=2,sort.by = "t", number = Inf)
HealthyMet <- topTable(eBfit,coef=3,sort.by = "t", number = Inf)
HealthyTumorAge <- topTable(eBfit,coef=5,sort.by = "t", number = Inf)
HealthyMetAge <- topTable(eBfit,coef=6,sort.by = "t", number = Inf)

##############
#Healthy vs Tumor
#############

#now for plotting --- a fancier version of this
#volcanoplot(eBfit,coef=2,style="B-statistic")

#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyTumor$adj.P.Val

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyTumor$logFC, HealthyTumor$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyTumor), plot_ht) 
limma.FC <- cbind(gene=rownames(limma.res.pval.FC), limma.res.pval.FC) 
interm_fc<-(match(plot_ht$gene,limma.FC$gene, nomatch = NA))
gene_labels <- data.frame(plot_ht$gene[teste])
plot_ht <- cbind(glabels = gene_labels, plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Volcano Plot, Healthy vs Tumor(Primary)") +
  geom_text_repel(data=head(plot_ht, 30), aes(label=plot_ht.gene.teste.)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 

################
#Healthy vs Metastasis
################


#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyMet$P.Value

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyMet$logFC, HealthyMet$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyMet), plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Volcano Plot, Healthy vs Metastasis") +
  geom_text_repel(data=head(plot_ht, 20), aes(label=gene)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 


################
#Healthy vs Tumor vs Age
################

#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyTumorAge$adj.P.Val

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyTumorAge$logFC, HealthyTumorAge$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyTumorAge), plot_ht) 
limma.FC <- cbind(gene=rownames(limmaage.res.pval.FC), limmaage.res.pval.FC) 
interm_fc<-(match(plot_ht$gene,limma.FC$gene, nomatch = NA))
gene_labels <- data.frame(limma.FC$gene[interm_fc])
plot_ht <- cbind(glabels = gene_labels, plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Volcano Plot, Healthy vs Tumor vs Age") +
  geom_text_repel(data=head(plot_ht, 98), aes(label=limma.FC.gene.interm_fc.)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 

################
#Healthy vs Met vs Age
################


#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyMetAge$P.Value

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}

significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyMetAge$logFC, HealthyMetAge$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyMetAge), plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Volcano Plot, Healthy vs Metastasis vs Age") +
  geom_text_repel(data=head(plot_ht, 20), aes(label=gene)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 

##############################################
#and now, join tumor and metastasis to see if results differ
########################################################

new_clinical2 = ClinicalInfo[FALSE,]
i_set<- seq(1, length(colnames_t))
for (i in i_set) {
  new_clinical2[i,]<-ClinicalInfo[colnames_index[i],]
}

#insert column in clinical annotations with tumor type
aa <- str_sub(colnames_lib, start=14)
for (i in i_set){
  if (aa[i]=="01") #tumor
  {aa[i]<- as.numeric(1) }
  else if (aa[i]=="06") #metasthasis
  {aa[i]<- as.numeric(1)}
  else {aa[i]<- as.numeric(0)} #healthy
}

#add type of tumor 
new_clinical2 <- add_column(new_clinical2, aa, .after = "Patient.ID")


#building design matrix, with tumor type and age, considering intersept
design.matrix <- model.matrix(~aa*Age.at.diagnosis..years., data = new_clinical2)

#remember that 0-healthy, 1-tumor. First column is baseline (healthy), 
#2nd is tumor, etc.

dge2 <- DGEList(counts=data.matrix(ReadCounts))
dge2 <- calcNormFactors(dge2)
v2 <- voom(dge2,design.matrix,plot=TRUE)
#o lmfit faz-te o ajuste linear dos teus dados, 
#se lhe deres já as counts normalizadas pelo voom e a tua design matrix
linearfit = lmFit(v2$E,design.matrix) 
#o eBayes dá-te uma data de testes estatísticos para esse fit linear. 
#podes ver quais são se fizeres names(EBFit)
eBfit = eBayes(linearfit)

#o que nos interessa realmente é a função topTable, 
#que como o nome indica dá-te os top número que lá 
#puseres genes, de acordo com uma condição que eu quero testar

#por exemplo,  os genes diferencialmente expressos em tumores, 
#sabendo que a minha baseline é normal
HealthyTumor2 <- topTable(eBfit,coef=2,sort.by = "t", number = Inf)
HealthyTumorAge2 <- topTable(eBfit,coef=4,sort.by = "t", number = Inf)


##############
#Healthy vs Tumor
#############

#now for plotting --- a fancier version of this
#volcanoplot(eBfit,coef=2,style="B-statistic")

#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyTumor2$P.Value

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyTumor2$logFC, HealthyTumor2$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyTumor2), plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Volcano Plot, Healthy vs Tumor(Primary and Metastasis)") +
  geom_text_repel(data=head(plot_ht, 20), aes(label=gene))#plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 

################
#Healthy vs Tumor vs Age
################


#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyTumorAge2$P.Value

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyTumorAge2$logFC, HealthyTumorAge2$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyTumorAge2), plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Volcano Plot, Healthy vs Tumor vs Age (Primary and Metastasis)") +
  geom_text_repel(data=head(plot_ht, 20), aes(label=gene)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 

################################

#HealthyTumor <- topTable(eBfit,coef=2,sort.by = "t", number = Inf)

GoAnalysis <- data.frame (HealthyTumor$t)
colnames(GoAnalysis)<- c("t")
GoAnalysis <- cbind(gene = rownames(HealthyTumor), GoAnalysis)

GoAnalysis_s <- sort(GoAnalysis$t, decreasing=TRUE)
GoAnalysis_s <-GoAnalysis[order(-GoAnalysis$t) ,]

write.table(GoAnalysis_s, "C:/Users/Carolina/Ambiente de Trabalho/BC/Lab 5 data/GEASe.txt",sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE) 


###############
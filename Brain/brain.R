## Make the Brain folder the working directory.

## install NbClust, cidr, tsne
library(NbClust)
library(cidr)
library(tsne)

## Read in data
load("Data/brain.RData")

info <- read.csv("Data/SraRunTable.txt",sep="\t")

cell_type <- info$cell_type_s[match(colnames(brain),info$Sample_Name_s)]

hyprid <- which(cell_type=="hybrid")

brain <- brain[,-hyprid]
cell_type <- info$cell_type_s[match(colnames(brain),info$Sample_Name_s)]
cell_type <- factor(cell_type)

types <- levels(cell_type)

# Set colours
col_2 <- rgb(0/255,73/255,73/255,1)
col_4 <- rgb(255/255,109/255,182/255,1)
col_6 <- rgb(73/255,0/255,146/255,1)
col_7 <- rgb(0/255,109/255,219/255,1)
col_8 <- rgb(182/255,109/255,255/255,1)
col_10 <- rgb(109/255,219/255,255/255,1)
col_11 <- rgb(146/255,0/255,0/255,1)
col_12 <- rgb(146/255,73/255,0/255,1)
col_13 <- rgb(219/255,209/255,0/255,1)
col_14 <- rgb(36/255,255/255,36/255,1)
col_15 <- rgb(255/255,255/255,109/255,1)

scols <- c(col_2, col_4, col_6, col_7, col_8, col_10, col_11, col_14)

adj_rand_ind_col <- "#e0e0e0" # light grey


cols <- rep(NA,length(cell_type))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cell_type[i])]
}

## Simple pre-processing
brain <- brain[rowSums(brain)>0,]

priorTPM <- 1
brain_lcpm <- log2(t(t(brain)/colSums(brain))*1000000+priorTPM)

brain10 <- brain[rowSums(brain)>10,]
brain10_lcpm <- log2(t(t(brain10)/colSums(brain10))*1000000+priorTPM)

## legend
types[3] <- "fetal quiescent neurons"
types[4] <- "fetal replicating neurons"
types[8] <- "oligodendrocyte precursor cells"

#### CIDR ######
################

start <- Sys.time()

scBrain <- scDataConstructor(as.matrix(brain))

scBrain <- determineDropoutCandidates(scBrain)
scBrain <- wThreshold(scBrain)
scBrain <- scDissim(scBrain)
scBrain <- scPCA(scBrain)
scBrain <- nPC(scBrain)
scBrain<- scCluster(scBrain)

Sys.time() - start

plot(scBrain@PC[,c(1,2)],col=cols,pch=scBrain@clusters,main="CIDR",xlab="PC1", ylab="PC2")
cidrM <- intersect(which(scBrain@clusters==4), which(cell_type!="neurons"))

text(scBrain@PC[cidrM,1], scBrain@PC[cidrM,2],label=1:6,pos=1)

nCluster(scBrain)
ARI_cidr <- adjustedRandIndex(scBrain@clusters,cols)
ARI_cidr
## 0.8977449

## CIDR_Full

start <- Sys.time()

scBrain_F <- scDataConstructor(as.matrix(brain))

scBrain_F <- determineDropoutCandidates(scBrain_F)
scBrain_F <- wThreshold(scBrain_F)
scBrain_F <- scDissim(scBrain_F, useStepFunction=FALSE)
scBrain_F <- scPCA(scBrain_F)
scBrain_F <- nPC(scBrain_F)
scBrain_F <- scCluster(scBrain_F)

Sys.time() - start

plot(scBrain_F@PC[,c(1,2)],col=cols,pch=scBrain_F@clusters,main="CIDR",xlab="PC1", ylab="PC2")

nCluster(scBrain_F)
ARI_cidr_F <- adjustedRandIndex(scBrain_F@clusters,cols)
ARI_cidr_F
## 0.8756049


##  prcomp
start <- Sys.time()
pca <- prcomp(t(brain10_lcpm))
variation <- summary(pca)[[6]][2,]
nPCs <- calc_npc(variation)

CH_prcomp <- NbClust(pca$x[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc=nPCs*3+3)$All.index
l <- length(CH_prcomp)
a <- as.vector(CH_prcomp[-c(1,l-1,l)]+CH_prcomp[-c(1:3)] - 2*CH_prcomp[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_prcomp <- b+c+2
} else {
  con_prcomp <- b+2
}

clusters_prcomp <- NbClust(pca$x[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_prcomp, max.nc=con_prcomp)$Best.partition

Sys.time()-start

ARI_prcomp <- adjustedRandIndex(clusters_prcomp,cols)
ARI_prcomp
# 0.4802283 

plot(pca$x[,c(1,2)],col=cols,pch=clusters_prcomp,xlab="PC1",ylab="PC2",main="prcomp")
prcompM <- intersect(which(clusters_prcomp%in%c(4,6,8)),which(cell_type!="neurons"))
prcompM
text(pca$x[prcompM,1],pca$x[prcompM,2],label=1:2,pos=1)

#############################################
#############################################

## t-SNE

nPCs <- scBrain@nPC

library(tsne)
set.seed(20)
start<-Sys.time()
y_tsne2 <- tsne(dist(t(brain10_lcpm)), k=nPCs, perplexity=10)
CH_tsne <- NbClust(y_tsne2[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc = nPCs*3+3)$All.index
l <- length(CH_tsne)
a <- as.vector(CH_tsne[-c(1,l-1,l)]+CH_tsne[-c(1:3)] - 2*CH_tsne[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_tsne <- b+c+2
} else {
  con_tsne <- b+2
}

clusters_tsne<- NbClust(y_tsne2[,c(1:nPCs)],method="ward.D2",index="ch",min.nc=con_tsne,max.nc = con_tsne)$Best.partition

Sys.time()-start

# plot
plot(y_tsne2[,c(1,2)],col=cols,pch=clusters_tsne,main="t-SNE",xlab="PC1", ylab="PC2")
tsneM <- c(intersect(which(clusters_tsne==2), which(cell_type!="neurons")))
tsneM[c(15,31,32)]
tsneM <- c(62,133,136)
text(y_tsne2[tsneM,1], y_tsne2[tsneM,2], label=1:3, pos=c(1,2,3))
ARI_tsne <- adjustedRandIndex(clusters_tsne,cols)
ARI_tsne
##  0.5692866

## ZIFA ##################
## forZifa <- t(log2(brain10+1))
## write.csv(forZifa,file="Results/forZifa.csv",row.names = FALSE)
## run zifa_brain.py
## 1:08:43.444764
## 68.7 mins
## 1.1 hours

Z4 <- read.csv(file="Results/Z4.csv", header = FALSE)

CH_zifa <- NbClust(Z4[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc=nPCs*3+3)$All.index
l <- length(CH_zifa)
a <- as.vector(CH_zifa[-c(1,l-1,l)]+CH_zifa[-c(1:3)] - 2*CH_zifa[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_zifa <- b+c+2
} else {
  con_zifa <- b+2
}
clusters_zifa<- NbClust(Z4[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_zifa, max.nc=con_zifa)$Best.partition

plot(Z4[,c(1,2)],col=cols,pch=clusters_zifa,main="ZIFA",xlab="PC1", ylab="PC2")
zifaM <- intersect(which(clusters_zifa%in%c(4,6)), which(cell_type!="neurons"))
zifaM
text(Z4[zifaM,1], Z4[zifaM,2], label=1:2, pos=4)


ARI_zifa <- adjustedRandIndex(clusters_zifa,cols)
ARI_zifa
##  0.5349539

## file.remove("Results/forZifa.csv")

#################
## RaceID #######
#################

source("RaceID_class.R")

set.seed(5)
start <- Sys.time()
# initialize SCseq object with transcript counts
sc <- SCseq(brain10)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000)
# k-means clustering
sc <- clustexp(sc, clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,rseed=17000)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75)
##
Sys.time()-start
## Time difference of 1.503236 mins
## 90.2 secs

######
clusters_raceid <- sc@kmeans$kpart
plot(y_tsne2[,c(1,2)],col=cols,pch=clusters_raceid,main="RaceID",xlab="PC1",ylab="PC2")

raceidM <- intersect(which(clusters_raceid%in%c(2, 3)),which(cell_type!="neurons"))
raceidM[c(32,24,43)]
## [1] 216 190 290
raceidM <- c(216, 190, 290)
text(y_tsne2[raceidM,1], y_tsne2[raceidM,2],label=1:3,pos=c(4,4,1))

ARI_raceID <- adjustedRandIndex(sc@kmeans$kpart,cols)
ARI_raceID
## 0.3854667

barplot(c(ARI_prcomp,ARI_tsne,ARI_zifa,ARI_raceID,ARI_cidr),col="BLUE",
        names=c("prcomp", "t-SNE","ZIFA","RaceID","CIDR"),
        main = "Adjusted Rand Index",ylim=c(0,1))

clusters_raceid <- sc@kmeans$kpart
clusters_cidr <- scBrain@clusters

clusters_brain <- list(clusters_prcomp, clusters_tsne, clusters_zifa, 
                    clusters_raceid, clusters_cidr)
names(clusters_brain) <- c("prcomp","tsne","zifa","raceid","cidr")


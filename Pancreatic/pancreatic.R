## Make the Pancreatic folder the working directory

## Install packages NbClust, cidr, tsne
library(NbClust)
library(cidr)
library(tsne)

## Read in data and filtering

load("Data/counts.RData")
pan <- counts - 1

info <- read.csv("Data/SraRunTable.txt",sep="\t")

pan <- pan[,!is.na(info$assigned_cell_type_s[match(colnames(pan),info$Sample_Name_s)])]
pan <- pan[,(info$assigned_cell_type_s[match(colnames(pan),info$Sample_Name_s)])!="undefined"]
pan <- pan[,(info$assigned_cell_type_s[match(colnames(pan),info$Sample_Name_s)])!="islet"]

cell_type <- info$assigned_cell_type_s[match(colnames(pan),info$Sample_Name_s)]

cell_type <- factor(cell_type)

types <- levels(cell_type)

## Set colours
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
# currently 6 types - adjust if types change
scols <- c(col_2, col_4, col_6, col_7, col_11, col_14)
adj_rand_ind_col <- "#e0e0e0" # light grey
cols <- rep(NA,length(cell_type))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cell_type[i])]
}

## simple pre-processing
pan <- pan[rowSums(pan)>0,]

priorTPM <- 1
pan_lcpm <- log2(t(t(pan)/colSums(pan))*1000000+priorTPM)

pan10 <- pan[rowSums(pan)>10,]
pan10_lcpm <- log2(t(t(pan10)/colSums(pan10))*1000000+priorTPM)

###########
## CIDR ###
###########

start <- Sys.time()

scPan <- scDataConstructor(as.matrix(pan))
scPan <- determineDropoutCandidates(scPan)
scPan <- wThreshold(scPan)
scPan <- scDissim(scPan)
scPan <- scPCA(scPan)
scPan <- nPC(scPan)
scPan <- scCluster(scPan)

Sys.time() - start

nCluster(scPan)
ARI_CIDR <- adjustedRandIndex(scPan@clusters,cols)
ARI_CIDR
##  0.6830087 

plot(scPan@PC[,c(1,2)],col=cols,pch=scPan@clusters,main="CIDR",xlab="PC1", ylab="PC2")


###########
## CIDR_FULL ###
###########

start <- Sys.time()

scPan_F <- scDataConstructor(as.matrix(pan))
scPan_F <- determineDropoutCandidates(scPan_F)
scPan_F <- wThreshold(scPan_F)
scPan_F <- scDissim(scPan_F, useStepFunction=FALSE)
scPan_F <- scPCA(scPan_F)
scPan_F <- nPC(scPan_F)
scPan_F <- scCluster(scPan_F)

Sys.time() - start

nCluster(scPan_F)
ARI_CIDR_F <- adjustedRandIndex(scPan_F@clusters,cols)
ARI_CIDR_F
##  0.4161589

# prcomp
start <- Sys.time()
PC_lcpm <- prcomp(t(pan10_lcpm))
variation <- summary(PC_lcpm)[[6]][2,]
nPCs <- calc_npc(variation)
CH_prcomp <- NbClust(PC_lcpm$x[,c(1:nPCs)],method="ward.D2",index="ch",min.nc=1,max.nc=nPCs*3+3)$All.index
l <- length(CH_prcomp)
a <- as.vector(CH_prcomp[-c(1,l-1,l)]+CH_prcomp[-c(1:3)] - 2*CH_prcomp[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_prcomp <- b+c+2
} else {
  con_prcomp <- b+2
}
clusters_prcomp <- NbClust(PC_lcpm$x[,c(1:nPCs)],method="ward.D2",index="ch",min.nc=con_prcomp,max.nc = con_prcomp)$Best.partition

Sys.time()-start

ARI_prcomp <- adjustedRandIndex(clusters_prcomp,cols)
ARI_prcomp
# 0.2081782

plot(PC_lcpm$x[,c(1,2)],col=cols,pch=clusters_prcomp,xlab="PC1",ylab="PC2",main="prcomp")

## t-SNE
nPCs <- scPan@nPC

set.seed(10)
start<-Sys.time()
y_tsne <- tsne(dist(t(pan10_lcpm)), k=nPCs, perplexity=10)
CH_tsne <- NbClust(y_tsne[,c(1:4)], method="ward.D2", index="ch", min.nc=1, max.nc=nPCs*3+3)$All.index
l <- length(CH_tsne)
a <- as.vector(CH_tsne[-c(1,l-1,l)]+CH_tsne[-c(1:3)] - 2*CH_tsne[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_tsne <- b+c+2
} else {
  con_tsne <- b+2
}
clusters_tsne<- NbClust(y_tsne[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_tsne, max.nc=con_tsne)$Best.partition

Sys.time()-start

plot(y_tsne[,c(1,2)],col=cols,pch=clusters_tsne,main="t-SNE",xlab="PC1", ylab="PC2")

ARI_tsne <- adjustedRandIndex(clusters_tsne,cols)
ARI_tsne
## [1] 0.1959748


## ZIFA ##################
## forZifa <- t(log2(pan10+1))
## write.csv(forZifa,file="Results/forZifa.csv",row.names = FALSE)
## run zifa_pan.py

## 0:40:08.280524
## 40.1 mins

Z4 <- read.csv(file="Results/Z4.csv",header = FALSE)

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
clusters_zifa<- NbClust(Z4[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_zifa, max.nc = con_zifa)$Best.partition

plot(Z4[,c(1,2)],col=cols,pch=clusters_zifa,main="ZIFA",xlab="PC1", ylab="PC2")

ARI_zifa <- adjustedRandIndex(clusters_zifa,cols)
ARI_zifa
# 0.2046189

## file.remove("/Results/forZifa.csv")

#################
## RaceID #######
#################

source("RaceID_class.R")

set.seed(5)
start <- Sys.time()
# initialize SCseq object with transcript counts
sc <- SCseq(pan10)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000)
# k-means clustering
sc <- clustexp(sc, clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,rseed=17000)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75)
##
Sys.time()-start

plot(y_tsne[,c(1,2)],col=cols,pch=sc@kmeans$kpart,main="RaceID",xlab="PC1",ylab="PC2")

ARI_raceID <- adjustedRandIndex(sc@kmeans$kpart,cols)
ARI_raceID
## 0.2220675

barplot(c(ARI_prcomp,ARI_tsne,ARI_zifa,ARI_raceID,ARI_CIDR),col="BLUE",
        names=c("prcomp","t-SNE","ZIFA","RaceID","CIDR"),
        main = "Adjusted Rand Index", ylim=c(0,0.7))


## Make the Simulation folder the working directory.

## Install packages NbClust, cidr, minpack.lm, tsne

library(NbClust)
library(cidr)
library(minpack.lm)
library(tsne)

## simulate data
k <- 50
N<-3

sData <- scSimulator(k=k,N=N)

x <- sData$tags
rownames(x) <- 1:nrow(x)
colnames(x) <- 1:ncol(x)
tags <-data.frame(x)

## Colours

cols <- rep(NA,k*N)

# http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
# chose colors 2, 6, 7, 8, 11, 12, 13, 14 (use 9 for barplot)
col_2 <- rgb(0/255,73/255,73/255,1)
col_4 <- rgb(255/255,109/255,182/255,1)
col_6 <- rgb(73/255,0/255,146/255,1)
col_7 <- rgb(0/255,109/255,219/255,1)
col_8 <- rgb(182/255,109/255,255/255,1)
col_9 <- rgb(109/255,182/255,255/255,1)
col_11 <- rgb(146/255,0/255,0/255,1)
col_12 <- rgb(146/255,73/255,0/255,1)
col_13 <- rgb(219/255,209/255,0/255,1)
col_14 <- rgb(36/255,255/255,36/255,1)

adj_rand_ind_col <- "#e0e0e0"  # light grey

col <- c(col_2, col_11, col_6)

sim_cols <- col
for(i in 1:N){
  cols[((i-1)*k+1):(i*k)] <- rep(col[i],k)
}

###########
## CIDR ###
###########

start <- Sys.time()
scData <- scDataConstructor(as.matrix(sData$tags))

scData <- determineDropoutCandidates(scData)
scData <- wThreshold(scData)
scData <- scDissim(scData)
scData <- scPCA(scData)
scData <- nPC(scData)
scData <- scCluster(scData)

Sys.time()-start
## Time difference of 1.934089 secs

plot(scData@PC[,c(1,2)],col=cols,pch=scData@clusters,main="CIDR",xlab="PC1", ylab="PC2")

nCluster(scData)
scData@nCluster
ARI_CIDR <- adjustedRandIndex(scData@clusters, cols)
ARI_CIDR
## 0.9203693

###########
## CIDR_FULL ###
###########

start <- Sys.time()
scData_F <- scDataConstructor(as.matrix(sData$tags))

scData_F <- determineDropoutCandidates(scData_F)
scData_F <- wThreshold(scData_F)
scData_F <- scDissim(scData_F, useStepFunction=FALSE)
scData_F <- scPCA(scData_F)
scData_F <- nPC(scData_F)
scData_F <- scCluster(scData_F)

Sys.time()-start
## Time difference of 2.294245 secs

plot(scData_F@PC[,c(1,2)],col=cols,pch=scData_F@clusters,main="CIDR",xlab="PC1", ylab="PC2")

nCluster(scData_F)
scData_F@nCluster
ARI_CIDR_F <- adjustedRandIndex(scData_F@clusters, cols)
ARI_CIDR_F
## 0.9013364


##  prcomp
start <- Sys.time()
ltpm_c <- log2(t(t(sData$tags)/colSums(sData$tags))*1000000+1)
y_pca_dn <- prcomp(t(ltpm_c)) # y-value for PCA - Drop-outs & Noise
variation <- summary(y_pca_dn)[[6]][2,]
nPCs <- calc_npc(var=variation)
CH_prcomp <- NbClust(y_pca_dn$x[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc=nPCs*3+3)$All.index
l <- length(CH_prcomp)
a <- as.vector(CH_prcomp[-c(1,l-1,l)]+CH_prcomp[-c(1:3)] - 2*CH_prcomp[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_prcomp <- b+c+2
} else {
  con_prcomp <- b+2
}
clusters_prcomp <- NbClust(y_pca_dn$x[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_prcomp, max.nc=con_prcomp)$Best.partition
Sys.time()-start
# Time difference of 2.885882 secs

ARI_prcomp <- adjustedRandIndex(clusters_prcomp,cols)
ARI_prcomp
#  0.4767767

plot(y_pca_dn$x[,c(1,2)],col=cols,pch=clusters_prcomp,xlab="PC1",ylab="PC2",main="prcomp")

## t-SNE
nPCs <- scData@nPC

set.seed(20)
start <- Sys.time()
ltpm_c <- log2(t(t(sData$tags)/colSums(sData$tags))*1000000+1)
y_tsne2 <- tsne(dist(t(ltpm_c)), k=nPCs, perplexity=10)
CH_tsne <- NbClust(y_tsne2[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc = 3*nPCs+3)$All.index
l <- length(CH_tsne)
a <- as.vector(CH_tsne[-c(1,l-1,l)]+CH_tsne[-c(1:3)] - 2*CH_tsne[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_tsne <- b+c+2
} else {
  con_tsne <- b+2
}
clusters_tsne<- NbClust(y_tsne2[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_tsne, max.nc = con_tsne)$Best.partition
Sys.time()-start
# Time difference of 14.24438 secs

plot(y_tsne2[,c(1,2)],col=cols,pch=clusters_tsne,main="t-SNE",xlab="PC1", ylab="PC2")

ARI_tsne <- adjustedRandIndex(clusters_tsne,cols)
ARI_tsne
##  [1] 0.01596457

## ZIFA ##################
## forZifa <- t(log2(sData$tags+1))
## write.csv(forZifa,file="results/forZifa.csv",row.names = FALSE)

## now run zifa_simulation.py
# 0:32:03.046335
# 32.1 mins

Z2 <- read.csv(file="results/Z2.csv",header = FALSE)

CH_zifa <- NbClust(Z2[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc=nPCs*3+3)$All.index
l <- length(CH_zifa)
a <- as.vector(CH_zifa[-c(1,l-1,l)]+CH_zifa[-c(1:3)] - 2*CH_zifa[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_zifa <- b+c+2
} else {
  con_zifa <- b+2
}
clusters_zifa<- NbClust(Z2[,c(1:nPCs)],method="ward.D2",index="ch",min.nc=con_zifa,max.nc = con_zifa)$Best.partition
ARI_ZIFA <- adjustedRandIndex(clusters_zifa,cols)
ARI_ZIFA
## -0.002382703

## file.remove("results/forZifa.csv")

plot(Z2[,c(1,2)], col=cols, pch=clusters_zifa, main="ZIFA", xlab="PC1", ylab="PC2")

## RaceID
## install required packages (only at first time)
##install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit"))

source("RaceID_class.R")

set.seed(50)

start <- Sys.time()
# initialize SCseq object with transcript tags
sc <- SCseq(tags)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000)
# k-means clustering
sc <- clustexp(sc, clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,rseed=17000)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75)
Sys.time()-start
## Time difference of 20.6664 secs

######
plot(y_tsne2[,c(1,2)],col=cols,pch=sc@kmeans$kpart,main="RaceID",xlab="PC1",ylab="PC2")

ARI_raceID <- adjustedRandIndex(sc@kmeans$kpart,cols)
ARI_raceID
## 0

barplot(c(ARI_prcomp,ARI_tsne,ARI_ZIFA,ARI_raceID,ARI_CIDR),col="BLUE",
        names=c("prcomp","t-SNE","ZIFA","RaceID","CIDR"),
        main = "Adjusted Rand Index",ylim=c(0,1))


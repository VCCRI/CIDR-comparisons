## Make the folder brain_mouse the working directory. 

## Install packages NbClust, cidr, tsne.
library(NbClust)
library(cidr)
library(tsne)

## Read in data 

load("Data/brain_mouse.RData")
load("Data/info.RData")

rownames(info) <- info$tissue
info <- info[, -c(1:2)]
ct <- lapply(info[8,], FUN="as.character") 
cts <- rep(NA, length(ct))
for(i in 1:length(ct)){
    cts[i] <- ct[[i]]
}
id <- lapply(info[7,], FUN="as.character")
ids <- rep(NA, length(id))
for (i in 1:length(id)){
    ids[i] <- id[[i]]
}

cell_types <-cts[match(colnames(brain_mouse),ids)]

cell_types <- factor(cell_types)
 
types <- levels(cell_types)

# colors

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

scols <- c(col_2, col_4, col_6, col_7, col_8, col_10, col_11)

adj_rand_ind_col <- "#e0e0e0" # light grey

cols <- rep(NA,length(cell_types))
for (i in 1:length(cols)){
   cols[i] <- scols[which(types==cell_types[i])]
}

## simple pre-processing
brain <- brain_mouse[rowSums(brain_mouse)>0,]

priorTPM <- 1
brain_lcpm <- log2(t(t(brain)/colSums(brain))*1000000+priorTPM)

brain10 <- brain[rowSums(brain)>10,]
brain10_lcpm <- log2(t(t(brain10)/colSums(brain10))*1000000+priorTPM)

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
## Time difference of 57.94945 secs

plot(scBrain@PC[,c(1,2)],pch=scBrain@clusters,main="CIDR",xlab="PC1", ylab="PC2",col=cols)

nCluster(scBrain)

ARI_cidr <- adjustedRandIndex(scBrain@clusters, cell_types)
ARI_cidr
# 0.5226713
clusters_cidr <- scBrain@clusters

#### CIDR - Full ######
################

start <- Sys.time()

scBrain_F <- scDataConstructor(as.matrix(brain))

scBrain_F <- determineDropoutCandidates(scBrain_F)
scBrain_F <- wThreshold(scBrain_F)
scBrain_F <- scDissim(scBrain_F, useStepFunction=FALSE)
scBrain_F <- scPCA(scBrain_F)
scBrain_F <- nPC(scBrain_F)
scBrain_F <- scCluster(scBrain_F)

Sys.time() - start
## Time difference of 1.053636 mins

plot(scBrain_F@PC[,c(1,2)],pch=scBrain_F@clusters,main="CIDR",xlab="PC1", ylab="PC2",col=cols)

nCluster(scBrain_F)

ARI_cidr_F <- adjustedRandIndex(scBrain_F@clusters, cell_types)
ARI_cidr_F
# 0.3705627

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
clusters_prcomp <- NbClust(pca$x[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_prcomp, max.nc = con_prcomp)$Best.partition

Sys.time()-start

## Time difference of 3.214987 mins

ARI_prcomp <- adjustedRandIndex(clusters_prcomp, cell_types)
ARI_prcomp
# [1] 0.2626759

plot(pca$x[,c(1,2)],pch=clusters_prcomp,xlab="PC1",ylab="PC2",main="prcomp", col=cols)


## t-SNE
nPCs <- scBrain@nPC

library(tsne)
set.seed(10)
start<-Sys.time()
y_tsne <- tsne(dist(t(brain10_lcpm)), k=nPCs, perplexity=10)
## save(y_tsne, file="Results/y_tsne.RData")

load("Results/y_tsne.RData")
CH_tsne <- NbClust(y_tsne[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc=nPCs*3+3)$All.index
l <- length(CH_tsne)
a <- as.vector(CH_tsne[-c(1,l-1,l)]+CH_tsne[-c(1:3)] - 2*CH_tsne[-c(1,2,l)])
b <- which.min(a)
c <- which.min(a[-c(1:b)])
if ((3*a[b+c])<a[b]){
  con_tsne <- b+c+2
} else {
  con_tsne <- b+2
}
clusters_tsne<- NbClust(y_tsne[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=con_tsne, max.nc = con_tsne)$Best.partition

Sys.time()-start
## Time difference of 23.07392 mins

plot(y_tsne[,c(1,2)],pch=clusters_tsne,main="t-SNE",xlab="PC1", ylab="PC2",col=cols)

ARI_tsne <- adjustedRandIndex(clusters_tsne, cell_types)
ARI_tsne
## 0.6218051

## ZIFA ##################
## forZifa <- t(log2(brain10+1))
## write.csv(forZifa,file="Results/forZifa.csv",row.names = FALSE)
## run zifa_brain.py
## 1:50:34.885381
##  110.5814 mins
## 1.8 hours

Z4 <- read.csv(file="Results/Z4.csv", header = FALSE)

CH_zifa <- NbClust(Z4[,c(1:nPCs)], method="ward.D2", index="ch", min.nc=1, max.nc =nPCs*3+3)$All.index
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

plot(Z4[,c(1,2)],pch=clusters_zifa,main="ZIFA",xlab="PC1", ylab="PC2", col=cols)

ARI_zifa <- adjustedRandIndex(clusters_zifa,cell_types)
ARI_zifa
## 0.3228093

## file.remove("/Users/paulyLin/Dropbox/CIDR_paper/manuscript/Revision/Mouse_Brain/Results/forZifa.csv")

#################
## RaceID #######
#################

source("RaceID_class.R")

set.seed(5)
start <- Sys.time()
start
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

## Time difference of 2.503459 hours
## 150.2075 mins

## Warning messages:
##1: did not converge in 10 iterations 
##2: did not converge in 10 iterations 
##3: did not converge in 10 iterations 
##4: did not converge in 10 iterations 

######
load("Results/clusters_raceid.RData")
plot(y_tsne[,c(1,2)],pch=clusters_raceid,main="RaceID",xlab="PC1",ylab="PC2", col=cols)

ARI_raceid <- adjustedRandIndex(clusters_raceid,cell_types)
ARI_raceid
##  0.3683523

barplot(c(ARI_prcomp,ARI_tsne,ARI_zifa,ARI_raceid,ARI_cidr),col="BLUE",
        names=c("prcomp", "t-SNE","ZIFA","RaceID","CIDR"),
        main = "Adjusted Rand Index",ylim=c(0,1))

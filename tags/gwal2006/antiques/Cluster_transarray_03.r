## 2003-02-04
## Dirk.Repsilber@gmx.de
##
##
## call by:
##
## echo 'number.environments <- 1; number.timepoints <- 3; filename <- "test.dat"' | R --slave --save
## R --slave --restore < Cluster_transarray_03.r


## loading the clustering functions:
library(mva)

## function for getting the names of genes and cluster-tree-nodes:
name <- function(a){
  if(a<0)  value <- names[-a]
  if(a>=0) value <- toString(a)
  value
}

# length of single profile vector:
N <- number.environments*number.timepoints

## direct all plots to a nice ps-file:
## postscript("ClusteResults.ps")

## reading array data:
daten <- read.table(file=filename,header=FALSE,skip=1)
names <- as.vector(daten[,1])
daten <- as.matrix(daten[,(1+number.timepoints+1):(N+number.timepoints+1)])
no.messbar <- length(names)

## prepare graphics output:
par(mfrow=c(2,1))

## Pearson correlation distances:---------------------------------------------------
pears.sim <- matrix(nrow=no.messbar, ncol=no.messbar)
Phi <- c()
for(gene in 1:no.messbar) Phi[gene] <- sqrt(sum(daten[gene,]^2/N))
for(gene1 in 1:no.messbar)
  {
    for(gene2 in 1:gene1)
      {
        pears.sim[gene1,gene2] <- abs(1/N*sum(daten[gene1,]*daten[gene2,]/Phi[gene1]/Phi[gene2]))
      }
  }
p.dist <- as.dist(1-pears.sim)

## Clustering:
clust.res.pears <- hclust(p.dist, method = "average")
plot(clust.res.pears, hang=-1, labels=names, main=paste(filename,' : PEARS',sep=''))

## Output merge table with distances:
out.pears <- cbind(clust.res.pears$merge,clust.res.pears$height)
write(paste('results pearson:') ,file="")
write(t(apply(out.pears,c(1,2),name)), file="", ncolumns=3)

## Eukleadian distances:-----------------------------------------------------------
p.dist <- dist(daten,method = "euclidean",upper=TRUE)

## Clustering:
clust.res.eukl <- hclust(p.dist, method = "average")
plot(clust.res.eukl, hang=-1, labels=names, main=paste(filename,' : EUKL',sep=''))

## Output merge table with distances: 
out.eukl <- cbind(clust.res.eukl$merge,clust.res.eukl$height)
write(paste('\nresults euklid:'),file="")
write(t(apply(out.eukl,c(1,2),name)),file="",ncolumns=3)
            


## turn off graphics device:
##dev.off()











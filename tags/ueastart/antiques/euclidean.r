# stuff1d.dat
# number.environments <- 1
# number.timepoints <- 1
# filename <- "stuff1d.dat"
#
# stuff2d.dat 2 x 1
number.environments <- 2
number.timepoints <- 1
filename <- "stuff2d.dat"
#
# stuff2d.dat 1 x 2
# number.environments <- 1
# number.timepoints <- 2
# filename <- "stuff2d.dat"
#
# stuff3d.dat 1 x 3
# number.environments <- 1
# number.timepoints <- 3
# filename <- "stuff3d.dat"
#
# stuff4d.dat 2 x 2
# number.environments <- 2
# number.timepoints <- 2
# filename <- "stuff4d.dat"
#
# stuff4d.dat 1 x 4
# number.environments <- 1
# number.timepoints <- 4
# filename <- "stuff4d.dat"
#
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
# JTK: I think this is wrong...
# dmatrix <- as.matrix(daten[,(1+number.timepoints+1):(N+number.timepoints+1)])
# ... and this one is right. Dirk: please check
dmatrix <- as.matrix(daten[,1:N + 1])
no.messbar <- length(names)

## prepare graphics output:
par(mfrow=c(2,1))

## Pearson correlation distances:---------------------------------------------------
## JTK: thid doesn't work properly -- NaN values arise... so, commented out...
## pears.sim <- matrix(nrow=no.messbar, ncol=no.messbar)
## Phi <- c()
## for(gene in 1:no.messbar) Phi[gene] <- sqrt(sum(dmatrix[gene,]^2/N))
## for(gene1 in 1:no.messbar)
##   {
##     for(gene2 in 1:gene1)
##       {
##         pears.sim[gene1,gene2] <- abs(1/N*sum(dmatrix[gene1,]*dmatrix[gene2,]/Phi[gene1]/Phi[gene2]))
##       }
##   }
## p.dist <- as.dist(1-pears.sim)

## Clustering:
## clust.res.pears <- hclust(p.dist, method = "average")
## plot(clust.res.pears, hang=-1, labels=names, main=paste(filename,' : PEARS',sep=''))

## Output merge table with distances:
## out.pears <- cbind(clust.res.pears$merge,clust.res.pears$height)
## write(paste('results pearson:') ,file="")
## write(t(apply(out.pears,c(1,2),name)), file="", ncolumns=3)
## write(paste(''),file="")

## Eukleadian distances:-----------------------------------------------------------
p.dist <- dist(dmatrix,method = "euclidean",upper=TRUE)
p.dist
## Clustering:
clust.res.eukl <- hclust(p.dist, method = "average")
plot(clust.res.eukl, hang=-1, labels=names, main=paste(filename,' : EUKL',sep=''))

## Output merge table with distances: 
out.eukl <- cbind(clust.res.eukl$merge,clust.res.eukl$height)
## write(paste('results euklid:'),file="")
write(t(apply(out.eukl,c(1,2),name)),file="",ncolumns=3)
write(paste(''),file="")
            


## turn off graphics device:
##dev.off()"""

# R file with functioons suitable for conducting several analysis of the
# transsys simulator's results.

readSimData <- function(inputfile)
{
  d <- read.table(inputfile, header=TRUE);
  return(d)
}

readDistData <- function(inputfile)
{
  d <- read.table(inputfile, header=TRUE);
  return(d)
}

plotDistance <- function(distanceList)
{
  plot(distanceList[["eu_dist"]], type='l', xlab="Timesteps", ylab="Euclidean Distance")
}


# add the SVD-PCA analysis as i has been implemented before....


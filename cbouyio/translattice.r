# The .R plotting package for the lattice factor table.
# cbouyio, UEA, 28/09/2006

readTransLattice <- function(filename)
{
  x <- read.table(filename, header = TRUE);
  return(x);
}


# Returns a lattice data frame containing only the
# rows of the time steps specified by the timesteps parameter.

getTimeSlice <- function(latticeFrame, timesteps)
{
  # Take a subset of the data.frame where the timestep is specified.
  timestepFrame <- subset(latticeFrame, timestep==timesteps);
  return(timestepFrame);
}


# Return a vector containing the permissible values of
# x coordinates for the lattice frame

getXCoordinates <- function(latticeFrame)
{
  xCoords <- latticeFrame[["x"]];
  return(unique(xCoords));
}


getYCoordinates <- function(latticeFrame)
{
  yCoords <- latticeFrame[["y"]];
  return(unique(yCoords));
}


getTimeSteps <- function(latticeFrame)
{
  return(unique(latticeFrame[["timestep"]]));
}


getMaxTimestep <- function(latticeFrame)
{
  return(max(getTimeSteps(latticeFrame)));
}

getConcentrationRange <- function(latticeFrame, factorName)
{
  cr <- c(min(latticeFrame[[factorName]]), max(latticeFrame[[factorName]]));
  return(cr);
}

getMaximumConcentration <- function(latticeFrame, factorName)
{
  mc <- max(latticeFrame[[factorName]])
  return(mc)
}

getXSize <- function(latticeFrame)
{
  return(max(getXCoordinates(latticeFrame)));
}


getYSize <- function(latticeFrame)
{
  return(max(getYCoordinates(latticeFrame)));
}


# Construct a matrix containing the concentration values
# of the factor specified by factorName at timestep.

getFactorConcentrationMatrix <- function(latticeFrame, factorName, timestep)
{
  if (length(timestep) != 1)
  {
    stop("multiple timesteps specified");
  }
  # Populate a matrix of the disered dimensions with NAs.
  timestepLattice <- getTimeSlice(latticeFrame, timestep);
  m <- matrix(NA, ncol=getYSize(timestepLattice), nrow=getXSize(timestepLattice));
  # Populate the matrix with the factor concentrations.
  for (i in 1:nrow(timestepLattice))
  {
    x <- timestepLattice[["x"]][i];
    y <- timestepLattice[["y"]][i];
    m[x, y] <- timestepLattice[[factorName]][i];
  }
  # Check for the existance of NA's in the matrix. (validation check)
  if (any(is.na(m)))
  {
    stop("Error in factor_concentration matrix population");
  }
  return(m);
}


plotConcentrationMatrix <- function(concentrationMatrix, xCoordinates, yCoordinates, concentrationRange, main="Image Title", ...)
{
  # All the job is done by the image function.
  image(xCoordinates, yCoordinates, concentrationMatrix, zlim = concentrationRange, xlab="x Coordinate", ylab="y Coordinate", col = heat.colors(64), main=main, ...);
  legend("right", c(as.character(round(seq(concentrationRange[1], concentrationRange[2], length.out=16), digits=3))), fill=heat.colors(16), title="Gradient", ...);
#  grid(max(xCoordinates), max(yCoordinates));
}


plotConcentrationSeries <- function(latticeFrame, factorName, concentrationRange=c(0, getMaximumConcentration(latticeFrame, factorName)), timeframeEndFunction=hitReturn, ...)
{
  # Some check about number of timesteps should precede.
  for (i in getTimeSteps(latticeFrame))
  {
    plotConcentrationMatrix(getFactorConcentrationMatrix(latticeFrame, factorName, i), getXCoordinates(latticeFrame), getYCoordinates(latticeFrame), concentrationRange, main=sprintf("Image of %s concentration on timestep %d", factorName, as.integer(i)), ...);
    timeframeEndFunction(i);
  }
}


plotAllInstances <- function(dataFrame, factorName, concentrationRange=c(0, max(dataFrame[[factorName]])), ylim=concentrationRange, ...)
{
  # First make the template plot using the first instance.
  instance1 <- subset(dataFrame, x==1 & y==1);
  factor1 <- instance1[[factorName]];
  plot(getTimeSteps(dataFrame), factor1, type="l", ylim=ylim, xlab="Timesteps", ylab="Factor Concentration", main=sprintf("All cell trajectories of factor %s", factorName), ...);
  # Then draw the lines.
  for (i in 1:getXSize(dataFrame))
  {
    for (j in 1:getYSize(dataFrame))
    {
      instance <- subset(dataFrame, x==i & y==j);
      factor <- instance[[factorName]];
      lines(getTimeSteps(dataFrame), factor, ...);
    }
  }
}


oneSecondDelay <- function(timestep)
{
  print(sprintf("timestep %d", as.integer(timestep)));
  Sys.sleep(1);
}


hitReturn <- function(timestep)
{
  readline(sprintf("timestep %d -- hit return", as.integer(timestep)));
}


getManhattanDistance <- function(x1, y1, x2, y2, X, Y)
# Calculate the Manhattan Distance between two cells on the toroidal lattice.
{
  manhX <- min(abs(x2 - x1), abs(x1 + X - x2))
  manhY <- min(abs(y2 - y1), abs(y1 + Y - y2))
  manhattanDistance = manhX + manhY
  return(manhattanDistance)
}


spatialCorrelation <- function(dframe, factorName, timestep=getMaxTimestep(dframe))
# Calculate the spatial correlation distribution.
{
  m <- getFactorConcentrationMatrix(dframe, factorName, timestep) ;
  maxManhttanDist <- sum(dim(m) %/% 2);
  distanceSum <- double(length=maxManhttanDist);
  distanceCount <- integer(length=maxManhttanDist);
  # Begin the calculations.
  for (x1 in 1:dim(m)[1])
  {
    for (y1 in 1:dim(m)[2])
    {
      for (x2 in x1:dim(m)[1])
      {
        for (y2 in y1:dim(m)[2])
        {
          d <- getManhattanDistance(x1, y1, x2, y2, dim(m)[1], dim(m)[2]);
          if (d != 0)
          {
            distanceSum[d] <- distanceSum[d] + abs(m[x1, y1] - m[x2, y2]);
            distanceCount[d] <- distanceCount[d] + 1;
          }
        }
      }
    }
  }
  if (length(distanceSum) != length(distanceCount))
  {
    stop("Error in Length of vectors in Spatial Correlation function.");
  }
  # Calculate the mean.
  for (i in 1:length(distanceSum))
  {
    distanceSum[i] <- distanceSum[i] / distanceCount[i];
  }
  return(distanceSum);
}


plotSpatialCorrelation <- function(dframe, factorName, timestep=getMaxTimestep(dframe))
{
  barplot(spatialCorrelation(dframe, factorName, timestep), main=sprintf("Barplot of Spatial Correlation of %s on %d timestep", factorName, as.integer(timestep)), xlab="Manhattan Distance", ylab="Mean of difference in FC");
}


# Some runs.
#lframe <- readTransLattice("onegene.dat");
#m1 <- getFactorConcentrationMatrix(lframe, "f", 1);
#m2 <- getFactorConcentrationMatrix(lframe, "f", 2);
#plotConcentrationSeries(getTimeSlice(lframe, 3), "f", c(0, 1), hitReturn);
# plotConcentrationSeries(lframe, "f", c(0, 1), hitReturn);
# plotConcentrationSeries(lframe, "A", getConcentrationRange(lframe, "A"), oneSecondDelay)

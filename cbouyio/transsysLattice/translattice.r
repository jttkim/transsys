# The .R plotting package for the lattice factor table.

# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

# cbouyio, UEA, 28/09/2006


readTransLattice <- function(filename)
# Reads the file with the headers.
{
  x <- read.table(filename, header = TRUE, comment.char = "#");
  return(x);
}


getTimeSlice <- function(latticeFrame, timesteps)
# Returns a lattice data frame containing only the
# rows of the time steps specified by the timesteps parameter.
{
  # Returns a time slice of the data frame.
  timestepFrame <- subset(latticeFrame, timestep==timesteps);
  return(timestepFrame);
}


getXCoordinates <- function(latticeFrame)
# Returns the X coordinates of a lattice frame.
{
  xCoords <- latticeFrame[["x"]];
  return(unique(xCoords));
}


getYCoordinates <- function(latticeFrame)
# Returns the Y coordinates of a lattice frame.
{
  yCoords <- latticeFrame[["y"]];
  return(unique(yCoords));
}


getTimeSteps <- function(latticeFrame)
# Returns the whole series of timesteps.
{
  return(unique(latticeFrame[["timestep"]]));
}


getMaxTimestep <- function(latticeFrame)
# Returns the time length of the simulation
{
  return(max(getTimeSteps(latticeFrame)));
}


getConcentrationRange <- function(latticeFrame, factorName)
# Returns the factor concentration range.
# (get the whole or slices of the data frame)
{
  cr <- c(min(latticeFrame[[factorName]]), max(latticeFrame[[factorName]]));
  return(cr);
}


getMaximumConcentration <- function(latticeFrame, factorName)
# Returns the maximal factor concentration.
# (get the whole or slices of the data frame)
{
  mc <- max(latticeFrame[[factorName]])
  return(mc)
}


getXSize <- function(latticeFrame)
# Returns the X size of the structure
{
  return(max(getXCoordinates(latticeFrame)));
}


getYSize <- function(latticeFrame)
# Returns the Y size of the structure
{
  return(max(getYCoordinates(latticeFrame)));
}


getTimeseries <- function(latticeFrame, X, Y)
# Returns the timeseries corresponding to the defined instance coordinates.
{
  timeseries <- subset(latticeFrame, x==X & y==Y)
  return(timeseries)
}


getFactorConcentrationMatrix <- function(latticeFrame, factorName, timestep)
# Low level function. Returns a matrix of factor concentrations
# (it keeps the coordinates from the original structure)
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
# Low level function to produce the image of the factor concentration.
{
  # All the job is done by the image function.
  image(xCoordinates, yCoordinates, concentrationMatrix, zlim = concentrationRange, xlab="x Coordinate", ylab="y Coordinate", col = heat.colors(64), main=main, ...);
  legend("right", c(as.character(round(seq(concentrationRange[1], concentrationRange[2], length.out=16), digits=3))), fill=heat.colors(16), title="Gradient", ...);
#  grid(max(xCoordinates), max(yCoordinates));
}


plotConcentrationSeries <- function(latticeFrame, factorName, concentrationRange=c(0, getMaximumConcentration(latticeFrame, factorName)), timeframeEndFunction=function(x){}, ...)
# The main graphics function, plot the image of the factor concentration in
# each timestep.
{
  # Some check about number of timesteps should precede.
  for (i in getTimeSteps(latticeFrame))
  {
    plotConcentrationMatrix(getFactorConcentrationMatrix(latticeFrame, factorName, i), getXCoordinates(latticeFrame), getYCoordinates(latticeFrame), concentrationRange, main=sprintf("Image of %s concentration on timestep %d", factorName, as.integer(i)), ...);
    timeframeEndFunction(i);
  }
}


plotAllInstances <- function(dataFrame, factorName, concentrationRange=c(0, max(dataFrame[[factorName]])), ylim=concentrationRange, ...)
# Plots the timeseries of all the transsys instances on the structure.
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
# timeFrameEndFunction 
{
  print(sprintf("timestep %d", as.integer(timestep)));
  Sys.sleep(1);
}

hitReturn <- function(timestep)
# timeFrameEndFunction
{
  readline(sprintf("timestep %d -- hit return", as.integer(timestep)));
}


getManhattanDistance <- function(x1, y1, x2, y2, X, Y)
# Calculate the Manhattan Distance between two cells on the toroidal lattice.
{
  manhX <- min(abs(x2 - x1), abs(x1 + X - x2), abs(x2 + X - x1))
  manhY <- min(abs(y2 - y1), abs(y1 + Y - y2), abs(y2 + Y - y1))
  manhattanDistance = manhX + manhY
  return(manhattanDistance)
}


spatialCorrelation <- function(mat)
# Returns the spatial correlation distribution of a matrix.
{
  # Calculate the maximum Manhattan Distance.
  maxManhDistance <- sum(dim(mat) %/% 2);
  # The data structure of the correlation distribution.
  spatialDist <- vector("list", maxManhDistance);
  # Begin the calculations.
  for (x1 in 1:dim(mat)[1])
  {
    for (y1 in 1:dim(mat)[2])
    {
      for (x2 in x1:dim(mat)[1])
      {
        for (y2 in y1:dim(mat)[2])
        {
          d <- getManhattanDistance(x1, y1, x2, y2, dim(mat)[1], dim(mat)[2]);
          if (d != 0)
          {
            spatialDist[[d]] <- append(spatialDist[[d]], abs(mat[x1, y1] - mat[x2, y2]));
          }
        }
      }
    }
  }
  return(spatialDist);
}


barplotSP <- function(spatialDist, ...)
# Plots the barplot of the spatial correlation means.
{
  # First calculate the mean (normalize) .
  for (i in 1:length(spatialDist))
  {
    spatialDist[[i]] <- mean(spatialDist[[i]]);
  }
  # Coerce to double.
  spatialDist <- as.double(spatialDist);
  # Plot
  barplot(spatialDist, main="Spatial Correlation Barplot", xlab="Manhattan Distance", ylab="Means of differences in FC", ...);
}


boxplotSP <- function(spatialDist, ...)
# Plots the box plot of the spatial correlation distribution.
{
  boxplot(spatialDist, main="Spatial Correlation Distribution Boxplot", xlab="Manhattan Distance", ylab="Differences in FC", ...);
}



## SVD/PCA analyisis part.

getFrameSlice <- function(latticeFrame, timestep=getMaxTimestep(latticeFrame))
# Returns the factors concentration matrix of all the lattice on the specified
# timestep.
{
  timeSlice <- getTimeSlice(latticeFrame, timestep);
  frameSlice <- timeSlice[,4:length(latticeFrame)];
  return(frameSlice);
}


centeringData <-function(frameSlice)
# Center the values of each expression profile to have mean 0 and SD 1.
{
  for (j in 1:ncol(frameSlice))
  {
    mu <- mean(frameSlice[[j]]);
    sd <- sd(frameSlice[[j]]);
    for (i in 1:nrow(frameSlice))
    {
      frameSlice[i, j] <- (frameSlice[i, j] - mu) / sd ;
    }
  }
  centeredFrameSlice <- frameSlice;
  return(centeredFrameSlice);
}



relativeVariance <- function(frameSlice)
# Draw the relative variance plot of the singular values (i.e. the percentage of
# the variance captured by each of the singular values.)
{
  m <-  centeringData(frameSlice);
  rv <- svd(m)$d**2 / sum(svd(m)$d**2);
  barplot(rv, main='Relative Variance Plot', xlab='Singular Values', ylab='Relative Variance %', ylim=c(0, 1.0));
}


projectComponents <- function(frameSlice)
# Plot the scores of the first two components.
{
  plot(princomp(frameSlice)$scores, main='Scores of the first two Principal Components');
  # plot(predict(princomp(frameSLice, main='Predict of the first two Principal Components')));
}


scatterplotSVD <- function(frame1, frame2, timestep=getMaxTimestep(latticeFrame))
# Draw the scatterplot of the projection of the two first eigengenes.
{
  matrix1 <- centeringMatrix(getMatrix(frame1));
  matrix2 <- centeringMatrix(getMatrix(frame2));
  svd1 <- svd(matrix1);
  svd2 <- svd(matrix2);
  # Calulate the projection of the eigengenes ######### NEEDS FURTHER
  # EXPLANATION....
  xv1 <- svd1$u %*% diag(svd1$d);
  xv2 <- svd2$u %*% diag(svd2$d);

  # Plot the projection of the first two eigengenes.
  plot (xv1[,1], xv1[,2], pch=0)
  points(xv2[,1], xv2[,2], pch=16)
 
}


# Some runs.
#lframe <- readTransLattice("onegene.dat");
#m1 <- getFactorConcentrationMatrix(lframe, "f", 1);
#m2 <- getFactorConcentrationMatrix(lframe, "f", 2);
#plotConcentrationSeries(getTimeSlice(lframe, 3), "f", c(0, 1), hitReturn);
# plotConcentrationSeries(lframe, "f", c(0, 1), hitReturn);
# plotConcentrationSeries(lframe, "A", getConcentrationRange(lframe, "A"), oneSecondDelay)

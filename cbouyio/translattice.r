# The .R plotting package for the lattice factor table.
# cbouyio, UEA, 28/09/2006

readTransLattice <- function(filename)
{
  x <- read.table(filename, header = TRUE);
  return(x);
}


# Returns a lattice data frame containing only the
# rows of the time steps specified by the timesteps parameter.

getTimeSteps <- function(latticeFrame, timesteps)
{
  # Take a subset of the data.frame where the timestep is specified.
  timestepFrame <- subset(latticeFrame, timestep==timesteps);
  return(timestepFrame);
}


getXCoordinates <- function(latticeFrame)
{
  xCoords <- latticeFrame$x;
  return(xCoords);
}


getYCoordinates <- function(latticeFrame)
{
  yCoords <- latticeFrame$y;
  return(xCoords);
}


getConcsntrationRange <- function(latticeFrame, factorName)
{
  cr <- c(min(latticeFrame$factorName), max(latticeFrame$factorName));
  return(cr)
}


# Construct a matrix containing the concentration values
# of the factor specified by factorName.
# Notice that this function cannot properly work if the
# latticeFrame contains multiple time steps.

getFactorConcentrationMatrix <- function(latticeFrame, factorName)
{
  # Populate a matrix of the disered dimensions with NAs.
  m <- matrix(NA, ncol=getYCoordinates(latticeFrame), nrow=getXCoordinates(latticeFrame), byrow=TRUE);
  fc <- latticeFrame$factorName;
  # Populate the matrix with the factor concentrations.
  for (k in fc);
  {
    for (i in getXCoordinates(latticeFrame));
    {
      for (j in getYCoordinates(latticeFrame));
      {
        m[i, j]  <- k;
      }
    }
  }
  # Check for the existance of NA's in the matrix. (validation check)
  for (e in as.double(m));
  {
    if (is.na(e));
    {
      print('Error in factor_concentration matrix population.');
      # The program should exit in this case
    }
  }
  return(m);
}


plotConcentrationMatrix <- function(concentrationMatrix, xCoordinates, yCoordinates, concentrationRange)
{
  # All the job is done by the image function.
  image(concentrationMatrix, axes = FALSE, xlab="x", ylab="y");
  axis(1, 0:(max(xCoordinates))); # Some problems have been observed
  axis(2, 0:(max(yCoordinates))); # in axis drawing
  box();
  gid(max(xCoordinates), max(yCoordinates), col="black", lty="solid");
  # Should also decide the way the gradient will be reproduced.
  # The concentrationRange faunction will be used for that.
}


plotConcentrationSeries <- function(latticeFrame, factorName, concentrationRange) # The concentration range can be obtained from the function above, maybe we can incorporate a timesteps and a delay option as well.
{
  # Some check about number of timesteps should precede.
  delay <- 1
  for i in (0:max(latticeFrame$timestep));
  {
    tframe <- getTimeSteps(latticeFrame, i);
    plotConcentrationMatrix(getFactorConcentrationMatrix(tframe, factorName), getXCoordinates(tframe), getYCoordinates(tframe));
    Sys.sleep(delay);
  }
}


# ./latticeSimulator onegene.tra  onegene.dat
# lframe <- readTransLattice("onegene.dat");
# plotConcentrationMatrix(lframe, "f", c(0, 1));


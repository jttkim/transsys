readTransLattice <- function(filename)
{
  x <- read.table(filename, header = TRUE);
  return(x);
}


# Returns a lattice data frame containing only the
# rows of the time steps specified by the timesteps parameter.

getTimeSteps <- function(latticeFrame, timesteps)
{  
}


getXCoordinates <- function(latticeFrame)
{
  xCoords <- latticeFrame[2];
  return(xCoords);
}


getYCoordinates <- function(latticeFrame)
{
  yCoords <- latticeFrame[3];
  return(xCoords);
}


# Construct a matrix containing the concentration values
# of the factor specified by factorName.
# Notice that this function cannot properly work if the
# latticeFrame contains multiple time steps.

getFactorConcentrationMatrix <- function(latticeFrame, factorName)
{
  m <- matrix(NA, ncol= , nrow= , byrow=TRUE)
}


plotConcentrationMatrix <- function(concentrationMatrix, xCoordinates, yCoordinates, concentrationRange)
{
}


plotConcentrationSeries <- function(latticeFrame, factorName, concentrationRange)
{
}

# ./latticeSimulator onegene.tra  onegene.dat
# lframe <- readTransLattice("onegene.dat");
# plotConcentrationMatrix(lframe, "f", c(0, 1));


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
}


getYCoordinates <- function(latticeFrame)
{
}


# Construct a matrix containing the concentration values
# of the factor specified by factorName.
# Notice that this function cannot properly work if the
# latticeFrame contains multiple time steps.

getFactorConcentrationMatrix <- function(latticeFrame, factorName)
{
}


plotConcentrationMatrix <- function(concentrationMatrix, xCoordinates, yCoordinates, concentrationRange)
{
}


plotConcentrationSeries <- function(latticeFrame, factorName, concentrationRange)
{
}


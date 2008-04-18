# An .R script implementing some usefull functions for the analysis of the
# optimisation experiments output.


plotBestObjective <- function(data)
# Plot the optimisation objective score of r the current best transsys program.
{
plot(data$OptCycle, data$BestObj, type= "l", ylim = c(min(data$BestObj), max(data$BestObj)), main = "Current Best's optimisation objective score", xlab = "Optimisation Cycles", ylab = "Optimisation Objective");
}


plotAltObjective <- function(data)
# Plot the optimisation objective score of r the current alternative transsys
# program.
{
plot(data$OptCycle, data$AltObj, type= "l", ylim = c(min(data$AltObj), max(data$AltObj)), main = "Current Alternative's optimisation objective score", xlab = "Optimisation Cycles", ylab = "Optimisation Objective");
}


plotObjectives <- function(data)
# Function to plot the Objective Score of the current best transsys program
# (red) and the objective of the current alternative (blue). The actual
# optimisation steps (i.e. the cycle where the current alternative is better
# thatn the curren best) are designated by "|" on the plot.
{
  plot(data$OptCycle, data$BestObj, type = "l", col = "red", ylim = c(min(data$BestObj, data$AltObj), max(data$BestObj, data$AltObj)), main = "Optimisation Objective scores of Current Best and Alternative.", xlab = "Optimisation Cycles", ylab = "Optimisation Objective");
  lines(data$OptCycle, data$AltObj, col = "blue");
  dOpt <- subset(data, OptFlag == TRUE);
  points(dOpt$OptCycle, (dOpt$OptFlag - abs(min(data$BestObj, data$AltObj))), pch = "|");
}


extractPostscript <- function(filename, data, h=8, w=10)
# Function to produce a .ps of the plot that is the output of the plotObjective
# function. Takes the filename as an argument.
{
  postscript(filename, horizontal = FALSE, paper = "special", height = h, width = w, onefile = FALSE);
  plotObjectives(data);
  dev.off();
}


bimodalitiesTabular <- function(data)
# Produce a tabular representation of the summary statistics for all the
# calcuated bimodalities.
{
  a <- data.frame(data$AltObj, data$AltLatBM, data$AltCtrlBM, data$BestObj, data$BestLatBM, data$BestCtrlBM);
  return(summary(a))
}


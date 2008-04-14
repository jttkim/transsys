# An .R script implementing some usefull functions for the analysis of the
# optimisation experiments output.


plotObjective <- function(data)
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

extractPostscript <- function(filename, data)
# Function to produce a .ps of the plot that is the output of the plotObjective
# function. Takes the filename as an argument.
{
  postscript(filename, horizontal = FALSE, paper = "special", height = 8, width = 10, onefile = FALSE);
  plotObjective(data);
  dev.off();
}


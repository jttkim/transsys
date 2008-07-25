# An .R script providing some usefull tools for the analyses of the
# optimisation experiments results.


plotOptSteps <- function(data, symbol = 6, ...)
# Low level function to point the actual optimisation steps on a plot.
{
  dOpt <- subset(data, OptFlag == TRUE);
  points(dOpt$OptCycle, dOpt$BestObj, pch = symbol);
}


plotBestObjective <- function(data, title = "Current Best optimisation objective score", ylim = c(min(data$BestObj), max(data$BestObj)), ...)
# Plot the optimisation objective score of the current best transsys program.
{
  plot(data$OptCycle, data$BestObj, type= "l", ylim = ylim, main = title, xlab = "Optimisation Cycles", ylab = "Optimisation Objective");
  plotOptSteps(data, ...);
}


plotAltObjective <- function(data, title = "Current Alternative optimisation objective score", ylim = c(min(data$AltObj), max(data$AltObj)), ...)
# Plot the optimisation objective score of the current alternative transsys
# program.
{
  plot(data$OptCycle, data$AltObj, type= "l", ylim = ylim, main = title, xlab = "Optimisation Cycles", ylab = "Optimisation Objective");
}


plotObjectives <- function(data, title = "Optimisation Objective scores of
Current Best and Alternative.", ylim = c(min(data$BestObj, data$AltObj), max(data$BestObj, data$AltObj)), ...)
# Function to plot the Objective Score of the current best transsys program
# (red) and the objective of the current alternative (blue). The actual
# optimisation steps (i.e. the cycle where the current alternative is better
# thatn the curren best) are designated by "|" on the plot.
{
  plot(data$OptCycle, data$BestObj, type = "l", col = "red", ylim = ylim, main = title, xlab = "Optimisation Cycles", ylab = "Optimisation Objective");
  lines(data$OptCycle, data$AltObj, col = "blue");
  plotOptSteps(data, ...);
}

postscriptObjectives <- function(filename, data, h=8, w=10, ...)
# Function to generate a .ps of the plot of both the objective function curves,
# that is the output of the plotObjectives() function.
# Takes the filename as an argument.
{
  postscript(filename, horizontal = FALSE, paper = "special", height = h, width = w, onefile = FALSE);
  plotObjectives(data, ...);
  dev.off();
}


postscriptBestObjective <- function(filename, data, h=8, w=10, ...)
# Function to generate a .ps of the plot of the transsys Best objective score,
# that is the output of the plotBestObjective() function.
# Takes the filename as an argument.
{
  postscript(filename, horizontal = FALSE, paper = "special", height = h, width = w, onefile = FALSE);
  plotBestObjective(data, ...);
  dev.off();
}


postscriptAltObjective <- function(filename, data, h=8, w=10, ...)
# Function to generate a .ps of the plot of the transsys Alternative objective
# score, that is the output of the plotAltObjective() function.
# Takes the filename as an argument.
{
  postscript(filename, horizontal = FALSE, paper = "special", height = h, width = w, onefile = FALSE);
  plotAltObjective(data, ...);
  dev.off();
}


summaryBimodalities <- function(data, ...)
# Produce a tabular representation of the summary statistics for all the
# calcuated bimodalities.
{
  a <- data.frame(data$AltObj, data$AltLatBM, data$AltCtrlBM, data$BestObj, data$BestLatBM, data$BestCtrlBM);
  Minimum <- c(min(data$AltObj), min(data$AltLatBM), min(data$AltCtrlBM), min(data$BestObj), min(data$BestLatBM), min(data$BestCtrlBM));
  FirstQuartile <- c(quantile(data$AltObj, probs=0.25), quantile(data$AltLatBM, probs=0.25), quantile(data$AltCtrlBM, probs=0.25), quantile(data$BestObj, probs=0.25), quantile(data$BestLatBM, probs=0.25), quantile(data$BestCtrlBM, probs=0.25));
  Median <- c(median(data$AltObj), median(data$AltLatBM), median(data$AltCtrlBM), median(data$BestObj), median(data$BestLatBM), median(data$BestCtrlBM));
  ThirdQuartile <- c(quantile(data$AltObj, probs=0.75), quantile(data$AltLatBM, probs=0.75), quantile(data$AltCtrlBM, probs=0.75), quantile(data$BestObj, probs=0.75), quantile(data$BestLatBM, probs=0.75), quantile(data$BestCtrlBM, probs=0.75));
  IQR <-  c(IQR(data$AltObj), IQR(data$AltLatBM), IQR(data$AltCtrlBM), IQR(data$BestObj), IQR(data$BestLatBM), IQR(data$BestCtrlBM));
  Maximum <- c(max(data$AltObj), max(data$AltLatBM), max(data$AltCtrlBM), max(data$BestObj), max(data$BestLatBM), max(data$BestCtrlBM));
  Mean <- mean(a)
  SD <- sd(a)
  Variance <- sd(a)**2
  d <- data.frame(Minimum, FirstQuartile, ThirdQuartile, Maximum, Median, IQR, Mean, SD, Variance, row.names=c("AltObj:","AltLatBM:","AltCtrlBM","BestObj:","BestLatBM:","BestCtrlBM:"));
  return(d)
}


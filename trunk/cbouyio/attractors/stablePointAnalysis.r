# An R script to reproduce the attractor analysis of the Cherry and Adler
# 2000 paper.

# $Id$

printTranssysProgram <- function(Ka, Kb, Da, Db, n, fileName, ...)
# Print the transsys program.
{
  Apower <- "A";
  Bpower <- "B";
  if (n > 1)
  {
    for (i in c(2:n))
    {
      Apower <- paste(Apower, "* A", sep=" ");
      Bpower <- paste(Bpower, "* B", sep=" ");
    }
  }

  tp <- paste("transsys cherryTests

# A transsys program implementing the Hill's equations described at
# Cherry&Adler2001 paper.

{
  factor A
  {
    decay: ", Da, ";
    diffusibility: 0;
  }

  factor B
  {
    decay: ", Db, ";
    diffusibility: 0;
  }

  gene gene_a
  {
    promoter
    {
      constitutive: (", Ka, " / (1 + (", Bpower, ")));
    }
    product
    {
      default: A;
    }
  }

  gene gene_b
  {
    promoter
    {
      constitutive: (", Kb, " / (1 + (", Apower, ")));
    }
    product
    {
      default: B;
    }
  }
}
", sep="");

  cat(tp, file=fileName);

}


plotFGx <- function(Ka, Kb, Da, Db, n, range, ...)
# Plots the F(G(x)) function fro both factors as well as the diagonal. 
# From Cherry&Adler 2000.
{
  curve(Ka/(Da) / (1 + ((Kb/Db)/(1 + x**n))**n), range[1], range[2], n=10000, lwd=1.5, xlim=range, ylim=range, ylab="F(G(x))", xlab="x", col="blue", ...);
  lines(0:(range[2] + 1), 0:(range[2] + 1), lty="dotted", ...);
}


plotGFy <- function(Ka, Kb, Da, Db, n, range, ...)
# Plots the G(F(y)) function fro both factors as well as the diagonal. 
# transformed from Cherry&Adler 2000.
{
  curve(Kb/(Db) / (1 + ((Ka/Da)/(1 + x**n))**n), range[1], range[2], n=10000,
  lwd=1.5, col="green", ylab="G(F(y))", xlab="Y", xlim=range, ylim=range, ...);
  lines(0:(range[2] + 1), 0:(range[2] + 1), lty="dotted", ...);
}


plotPhaseSpace <- function(Ka, Kb, Da, Db, n, range, ...)
# Plots the phase space curves of the functions.
{
  xy <- seq(0, range[2], by=0.0001);
  plot(xy, (Kb/Db) / (1 + xy**n), xlim=range, ylim=range, xlab="A", ylab="B",
  type="l", col="blue");
  lines((Ka/Da) / (1 + xy**n), xy, type="l", col="green");
}

analyticalSolution <- function(Ka, Kb, Da, Db)
# Define the intersection points analyticaly. (only for Hill exponent 
# equals to two (2))
{
  sA <- polyroot(c(1, -(Ka/Da), 2, -(2*(Ka/Da)), (1 + (Kb/Db)**2), -(Ka/Da)));
  sB <- polyroot(c(1, -(Kb/Db), 2, -(2*(Kb/Db)), (1 + (Ka/Da)**2), -(Kb/Db)));
  s <- rbind(sA, sB)
  rownames(s) <- c("A", "B");
  return(s);
}

plotFGxActivate <- function(Ka, Kb, Da, Db, n, range=c(0, 25), ...)
# Plot the FGx in the Hills activation case.
{
  curve((Ka/Da) / (1 + (((Kb/Db)*x)/(1 + x**n))**n), range[1], range[2], n=1000,
  lwd=1.5, xlim=range, ylim=range, ylab="F(G(x))_activation", xlab="X", col="blue", ...);
  lines(0:(range[2] + 1), 0:(range[2] + 1), lty="dotted", ...);
}


plotGFyActivate <- function(Ka, Kb, Da, Db, n, range=c(0, 10), ...)
# Plot the GFy in the Hills activation case.
{
  curve(((Kb/Db)*(Ka/Da)/(1+ x**n)) / (1 + ((Ka/Da)/(1 + x**n))**n), range[1],
  range[2], n=1000, lwd=1.5, xlim=range, ylim=range, ylab="G(F(y))_activation",
  xlab="Y", col="green", ...);
  lines(0:(range[2] + 1), 0:(range[2] + 1), lty="dotted", ...);
}


drawPrint <-function(Ka, Kb, Da, Db, n, range=c(0,4), fileName, ...)
# Wraps all the functions
{
  plotFG(Ka, Kb, Da, Db, n, range=range, ...);
  sprintf("transsys programs, untill now, can accept only integer Hill
  exponents the parameter n will coerce to integer.");
  n <- as.integer(n);
  printTranssysProgram(Ka, Kb, Da, Db, n, fileName=fileName, ...);
  if (n==2)
  {
    analyticalSolution(Ka, Kb, Da, Db);
  }
}


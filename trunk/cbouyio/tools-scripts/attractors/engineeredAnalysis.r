# R script to conduct printing and analysis of engineered transsys programs.

# $Id$

generateTranssysProgram <- function(Da, Db, Pa, Pb, d, Ha, Hb, La, Lb, Dfa=0, Dfb=0, fileName, ...)
# Export a transsys program to the specified file.
# Da, Db, Decay rates of FactorA, FactorB.
# Pa, Pb, Control point coordinates.
# d, Control circle threshold.
# Ha, Hb, high expression rate of FactorA, FactorB.
# La, Lb, low expression rate of FactorA, FactorB.

{
  tp <- paste("transsys engineered
{
  factor FactorA
  {
    decay: ", Da, ";
    diffusibility: ", Dfa, ";
  }

  factor FactorB
  {
    decay: ", Db, ";
    diffusibility: ", Dfb, ";
  }

  gene geneA
  {
    promoter
    {
      constitutive: ", La, " + (", Ha, "  - ", La, ") * ((FactorA - ", Pa, ") * (FactorA - ", Pa, ") + (FactorB - ", Pb, ") * (FactorB - ", Pb, ") <= ", d, " * ", d, ");
    }
    product
    {
      default: FactorA;
    }
  }

  gene geneB
  {
    promoter
    {
      constitutive: ", Lb, " + (", Hb, "  - ", Lb, ") * ((FactorA - ", Pa, ") * (FactorA - ", Pa, ") + (FactorB - ", Pb, ") * (FactorB - ", Pb, ") <= ", d, " * ", d, ");
    }
    product
    {
      default: FactorB;
    }
  }
}
", sep="");

  cat(tp, file=fileName, ...);
}


drawPhaseSpace <- function(Da, Db, Pa, Pb, d, Ha, Hb, La, Lb, ...)
{
   a <- max(Ha/Da, Hb/Db, Pa+d, Pb+d);
   x <- c(0, a + 0.05*a);
   y <- c(0, a + 0.05*a);
  plot(La/Da, Lb/Db, pch=20, xlim=x, ylim=y, xlab="FactorA", ylab="FactorB", ...);
  abline(h=(y[1]:y[2]), v=(x[1]:x[2]), lty="dotted", col="lightgrey");
  points(Pa, Pb, pch=19, col="green");
  points(Ha/Da, Hb/Db, pch=19, col="red");
  points(La/Da, Lb/Db, pch=19, col="blue");
  circle(Pa, Pb, d);
}


circle <- function(x, y, r, ...)
{
  ang <- seq(0, 2*pi, length = 1000);
  xx <- x + r * cos(ang);
  yy <- y + r * sin(ang);
  polygon(xx, yy, ...);
}

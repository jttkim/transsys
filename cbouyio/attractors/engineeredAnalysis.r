# R script to conduct printing and analysis of engineered transsys programs.

# $Id$

generateTranssysProgram <- function(Da, Di, Pa, Pi, d, Ha, Hi, La, Li, Df1=0, Df2=0, fileName, ...)
# Export a transsys program to the specified file.
# Da, Di, Decay rates of activator, inhibitor (a, i).
# Pa, Pi, control point coordinates.
# d, threshold
# Ha, Hi, high expression rate of a, i.
# La, Li, low expression rate of a, i.

{
  
  tp <- paste("transsys engineered

# A transsys program designed to have a Low and a High equilibrium points.

{
  factor activator
  {
    decay: ", Da, ";
    diffusibility: ", Df1, ";
  }

  factor inhibitor
  {
    decay: ", Di, ";
    diffusibility: ", Df2, ";
  }

  gene activatorgene
  {
    promoter
    {
      constitutive: ", La, " + (", Ha, "  - ", La, " ) * ((activator - ", Pa, ") * (activator - ", Pa, ") + (inhibitor - ", Pi, ") * (inhibitor - ", Pi, ") <= ", d, " * ", d, ");
    }
    product
    {
      default: activator;
    }
  }

  gene inhibitorgene
  {
    promoter
    {
      constitutive: ", Li, " + (", Hi, "  - ", Li, " ) * ((activator - ", Pa, ") * (activator - ", Pa, ") + (inhibitor - ", Pi, ") * (inhibitor - ", Pi, ") <= ", d, " * ", d, ");
    }
    product
    {
      default: inhibitor;
    }
  }
}
", sep="");

  cat(tp, file=fileName, ...);
}


drawPhaseSpace <- function(Da, Di, Pa, Pi, d, Ha, Hi, La, Li, ...)
{
   a <- max(Ha/Da, Hi/Di, Pa+d, Pi+d);
   x <- c(0, a + 0.1*a);
   y <- c(0, a + 0.1*a);
#  symbols(Pa, Pi, circles=d, inches=FALSE, xlim=x, ylim=y, xlab="Factor activator", ylab="Factor inhibitor");
  plot(La/Da, Li/Di, pch=20, xlim=x, ylim=y, xlab="Factor A", ylab="Factor B", ...);
  abline(h=(y[1]:y[2]), v=(x[1]:x[2]), lty="dotted", col="lightgrey");
  points(Pa, Pi, pch=20, col="green");
  points(Ha/Da, Hi/Di, pch=20, col="red");
  points(La/Da, Li/Di, pch=20, col="blue");
  circle(Pa, Pi, d);
}


circle <- function(x, y, r, ...)
{
  ang <- seq(0, 2*pi, length = 1000);
  xx <- x + r * cos(ang);
  yy <- y + r * sin(ang);
  polygon(xx, yy, ...);
}

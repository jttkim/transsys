# This is an .R script to perform several attractor analytical studies.

# Â$Id$

solveQuad <- function(a, b, c)
{
# Check for real roots.
  d <- (b)**2 - 4*(a)*(c);
  if (d < 0)
  {
    stop("The equation has no real roots");
  }
  y1 <- (-(b) + sqrt(d)) / (2*a);
  y2 <- (-(b) - sqrt(d)) / (2*a);
  return(c(y1, y2));
}

plotQuad <- function(a, b, c, ...)
{
  x <- seq(-100, 100, by=0.1);
  y <- vector();
  for (i in x)
  {
    y <- append(y, a*(i**2) + b*i + c)
  }
  plot(x, y, type="l", xlim=c(-5, 5), ylim=c(-5, 5), ...);
#  lines(x**2 + x, x);
}

drawArea <- function(a=-20, b=20)
{
  d <- length(c(a:b));
  plot(rep(0, d), c(a:b), xlim=c(a, b), ylim=c(a, b), type="l", xlab="Factor1", ylab="Factor2");
  abline(h=(a:b), v=(a:b), lty="dotdash", col="lightgrey");
  abline(h=0, v=0);
}


plotFactors <- function(C1, C2, Km1, Km2, Vmax1, Vmax2, D1, D2)
{
  # Draw the plot arrea.
  plot.new();
  a <- -10;
  b <- 10;
  d <- length(c(-10:10));
  plot(rep(0, d), c(a:b), xlim=c(a, b), ylim=c(a, b), type="l", xlab="FactorA",
  ylab="FactorB");
  abline(h=(a:b), v=(a:b), lty="dotdash", col="lightgrey");
  abline(h=0, v=0);

  # Draw the first factor as a function of the second.
  Cf1 <- vector();
  Cf2 <- seq(-10, 10, by=0.1);
  for (x in Cf2)
  {
    Cf1 <- append(Cf1, ((C1*Km2 + C1*x + (Vmax2)*x) / (D1*Km2 + D1*x)))
  }
  lines(Cf1, Cf2);
  abline(h=-Km2, lty="dotdash");

  # Draw the second factor as a function of the first.
  Cf2 <- vector();
  Cf1 <- seq(-10, 10, by=0.1);
  for (y in Cf1)
  {
    Cf2 <- append(Cf2, ((C2*Km1 + C2*y +(Vmax1*y)) / (D2*Km1 + D2*y)))
  }
  lines(Cf1, Cf2, col="blue");
  abline(v=-Km1, lty="dotdash");
}

attractorPoints <- function(C1, C2, Km1, Km2, Vmax1, Vmax2, D1, D2)
{
  # for factor 1.
  a1 <- D1*(D2*Km2 + C2 + Vmax1) ;
  b1 <- D1*D2*Km1*Km2 + C2*D1*Km1 - C1*D2*Km2 - C1*C2 - C1*Vmax1 - C2*Vmax2 - Vmax1*Vmax2 ;
  c1 <- -Km1*(C1*D2*Km2 + C1*C2 + C2*Vmax2) ;
  f1 <- solveQuad(a1, b1, c1);

  # for factor 2.
  a2 <- D2*(D1*Km1 + C1 + Vmax2) ;
  b2 <- D1*D2*Km1*Km2 + C1*D2*Km2 + -C2*D1*Km1 - C1*C2 - C2*Vmax2 -C1*Vmax1 - Vmax1*Vmax2 ;
  c2 <- -(C2*D1*Km1*Km2 + C1*C2*Km2 + C1*Km2*Vmax1) ;
  f2 <- solveQuad(a2, b2, c2);

  points <- rbind(c(f1[1], f2[1]), c(f1[2], f2[2])) ;
  colnames(points) <- c("A", "B") ;
  rownames(points) <- c("p1", "p2") ;
  return(points) ;
}

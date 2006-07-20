"transsysScatter" <- function(transsysFrame, factors = NULL)
{
  if (class(transsysFrame)[1] != "transsysFrame")
  {
    stop("not a transsysFrame");
  }
  if (is.null(factors))
  {
    factors = getAllFactors(transsysFrame);
  }
  nFactors <- length(factors);
  for (i in 1:nFactors)
  {
    if (i < nFactors)
    {
      f1 <- factors[i];
      for (j in (i + 1):nFactors)
      {
        f2 <- factors[j];
        plot(transsysFrame[[f1]], transsysFrame[[f2]]);
        invisible(readline(sprintf("%s vs. %s, hit return", f1, f2)));
      }
    }
  }
}


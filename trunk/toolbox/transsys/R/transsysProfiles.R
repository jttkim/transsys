"transsysProfiles" <- function(transsysFrame, factors = NULL)
{
  if (class(transsysFrame)[1] != "transsysFrame")
  {
    stop("not a transsysFrame");
  }
  if (is.null(factors))
  {
    factors = getAllFactors(transsysFrame);
  }
  for (f in factors)
  {
    plot(transsysFrame[["time"]], transsysFrame[[f]], type = "l", xlab = "time", ylab = f);
    invisible(readline(sprintf("factor %s, hit return", f)));
  }
}

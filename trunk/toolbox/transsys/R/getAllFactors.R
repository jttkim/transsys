"getAllFactors" <- function(transsysFrame)
{
  if (class(transsysFrame)[1] != "transsysFrame")
  {
    stop("not a transsysFrame");
  }
  x <- colnames(transsysFrame);
  return(sort(x[2:length(x)]));
}


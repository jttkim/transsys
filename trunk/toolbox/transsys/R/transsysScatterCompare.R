"transsysScatterCompare" <- function(transsysFrame1, transsysFrame2, factors = NULL)
{
  if (is.null(factors))
  {
    factors <- getAllFactors(transsysFrame1);
  }
  if (!all(getAllFactors(transsysFrame1) %in% factors) || !all(getAllFactors(transsysFrame2) %in% factors))
  {
    stop("factor list not consistent with transsys frames");
  }
  if (!all(transsysFrame1[["time"]] == transsysFrame2[["time"]]))
  {
    stop("incompatible time series");
  }
  for (f in factors)
  {
    x <- transsysFrame1[[f]];
    y <- transsysFrame2[[f]];
    cc <- cor(x, y);
    plot(x, y, main = f, sub = sprintf("corr = %f", cc), xlab = sprintf("%s, set1", f), ylab = sprintf("%s, set 2", f));
    invisible(readline(sprintf("factor %s, hit return", f)));
  }
}


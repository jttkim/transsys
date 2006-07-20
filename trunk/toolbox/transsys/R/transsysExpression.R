"transsysExpression" <- function(transsysProgram, num.timesteps)
{
  d <- runTransexpr(transsysProgram, num.timesteps);
  n <- c("time", colnames(d)[grep("\\.avg$", colnames(d))]);
  d <- d[, n];
  colnames(d) <- sub("\\.avg$", "", n);
  class(d) <- c("transsysFrame", class(d));
  return(d);
}


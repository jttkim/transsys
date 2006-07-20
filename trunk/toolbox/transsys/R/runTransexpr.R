"runTransexpr" <- function(transsysProgram, num.timesteps)
{
  if (length(num.timesteps) != 1)
  {
    stop("num.timesteps must an one-element numeric");
  }
  num.timesteps <- as.integer(num.timesteps);
  cmd = sprintf("transexpr -n %d", num.timesteps);
  x <- xpipe(cmd, transsysProgram);
  tc <- textConnection(x);
  d <- read.table(tc, header = TRUE);
  close(tc);
  return(d);
}


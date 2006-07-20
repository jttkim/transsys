"readTranssysProgram" <- function(f)
{
  tp <- readLines(f);
  class(tp) <- c("transsysProgram", class(tp));
  return(tp);
}


\name{readTranssysProgram}
\alias{readTranssysProgram}
\title{Read a transsys program}
\description{
  Read a \code{transsys} program from a text file.
}
\usage{
readTranssysProgram(f)
}
\arguments{
  \item{f}{A string specifying a file name, or a file object (anything
    that works with \code{read.table}}
}
\details{
  \code{transsys} programs are currently represented by character
  variables containing the program's source code.
}
\value{
  The \code{transsys} program.
}
\references{\url{http://www.cmp.uea.ac.uk/~jtk/transsys/}}
\author{Jan T. Kim, \email{jtk@cmp.uea.ac.uk}}
\examples{
prgtext <- c("transsys example",
            "{",
            "  factor x { decay: 0.1; diffusibility: 0.0; }",
            "}");
f <- textConnection(prgtext);
tp <- readTranssysProgram(f);
close(f);
print(tp);
print(class(tp));
}
\keyword{misc}% at least one, from doc/KEYWORDS

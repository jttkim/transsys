\name{runTransexpr}
\alias{runTransexpr}
\title{Low level interface to the transexpr program}
\description{
  Low level interface to the transexpr program
}
\usage{
runTransexpr(transsysProgram, num.timesteps)
}
\arguments{
  \item{transsysProgram}{the transsys program}
  \item{num.timesteps}{number of time steps to simulate}
}
\details{
  This function invokes the \code{transexpr} program of the
  \code{transsys} system. It creates a \code{transexpr} subprocess,
  feeds the program into the standard input of this process and reads
  the resulting data from the subprocess's standard output. It depends
  on the \code{xpipe} package to do this.

  This function should not be used without specific reason.
  \link{transsysExpression} is a higher level function that should be preferred.
}
\value{
  A data frame containing the expression results. Rows correspond to
  time steps, columns to factors.
}
\references{\url{http://www.cmp.uea.ac.uk/~jtk/transsys/}}
\author{Jan T. Kim, \email{jtk@cmp.uea.ac.uk}}
\examples{
\dontrun{
runTransexpr(tp, 20);
}
}
\keyword{misc}

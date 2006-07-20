\name{transsysExpression}
\alias{transsysExpression}
\title{Generate a synthetic gene expression dataset from a transsys program}
\description{
  Generate a synthetic gene expression dataset from a transsys program
}
\usage{
transsysExpression(transsysProgram, num.timesteps)
}
\arguments{
  \item{transsysProgram}{the \code{transsys} program}
  \item{num.timesteps}{number of timesteps}
}
\details{
  This function produces a subclass of \code{data.frame} which currently
  consists of a column containing the time steps and one column per
  factor in the \code{transsys} program.
}
\value{
  The data frame containing the synthetic gene expression set.
}
\references{\url{http://www.cmp.uea.ac.uk/~jtk/transsys/}}
\author{Jan T. Kim, \email{jtk@cmp.uea.ac.uk}}
\examples{
\dontrun{
transsysExpression(tp, 20);
}
}
\keyword{misc}

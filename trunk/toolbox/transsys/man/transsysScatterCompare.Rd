\name{transsysScatterCompare}
\alias{transsysScatterCompare}
\title{Compare two synthetic gene expression sets by scatter plotting}
\description{
  Compare two synthetic gene expression sets by scatter plotting
}
\usage{
transsysScatterCompare(transsysFrame1, transsysFrame2, factors = NULL)
}
\arguments{
  \item{transsysFrame1}{first transsys frame to be compared}
  \item{transsysFrame2}{second transsys frame to be compared}
  \item{factors}{character vector containing the factors for which
    scatter plots are desired, \code{NULL} means all factors.}
}
\details{
  This function produces a series of scatter plots, each plotting the
  expression profile of a given factor in one data set vs. the profile
  of the same factor in the other set. After each plot, the user is
  requested to hit return in the controlling terminal. This function is
  thus intended for interactive use.

  It is required that (1) the set of factor names is identical in both
  data sets, or at least that the \code{factors} parameter specifies a
  subset common to both datasets and (2) that the time steps for which
  the sets contain data are identical.
}
\value{
  None, the function is used for its side effect.
}
\references{\url{http://www.cmp.uea.ac.uk/~jtk/transsys/}}
\author{Jan T. Kim, \email{jtk@cmp.uea.ac.uk}}
\examples{
\dontrun{
transsysScatterCompare(tp1, tp2);
}
}
\keyword{misc}

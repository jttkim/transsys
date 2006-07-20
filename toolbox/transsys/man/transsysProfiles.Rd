\name{transsysProfiles}
\alias{transsysProfiles}
\title{Draw a series of profile plots of transsys gene expression data}
\description{
  Draw a series of profile plots of transsys gene expression data
}
\usage{
transsysProfiles(transsysFrame, factors = NULL)
}
\arguments{
  \item{transsysFrame}{The \code{transsysFrame} containing the data to
    plot profiles of}
  \item{factors}{Character variable with the factors to plot,
    \code{NULL} means all factors}
}
\details{
  The function generates a series of gene expression profiles, one for
  each factor,
  prompting the user to hit return in the controlling
  terminal after each plot. This is intended for use with an interactive
  graphics display.
}
\value{
  None, function is used for its side effect.
}
\references{\url{http://www.cmp.uea.ac.uk/~jtk/transsys/}}
\author{Jan T. Kim, \email{jtk@cmp.uea.ac.uk}}
\examples{
\dontrun{
transsysProfiles(tp);
}
}
\keyword{misc}

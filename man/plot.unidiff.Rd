\name{plot.unidiff}
\alias{plot.unidiff}
\title{
  Plot Layer Coefficients From a UNIDIFF Model
}
\description{
  Plots the layer coefficient estimates from a UNIDIFF model, together with confidence bars
  based on quasi-standard errors or \dQuote{traditional} standard errors.
}
\usage{
\method{plot}{unidiff}(x,  reference=c("one", "phi"),
     exponentiate = TRUE, se.type = c("quasi.se", "se"),
     conf.int = 0.95, numeric.auto = TRUE, type = "o",
     xlab = names(dimnames(x$data))[3], ylab = "Layer coefficient",
     add = FALSE, ...)
}
\arguments{
  \item{x}{an object resulting from a call to \code{\link{unidiff}}}
  \item{reference}{\dQuote{one} to set the first layer to 1, \dQuote{phi} to set it
      to the value of Goodman's phi association coefficient (see \code{\link{goodman.phi}}).}
  \item{exponentiate}{whether the exponential of the layer coefficients should be plotted
      rather than their value in log scale (see \dQuote{Details} below).}
  \item{se.type}{whether to use quasi-standard errors or \dQuote{traditional} standard errors
      to compute confidence intervals.}
  \item{conf.int}{the confidence level to retain for confidence bars.}
  \item{numeric.auto}{whether layer names should be converted to numeric values when possible
      (see \dQuote{Details} below).}
  \item{type}{what type of plot should be drawn: see \code{\link{plot}}.}
  \item{xlab}{a title for the x axis: see see \code{\link{title}}.}
  \item{ylab}{a title for the y axis: see see \code{\link{title}}.}
  \item{add}{whether to create a new plot using \code{\link{plot}}, or draw over the existing
      plot by calling \code{\link{points}} and \code{\link{segments}} directly.}
  \item{...}{Further arguments passed to \code{\link{plot}}.}
}
\details{
  If \code{exponentiate = TRUE}, coefficients are plotted in the \dQuote{natural}, multiplicative
  scale, rather than in log scale; the first layer coefficient will be equal to 1 rather than 0.
  Confidence intervals will not be symmetric.

  If \code{numeric.auto = TRUE} and layer names (issued from the \code{\link{dimnames}} of the third
  dimension of the original table) can be converted to numeric (i.e. they consist of figures),
  the position of points on the x axis will be determined by the value of the name. This makes most
  sense when layers represent years, especially when they are not regularly spaced. If this behaviour
  is disabled, layers will be placed regularly on the x axis, disregarding their possible interpretation
  as numeric values.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{unidiff}}, \code{\link{summary.unidiff}}
}
\examples{
  # See ?unidiff
}

\keyword{ models }
\keyword{ nonlinear }

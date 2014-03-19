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
\method{plot}{unidiff}(x, what = c("layer.coef", "phi", "maor"),
     se.type = c("quasi.se", "se"),
     conf.int = 0.95, numeric.auto = TRUE, type = "p",
     xlab = names(dimnames(x$data))[3], ylab = NULL,
     add = FALSE, ...)
}
\arguments{
  \item{x}{an object resulting from a call to \code{\link{unidiff}}}
  \item{what}{\dQuote{layer.coefficient} to plot the layer coefficients in the log odds ratio
      scale, with a reference of 1 for the first layer; \dQuote{phi} to plot the intrinsic
      association coefficient (on the log odds ratio scale); \dQuote{maor} to plot the
      mean absolute odds ratio (see \code{\link{maor}}).}
  \item{se.type}{whether to use quasi-standard errors or \dQuote{traditional} standard errors
      to compute confidence intervals.}
  \item{conf.int}{the confidence level to retain for confidence bars.}
  \item{numeric.auto}{whether layer names should be converted to numeric values when possible
      (see \dQuote{Details} below).}
  \item{type}{what type of plot should be drawn: see \code{\link{plot}}. Set to \dQuote{o} or
      \dQuote{b} join points with lines.}
  \item{xlab}{a title for the x axis: see see \code{\link{title}}.}
  \item{ylab}{a title for the y axis: see see \code{\link{title}}; if \code{NULL}, an
      appropriate default is used.}
  \item{add}{whether to create a new plot using \code{\link{plot}}, or draw over the existing
      plot by calling \code{\link{points}} and \code{\link{segments}} directly.}
  \item{...}{Further arguments passed to \code{\link{plot}}.}
}
\details{
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

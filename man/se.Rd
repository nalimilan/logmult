\name{se.assoc}
\alias{se}
\alias{se.assoc}
\alias{se.rc}
\alias{se.hmskew}
\alias{se.yrcskew}
\alias{se.rcL}
\title{Standard Errors for Association Models}
\description{Get standard errors for log-multiplicative association scores and intrinsic
             association coefficients.}
\usage{
  se(x, ...)
  \method{se}{assoc}(x, type = c("se", "quasi.se"), ...)
  \method{se}{rc}(x, type = c("se", "quasi.se"), ...)
  \method{se}{hmskew}(x, type = c("se", "quasi.se"), ...)
  \method{se}{yrcskew}(x, type = c("se", "quasi.se"), ...)
  \method{se}{rcL}(x, type = c("se", "quasi.se"), ...)
}
\arguments{
  \item{ass}{an \code{\link{assoc}} object with a non-null \code{covmat} component.}
  \item{model}{a \code{\link{rc}}, \code{\link{hmskew}}, \code{\link{yrcskew}} or \code{\link{rcL}} object
               fitted with the \code{se} argument different from \code{none}.}
  \item{type}{the type of standard errors to be computed (see \dQuote{Details} below).}
}
\details{
  Quasi-standard errors are computed using \code{\link{qvcalc}}. See the help page for this function for
  details and references about them.
}
\author{
  Milan Bouchet-Valat
}
\seealso{\code{\link{assoc}}, \code{\link{rc}}, \code{\link{hmskew}}, \code{\link{yrcskew}}, \code{\link{rcL}}}
\examples{
  # See ?rc about Wong (2010)
}


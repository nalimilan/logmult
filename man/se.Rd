\name{se.assoc}
\alias{se}
\alias{se.assoc}
\alias{se.rc}
\alias{se.hmskew}
\alias{se.yrcskew}
\alias{se.rcL}
\alias{se.hmskewL}
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
  \item{x}{an \code{\link{assoc}} object with a non-null \code{covmat} component (for \code{se.assoc});.
           or a \code{\link{rc}}, \code{\link{hmskew}}, \code{\link{hmskewL}}, \code{\link{yrcskew}},
           \code{\link{rcL}} or \code{\link{rcL.trans}} object fitted with the \code{se} argument different
           from \dQuote{none} (for other functions).}
  \item{type}{the type of standard errors to be computed (see \dQuote{Details} below).}
  \item{\dots}{currently unused.}
}
\details{
  Currently, only jackknife or bootstrap standard errors are supported, depending on the \code{se} argument
  passed when fitting the model. See \code{\link{summary.assocmod}} for how to ensure the standard errors are correct.

  Quasi-standard errors are computed using \code{\link[qvcalc]{qvcalc}}. See the help page for this function for
  details and references about them.
}
\value{
  An object of the same form as the \code{assoc} component of the model, but with standard errors rather than
  the corresponding coefficients.
}
\author{
  Milan Bouchet-Valat
}
\seealso{\code{\link{assoc}}, \code{\link{rc}}, \code{\link{hmskew}}, \code{\link{hmskewL}}, \code{\link{yrcskew}},
         \code{\link{rcL}}, \code{\link{rcL.trans}}}
\examples{
  # See ?rc about Wong (2010)
}


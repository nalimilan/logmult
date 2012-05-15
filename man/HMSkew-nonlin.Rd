\name{HMSkew}
\alias{HMSkew}
\title{Specify a Skew-Symmetric Association in a gnm Model Formula}
\description{
  A function of class \code{"nonlin"} to specify a van der Heijden & Mooijaart (1995) skew-symmetric
  association in the formula argument to \code{\link{gnm}}.
}
\usage{
  HMSkew(row, col, inst = NULL)
}
\arguments{
  \item{row}{the levels of the row variable}
  \item{col}{the levels of the column variable}
  \item{inst}{a positive integer specifying the instance number of the term}
}
\details{
  This function is used by \code{\link{hmskew}} to fit skew-symmetric models proposed by
  van der Heijden & Mooijaart (1995) and their variants. It can be used directly to fit custom
  variants of the model not supported by \code{hmskew}.

  This function combines its arguments in the following way:
  \deqn{ HMSkew(i, j) = \nu_i \mu_j - \mu_j \nu_i }
  where \eqn{HMSkew(i, j)} is the skew association for the cell at the intersection of row i and column j of
  the table. See reference for mathematical details.
}
\value{
  A list with the required components of a \code{"nonlin"} function:
  \item{predictors}{ the expressions passed to \code{Mult}}
  \item{term}{ a function to create a deparsed mathematical expression
    of the term, given labels for the predictors.}
  \item{call}{ the call to use as a prefix for parameter labels.}
}
\references{
  van der Heijden, P.G.M., and Mooijaart, A. (1995). Some new log bilinear models for the analysis
    of asymmetry in a square contingency table. \emph{Sociol. Methods and Research} 24, 7-29.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{hmskew}}
}
\examples{
  # See ?hmskew.
}

\keyword{ models }
\keyword{ nonlinear }

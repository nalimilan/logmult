\name{RCTrans}
\alias{RCTrans}
\alias{RCTransSymm}
\title{Specify a Row-Column Association With Transitional Layer Effect in a gnm Model Formula}
\description{
  A function of class \code{"nonlin"} to specify a log-multiplicative row-column association models with transitional
  layer effect with one or several dimensions in the formula argument to \code{\link{gnm}}. \code{RCTransSymm}
  allows specifying a variant with symmetric (homogeneous) row and column scores.
}
\usage{
  RCTrans(row, col, layer, inst = NULL)
  RCTransSymm(row, col, layer, inst = NULL)
}
\arguments{
  \item{row}{the levels of the row variable}
  \item{col}{the levels of the column variable}
  \item{layer}{the levels of the layer variable}
  \item{inst}{a positive integer specifying the instance number of the term}
}
\details{
  This function is used by \code{\link{rcL.trans}} to fit an experimental model.

  \code{RCTrans} combines its arguments in the following way:
  \deqn{ RCTrans(i, j, k) = (\mu^S_i + \psi_k \mu^V_i) (\nu^S_j + \psi_k \nu^V_j)}
  where \eqn{RCTrans(i, j, k)} is the skew association for the cell at the intersection of row i, column j
  and layer k of the table.

  \code{RCTransSymm} is similar, but forces \eqn{\mu^S_i} and \eqn{\nu^S_i} (respectively \eqn{\mu^V_i} and
  \eqn{\nu^V_i}) to be equal for identical values of i (diagonal cells).
}
\value{
  A list with the required components of a \code{"nonlin"} function:
  \item{predictors}{ the expressions passed to \code{Mult}}
  \item{term}{ a function to create a deparsed mathematical expression
    of the term, given labels for the predictors.}
  \item{call}{ the call to use as a prefix for parameter labels.}
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{rcL.trans}}
}

\keyword{ models }
\keyword{ nonlinear }

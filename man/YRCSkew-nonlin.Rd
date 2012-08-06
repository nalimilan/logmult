\name{YRCSkew}
\alias{YRCSkew}
\title{Specify a Skew-Symmetric Association in a gnm Model Formula}
\description{
  A function of class \code{"nonlin"} to specify a Yamaguchi (1990) skew-symmetric association in
  the formula argument to \code{\link{gnm}}.
}
\usage{
YRCSkew(mincat, maxcat, skew, inst = NULL)
}
\arguments{
  \item{mincat}{for each cell in the table, the category (row or column) with the lower level}
  \item{maxcat}{for each cell in the table, the category (row or column) with the higher level}
  \item{skew}{1 for cells above the diagonal, -1 for cells below the diagonal, and 0 for diagonal cells}
  \item{inst}{a positive integer specifying the instance number of the term}
}
\details{
  This function is used by \code{\link{yrcskew}} to fit the \dQuote{row-column-effect skew-symmetric
  association (logbilinear) model with full quasi-symmetry (QS+RC_SK)} proposed by Yamaguchi (1990).
  It can be used directly to fit custom variants of the model not supported by \code{yrcskew}.

  This function combines its arguments in the following way:
  \deqn{ YRCSkew(mincat, maxcat, skew) = skew * \nu_{mincat} * (\mu_{maxcat} - \nu_{mincat}) }
  When arguments are set according to what is suggested above, this amounts to the equation:
  \deqn{ YRCSkew_{ij} = \delta_{i<j} \nu_i (\nu_j - \nu_i) - \delta_{i>j} \nu_j (\nu_i - \nu_j)
              = (\delta_{i<j} - \delta_{i>j}) \nu_{min(i,j)} (\nu_{max(i,j)} - \nu_{min(i,j)})}
  where \eqn{YRCSkew_{ij}} is the skew association for the cell at the intersection of row i and column j of
  the table. See reference for mathematical details, and the code of \code{\link{yrcskew}} for real-world usage.
}
\value{
  A list with the required components of a \code{"nonlin"} function:
  \item{predictors}{ the expressions passed to \code{Mult}}
  \item{term}{ a function to create a deparsed mathematical expression
    of the term, given labels for the predictors.}
  \item{call}{ the call to use as a prefix for parameter labels.}
}
\references{
  Yamaguchi, K. (1990). Some Models for the Analysis of Asymmetric Association in Square Contingency Tables with
      Ordered Categories. \emph{Sociol. Methodology} 20, 181-212.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{yrcskew}}
}
\examples{
  # See ?yrcskew.
}

\keyword{ models }
\keyword{ nonlinear }

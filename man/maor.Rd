\name{maor}
\alias{maor}
\title{
  Mean Absolute Odds Ratio or Intrinsic Association Coefficient
}
\description{
  Compute the mean absolute odds ratio of a table, i.e. the (possibly weighted)
  geometric mean of the odds ratios or of their inverse when they are above one,
  which is also equal to the exponential of the intrinsic association coefficient.
  The latter coefficient was first devised by Goodman (1996) as the
  \dQuote{generalized contigency} with a logarithm link, and it is equal to the
  mean of the absolute value of log-linear two-way interaction parameters
  \eqn{\lambda} (in its original version it consists in the square root of
  the sum of squared parameters).

  This index is marginal-free (except of course in its marginal weighted version).
  It is closely related to the association modeled by UNIDIFF and RC models.
}
\usage{
maor(tab, log = FALSE,
     weighting = c("marginal", "uniform", "none"),
     norm = 2, row.weights = NULL, col.weights = NULL)
}
\arguments{
  \item{tab}{a two-way table without zero cells.}
  \item{log}{if \dQuote{TRUE}, return the intrinsic association coefficient
      rather than the mean absolute odds ratio.}
  \item{weighting}{what weights should be used when normalizing the scores.}
  \item{norm}{the norm to use to compute the mean of \eqn{\lambda} parameters,
      typically 1 for the mean of absolute values, or 2 for the square root
      of the sum of squared parameters (as in the original version).}
  \item{row.weights}{optional custom weights to be used for rows, e.g. to compute
      the phi coefficient for several tables using their overall marginal
      distribution. If specified, \code{weighting} is ignored.}
  \item{col.weights}{see \code{col.weights}.}
}
\details{
  See Goodman (1996), Equation 52 for the (marginal or other) weighted version;
  the unweighted version can be computed with unit weights. The coefficient is
  called \eqn{\tilde \lambda^2} in the original article, but to avoid the confusion
  with Goodman and Kruskal's lambda coefficient, it is here referred to as the
  intrinsic association coefficient, denoted as \eqn{\phi} in row-column association
  models.

  This index cannot be computed in the presence of zero cells since it is based
  on the logarithm of proportions. In these cases, 0.5 is added to all cells in
  the table (Agresti 2002, p. 249; Berkson 1955), and a warning is printed.
  Make sure this correction does not affect too much the results (especially
  with small samples) by manually adding different values before calling this function.
}
\value{
  The numeric value of the mean absolute odds ratio, or of the intrinsic association
  coefficient (if \code{log = TRUE}).
}
\references{
  Agresti, A. 2002. Categorical Data Analysis. New York: Wiley.

  Goodman, L. A. (1996). A Single General Method for the Analysis of Cross-Classified Data: Reconciliation
      and Synthesis of Some Methods of Pearson, Yule, and Fisher, and Also Some Methods of Correspondence
      Analysis and Association Analysis. \emph{J. of the Am. Stat. Ass.} 91(433):408-428.

  Berkson, J. (1955). Maximum Likelihood and Minimum chi2 Estimates of the Logistic Function.
      \emph{J. of the Am. Stat. Ass.} 50(269):130-162.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{unidiff}}, \code{\link{rc}}
}

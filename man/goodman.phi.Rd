\name{goodman.phi}
\alias{goodman.phi}
\title{
  Compute Goodman's (1996) Phi Association Coefficient
}
\description{
  Compute the value of the association coefficient first devised by Goodman (1996),
  called the \dQuote{generalized contigency} for the logarithm link.
  This coefficient is derived from log-linear two-way interaction parameters (\eqn{\lambda}),
  and as such marginal-free (except in its marginal weighted version). It is closely
  related to the association modeled by UNIDIFF and RC models.
}
\usage{
goodman.phi(tab,
            weighting = c("marginal", "uniform", "none"),
            row.weights = NULL, col.weights = NULL)
}
\arguments{
  \item{tab}{a two-way table without zero cells.}
  \item{weighting}{what weights should be used when normalizing the scores.}
  \item{row.weights}{optional custom weights to be used for rows, e.g. to compute
      the phi coefficient for several tables using their overall marginal
      distribution. If specified, \code{weighting} is ignored.}
  \item{col.weights}{see \code{col.weights}.}
}
\details{
  See Goodman (1996), Equation 52 for the (marginal or other) weighted version;
  the unweighted version can be computed with unit weights. The coefficient is
  called \eqn{\tilde \lambda^2} in the original article, but to avoid the confusion
  with Goodman and Kruskal's lambda coefficient, it is here referred to as \eqn{\phi}
  due to the homology with both Pearson's \dQuote{mean square contingency} coefficient
  and the \dQuote{intrinsic association coefficient} of RC models.

  This index cannot be computed in the presence of zero cells since it is based
  on the logarithm of proportions. Consider adding a small quantity like 0.5 to
  these cells (Agresti 2002, p. 249; Berkson 1955).
}
\value{
  The numeric value of the coefficient (i.e. the square root of the \eqn{\tilde \lambda^2}
  in Goodman (1996) Equation 52).
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

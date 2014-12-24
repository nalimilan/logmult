\name{maor}
\alias{maor}
\title{
  Mean Absolute Odds Ratio or Intrinsic Association Coefficient
}
\description{
  Compute the mean absolute odds ratio of a table, i.e. the (possibly weighted)
  geometric mean of the odds ratios or of their inverse when they are above one,
  which is also closely related to the the intrinsic association coefficient.
  The latter coefficient was first devised by Goodman (1996) as the
  \dQuote{generalized contigency} when a logarithm link is used, and it is equal
  to the mean of the absolute value of log-linear two-way interaction parameters
  \eqn{\lambda_{ij}} (in its original version it consists in the square root of
  the sum of squared parameters).
}
\usage{
maor(tab, phi = FALSE, cell = FALSE,
     weighting = c("marginal", "uniform", "none"),
     norm = 2, row.weights = NULL, col.weights = NULL)
}
\arguments{
  \item{tab}{a two- or three-way table without zero cells; for three-way tables,
      average marginal weighting is used when \dQuote{weighting = "marginal"}, and
      the MAOR is computed for each layer (third dimension).}
  \item{phi}{if \dQuote{TRUE}, return the intrinsic association coefficient
      rather than the Mean absolute odds ratio.}
  \item{cell}{if \dQuote{TRUE}, return the per-cell contributions
       (affected by the value of \code{phi}, see \dQuote{Details} below).}
  \item{weighting}{what weights should be used when normalizing the scores.}
  \item{norm}{the norm to use to compute the mean of \eqn{\lambda_{ij}} parameters,
      1 for the mean of absolute values, or 2 for the square root
      of the sum of squared parameters (as in the original version).}
  \item{row.weights}{optional custom weights to be used for rows, e.g. to compute
      the phi coefficient for several tables using their overall marginal
      distribution. If specified, \code{weighting} is ignored.}
  \item{col.weights}{see \code{row.weights}.}
}
\details{
  See Goodman (1996), Equation 52 for the (marginal or other) weighted version of
  the intrinsic association coefficient (\eqn{\tilde \phi}); the unweighted version can be
  computed with unit weights. The coefficient is called \eqn{\tilde \lambda^2} in the original
  article, but to avoid the confusion with Goodman and Kruskal's lambda coefficient,
  it is here denoted as \eqn{\phi}, as usual in row-column association models.
  The uniform-weighed version is defined as:
  \deqn{\phi = \left( \sum_{i = 1}^I \sum_{j = 1}^J \frac{1}{IJ} |\lambda_{ij}|^n \right)^{1/n}}
  The (marginal or other) weighted version is defined as:
  \deqn{\tilde \phi = \left( \sum_{i = 1}^I \sum_{j = 1}^J P_{i+} P_{+j} |\tilde \lambda_{ij}|^n \right)^{1/n}}
  with \eqn{n} the chosen norm (typically 1 or 2), \eqn{\sum_{i = 1}^I \lambda_{ij} = \sum_{j = 1}^J \lambda_{ij} = 0},
  \eqn{\sum_{i = 1}^I P_{i+} \tilde \lambda_{ij} = \sum_{j = 1}^J P_{+j} \tilde \lambda_{ij} = 0}.

  The uniform-weighted version of the mean absolute odds ratio (MAOR) is defined as:
  \deqn{MAOR = \exp \left[ \left( \frac{2 \phi^n}{IJ (I-1) (J-1)} \right)^{1/n} \right]}
  The (marginal or other) weighted version is defined as:
  \deqn{MAOR = \exp \left[ \left( \frac{2 \tilde \phi^n}{\sum_{i = 1}^I \sum_{j = 1}^J P_{i+}
      (1 -  P_{i+}) P_{+j} (1 - P_{+j})} \right) ^{1/n} \right]}

  Per-cell contributions \eqn{c_{ij}} are defined so that \eqn{\tilde \phi = \left( {\sum_{i = 1}^I \sum_{j = 1}^J c_{ij} } \right)^{1/n}} \cr
  and \eqn{MAOR = exp \left[ \left( {\sum_{i = 1}^I \sum_{j = 1}^J c_{ij}} \right)^{1/n} \right]}.

  This index cannot be computed in the presence of zero cells since it is based
  on the logarithm of proportions. In these cases, 0.5 is added to these cells
  (Agresti 2002, p. 249; Berkson 1955), and a warning is printed.
  Make sure this correction does not affect too much the results (especially
  with small samples) by manually adding different values before calling this function.
}
\value{
  The numeric value of the mean absolute odds ratio, or of the intrinsic association
  coefficient (if \code{phi = TRUE}), or the corresponding per-cell contributions
  (if \code{cell = TRUE}).
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

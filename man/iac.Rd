\name{iac}
\alias{iac}
\title{
  Intrinsic Association Coefficient
}
\description{
  Compute the intrisic association coefficient of a table.
  This coefficient was first devised by Goodman (1996) as the
  \dQuote{generalized contingency} when a logarithm link is used, and it is equal
  to the standard deviation of the log-linear two-way interaction parameters
  \eqn{\lambda_{ij}}.
  To obtain the Altham index, multiply the result by
  \code{sqrt(nrow(tab) * ncol(tab)) * 2} (see \dQuote{Examples} below).
}
\usage{
iac(tab, cell = FALSE,
    weighting = c("marginal", "uniform", "none"),
    component = c("total", "symmetric", "antisymmetric"),
    shrink = FALSE,
    row.weights = NULL, col.weights = NULL)
}
\arguments{
  \item{tab}{a two- or three-way table without zero cells; for three-way tables,
      average marginal weighting is used when \dQuote{weighting = "marginal"}, and
      the MAOR is computed for each layer (third dimension).}
  \item{cell}{if \dQuote{TRUE}, return the per-cell contributions
       (affected by the value of \code{phi}, see \dQuote{Details} below).}
  \item{weighting}{what weights should be used when normalizing the scores.}
  \item{component}{whether to compute the total association, or from symmetric
      or antisymmetric interaction coefficients only.}
  \item{shrink}{whether to use the empirical Bayes shrinkage estimator proposed by Zhou (2015)
      rather than the direct estimator.}
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
  The uniform-weighted version is defined as:
  \deqn{\phi = \sqrt{ \frac{1}{IJ} \sum_{i = 1}^I \sum_{j = 1}^J \lambda_{ij}^2 }}
  The (marginal or other) weighted version is defined as:
  \deqn{\tilde \phi = \sqrt{ \sum_{i = 1}^I \sum_{j = 1}^J \tilde \lambda_{ij}^2 P_{i+} P_{+j} }}
  with \eqn{\sum_{i = 1}^I \lambda_{ij} = \sum_{j = 1}^J \lambda_{ij} = 0} and
  \eqn{\sum_{i = 1}^I P_{i+} \tilde \lambda_{ij} = \sum_{j = 1}^J P_{+j} \tilde \lambda_{ij} = 0}.

  Per-cell contributions \eqn{c_{ij}} are defined so that:
  \eqn{\tilde \phi = \sqrt{ \sum_{i = 1}^I \sum_{j = 1}^J c_{ij} }}. For the unweighted case,
  \eqn{c_{ij} = \lambda_{ij}^2 / IJ}; for the weighted case,
  \eqn{\tilde c_{ij} = \tilde \lambda_{ij}^2 P_{i+} P_{+j}}.

  This index cannot be computed in the presence of zero cells since it is based
  on the logarithm of proportions. In these cases, 0.5 is added to all cells of the table
  (Agresti 2002, sec. 9.8.7, p. 397; Berkson 1955), and a warning is printed.
  Make sure this correction does not affect too much the results (especially
  with small samples) by manually adding different values before calling this function.
}
\value{
  The numeric value of the intrinsic association coefficient (if \code{cell = FALSE}),
  or the corresponding per-cell contributions (if \code{cell = TRUE}).
}
\references{
  Agresti, A. 2002. Categorical Data Analysis. New York: Wiley.

  Altham, P. M. E., Ferrie J. P., 2007. Comparing Contingency Tables Tools for Analyzing Data
      from Two Groups Cross-Classified by Two Characteristics. \emph{Historical Methods} 40(1):3-16.

  Bouchet-Valat, M. (in press). General Marginal-free Association Indices for Contingency Tables:
  From the Altham Index to the Intrinsic Association Coefficient. \emph{Sociological Methods & Research}.
  Online first.

  Berkson, J. (1955). Maximum Likelihood and Minimum chi2 Estimates of the Logistic Function.
      \emph{J. of the Am. Stat. Ass.} 50(269):130-162.

  Goodman, L. A. (1996). A Single General Method for the Analysis of Cross-Classified Data: Reconciliation
      and Synthesis of Some Methods of Pearson, Yule, and Fisher, and Also Some Methods of Correspondence
      Analysis and Association Analysis. \emph{J. of the Am. Stat. Ass.} 91(433):408-428.

  Zhou, X. (2015). Shrinkage Estimation of Log-Odds Ratios for Comparing Mobility Tables.
      \emph{Sociological Methodology} 45(1):33-63.
}
\author{
  Milan Bouchet-Valat
}
\examples{
  # Altham index (Altham and Ferrie, 2007, Table 1, p. 3 and commentary p. 8)
  tab1 <- matrix(c(260, 195, 158, 70,
                   715, 3245, 874, 664,
                   424, 454, 751, 246,
                   142, 247, 327, 228), 4, 4)
  iac(tab1, weighting="n") * sqrt(nrow(tab1) * ncol(tab1)) * 2

  # Zhou (2015)
  data(hg16)
  # Add 0.5 due to the presence of zero cells
  hg16 <- hg16 + 0.5
  # Figure 3, p. 343: left column then right column
  # (reported values are actually twice the Altham index)
  iac(hg16, weighting="n") * sqrt(nrow(hg16) * ncol(hg16)) * 2 * 2
  iac(hg16, weighting="n", shrink=TRUE) * sqrt(nrow(hg16) * ncol(hg16)) * 2 * 2
  # Table 4, p. 347: values are not exactly the same
  u <- unidiff(hg16)
  # First row
  cor(u$unidiff$layer$qvframe$estimate, iac(hg16, weighting="n"))
  cor(u$unidiff$layer$qvframe$estimate, iac(hg16, weighting="n"), method="spearman")
  # Second row
  cor(u$unidiff$layer$qvframe$estimate, iac(hg16, shrink=TRUE, weighting="n"))
  cor(u$unidiff$layer$qvframe$estimate, iac(hg16, shrink=TRUE, weighting="n"), method="spearman")
}
\seealso{
  \code{\link{unidiff}}, \code{\link{rc}}, \code{\link{maor}}
}

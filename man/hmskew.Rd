\name{hmskew}
\alias{hmskew}
\alias{print.hmskew}
\title{Fitting van der Heijen & Mooijaart Skew-Symmetric Association Model}
\description{
  \code{hmskew} fits a skew-symmetric association model proposed in
  van der Heijden & Mooijaart (1995) to describe asymmetry of square tables.
  This model can be combined with symmetric association models like a quasi-symmetry
  (the default) or symmetric (homogeneous) RC(M) models.
}
\usage{
hmskew(tab, nd.symm = NA, diagonal = FALSE,
       weighting = c("marginal", "uniform", "none"),
       se = c("none", "jackknife"), family = poisson,
       start = NA, etastart = NULL, tolerance = 1e-6,
       iterMax = 15000, trace = TRUE, ...)
}
\arguments{
  \item{tab}{ a two-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above two will be collapsed as appropriate.}
  \item{nd.symm}{ the number of dimensions to include in the \emph{symmetric} RC(M) association. Cannot exceed
    \code{2 * min(nrow(tab) - 1, ncol(tab) - 1)} (the model is already saturated at that point).
    If \code{NA}, a quasi-symmetric association is used instead of a RC(M) model.}
  \item{diagonal}{ should the model include parameters specific to each diagonal cell? This amounts to
    taking quasi-independence, rather than independence, as the baseline model.}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{se}{ whether to compute standard errors for scores. \dQuote{jackknife} is the only currently
     supported method, and is disabled by default because it increases computation time.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{start}{either \code{NA} (the default) to use reasonable starting values, \code{NULL} to use
     random starting values, or a vector of starting values for the parameters in the model;
     if a starting value is \code{NA}, the default starting value will be used.}
  \item{etastart}{ starting values for the linear predictor; set to \code{NULL} to use either default
     starting values (if \code{start = NA}), or random starting values (in all other cases).}
  \item{tolerance}{ a positive numeric value specifying the tolerance level for
     convergence; higher values will speed up the fitting process, but beware of numerical
     instability of estimated scores!}
  \item{iterMax}{ a positive integer specifying the maximum number of main iterations to perform;
     consider raising this value if your model does not converge.}
  \item{trace}{ a logical value indicating whether the deviance
     should be printed after each iteration.}
  \item{\dots}{ more arguments to be passed to \code{\link{gnm}}}
}
\details{
  The original model presented by van der Heijden & Mooijaart (1995), called \dQuote{quasi-symmetry plus
  skew-symmetry}, combines a skew-symmetric association with a quasi-symmetry baseline; it is the variant
  fitted by default by this function. If \code{nd.symm} is set to a positive integer value, though, variants
  using a RC(M) model to describe the \emph{symmetric association} are used, with our without
  diagonal-specific parameters (depending on the value of the \code{diagonal} argument).

  These models follow the equation:
  \deqn{ log F_{ij} = q_{ij} + \nu_i \mu_j - \mu_j \nu_i }
  where \eqn{F_ij} is the frequency observed in the cell at the intersection of row i and column j of
  \code{tab}, and \eqn{q_{ij}} a quasi-symmetric or a RC(M) association. See reference for detailed information
  about the degrees of freedom and the identification constraints applied to the scores.

  Another model presented in the paper, the \dQuote{symmetry plus skew-symmetry model} is not currently supported
  out of the box, but should be relatively straightforward to implement using the underlying \code{\link{assoc.hmskew}}
  function combined with a symmetric association model.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife standard errors, and plotting.
  The default starting values correspond to the parameters computed for the model without association
  parameters (\dQuote{base model}); association parameters start with random starting values. In some complex
  cases, using \code{start = NA} to get completely random starting values can be more efficient, but it is also
  less stable.
}
\value{
  A \code{hmskew} object, which is a subclass of an \code{rc.symm} object (see \code{\link{rc}}) unless
    \code{nd} is \code{NA}. In addition to this class, it contains a \code{assoc.hmskew} component holding
    information about the \emph{asymmetric} association:
  \item{phi }{The intrisic association parameters, one per dimension.}
  \item{row }{Row scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{row.weights }{The row weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{col.weights }{The column weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{covmat }{The variance-covariance matrix of the parameter's estimators, if \code{se}
    was set. It contains first phi coefficients, then normalized row and column scores, then
    sqrt(phi) * row/column scores.}
  \item{covtype }{The method used to compute the variance-covariance matrix. Currently, only
    \dQuote{jackknife} is supported.}
}
\references{
  van der Heijden, P.G.M., and Mooijaart, A. (1995). Some new log bilinear models for the analysis
    of asymmetry in a square contingency table. \emph{Sociol. Methods and Research} 24, 7-29.
}
\author{
Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.hmskew}}, \code{\link{gnm}}
}
\examples{
  ## van der Heijden & Mooijaart (1995), Table 2c, p. 23
  data(ocg1973)
  # 5:1 is here to take "Farmers" as reference category (angle 0)
  model <- hmskew(ocg1973[5:1, 5:1], weighting="uniform")
  ass <- model$assoc.hmskew
  # First column of the table
  round(ass$row[,,1] * sqrt(ass$phi[1,1]), d=2)[5:1,]
  \dontshow{
      # First score for Farmers is slightly different from the original article
      stopifnot(isTRUE(all.equal(round(ass$row[,,1] * sqrt(ass$phi[1,1]), d=2)[5:1,],
                                       matrix(c(-0.08, -0.2, -0.23, -0.11, 0.61,
                                                 0.34,  0.3, -0.13, -0.51, 0), 5, 2),
                                 check.attributes=FALSE)))
  }
  # Right part of the table
  round(ass$phi[1] * (ass$row[,2,1] \%o\% ass$row[,1,1] - ass$row[,1,1] \%o\% ass$row[,2,1]), d=3)[5:1, 5:1]
  # Plot
  plot(model, coords="cartesian")
}

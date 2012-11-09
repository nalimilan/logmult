\name{hmskew}
\alias{hmskew}
\alias{print.hmskew}
\title{Fitting van der Heijen & Mooijaart Skew-Symmetric Association Model}
\description{
  Fits a skew-symmetric association model proposed in van der Heijden & Mooijaart (1995)
  to describe asymmetry of square tables. Skew-symmetric association can be combined with
  quasi-symmetry (the default), quasi-independence, or symmetric (homogeneous) RC(M) associations.
}
\usage{
hmskew(tab, nd.symm = NA, diagonal = FALSE,
       weighting = c("marginal", "uniform", "none"),
       se = c("none", "jackknife", "bootstrap"),
       nreplicates = 100, ncpus = getOption("boot.ncpus"),
       family = poisson, weights = NULL,
       start = NA, etastart = NULL, tolerance = 1e-8,
       iterMax = 5000, trace = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{tab}{ a square two-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above two will be collapsed.}
  \item{nd.symm}{ the number of dimensions to include in the \emph{symmetric} RC(M) association. Cannot exceed
    \code{2 * min(nrow(tab) - 1, ncol(tab) - 1)} (the model is already saturated at that point).
    If \code{NA} (the default), a full quasi-symmetric association is used instead of a RC(M) model; if \code{0},
    quasi-independence is used.}
  \item{diagonal}{ should the model include parameters specific to each diagonal cell? This amounts to
    taking quasi-independence, rather than independence, as the baseline model.}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{se}{ which method to use to compute standard errors for parameters.}
  \item{nreplicates}{ the number of bootstrap replicates, if enabled.}
  \item{ncpus}{ the number of processes to use for jackknife or bootstrap parallel computing. Defaults to
    the number of cores (see \code{\link[parallel]{detectCores}}), with a maximum of 5, but falls back to 1
    (no parallelization) if packages \code{parallel} or \code{snow} are not available.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{weights}{ an optional vector of weights to be used in the fitting process.}
  \item{start}{either \code{NA} (the default) to use optimal starting values, \code{NULL} to use
     random starting values, or a vector of starting values for the parameters in the model.}
  \item{etastart}{ starting values for the linear predictor; set to \code{NULL} to use either default
     starting values (if \code{start = NA}), or random starting values (in all other cases).}
  \item{tolerance}{ a positive numeric value specifying the tolerance level for
     convergence; higher values will speed up the fitting process, but beware of numerical
     instability of estimated scores!}
  \item{iterMax}{ a positive integer specifying the maximum number of main iterations to perform;
     consider raising this value if your model does not converge.}
  \item{trace}{ a logical value indicating whether the deviance
     should be printed after each iteration.}
  \item{verbose}{ a logical value indicating whether progress indicators should be printed,
     including a diagnostic error message if the algorithm restarts.}
  \item{\dots}{ more arguments to be passed to \code{\link{gnm}}}
}
\details{
  The original model presented by van der Heijden & Mooijaart (1995), called \dQuote{quasi-symmetry plus
  skew-symmetry}, combines a skew-symmetric association with a quasi-symmetry baseline; it is the variant
  fitted by default by this function. If \code{nd.symm} is set to a positive integer value, though, variants
  using a RC(M) model to describe the \emph{symmetric association} are used, with our without
  diagonal-specific parameters (depending on the value of the \code{diagonal} argument).

  These models follow the equation:
  \deqn{ log F_{ij} = q_{ij} + \phi (\nu_i \mu_j - \mu_j \nu_i) }
  where \eqn{F_{ij}} is the expected frequency for the cell at the intersection of row i and column j of
  \code{tab}, and \eqn{q_{ij}} a quasi-symmetric specification, with either full interaction parameters, or
  a RC(M) association. See reference for detailed information about the degrees of freedom and the identification
  constraints applied to the scores.

  Another model presented in the paper, the \dQuote{symmetry plus skew-symmetry model} is not currently supported
  out of the box, but should be relatively straightforward to implement using the underlying \code{\link{assoc.hmskew}}
  function combined with a symmetric association model.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife or bootstrap standard errors, and plotting.
  The default starting values for skew association parameters are computed using an eigen value decomposition from the
  results of the model without skew association component (\dQuote{base model}); if \code{nd.symm} is not \code{NA} and
  strictly positive, random starting values are used. In some complex cases, using \code{start = NULL} to start with
  random values can be more efficient, but it is also less stable and can converge to non-optimal solutions.
}
\value{
  A \code{hmskew} object, which is a subclass of an \code{rc.symm} object (see \code{\link{rc}}) if
    \code{nd.symm} is strictly positive. In addition to this class, it contains a \code{assoc.hmskew} component
    holding information about the \emph{skew-symmetric} association:
  \item{phi }{The intrisic association parameters, one per dimension.}
  \item{row }{Row scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{row.weights }{The row weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{col.weights }{The column weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{covmat }{The variance-covariance matrix for phi coefficients and normalized row and column
    scores. Only present if \code{se} was not \dQuote{none}.}
  \item{adj.covmats }{An array stacking on its third dimension one variance-covariance matrix for
    the adjusted scores of each layer in the model (used for plotting). Only present if \code{se}
    was not \dQuote{none}.}
  \item{covtype }{The method used to compute the variance-covariance matrix (corresponding to the
    \code{se} argument.}
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
  model
  ass <- model$assoc.hmskew

  # First column of the table
  round(ass$row[,,1] * sqrt(ass$phi[1,1]), d=2)[5:1,]

  # Right part of the table
  round(ass$phi[1] * (ass$row[,2,1] \%o\% ass$row[,1,1] - ass$row[,1,1] \%o\% ass$row[,2,1]), d=3)[5:1, 5:1]

  # Plot
  plot(model, coords="cartesian")
}

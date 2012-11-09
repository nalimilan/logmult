\name{rcL}
\alias{rcL}
\alias{print.rcL}
\alias{print.rcL.symm}
\title{Fitting Row-Column Association Models With Layer Effect}
\description{
  Fit log-multiplicative row-column association models with layer effect, also called RC(M)-L models,
  with one or several dimensions. Supported variants include homogeneous or heterogeneous scores over
  the layer variable, and (for square tables) symmetric (homogeneous) row and column scores, possibly
  combined with separate diagonal parameters.
}
\usage{
rcL(tab, nd = 1,
    layer.effect = c("homogeneous.scores", "heterogeneous", "none"),
    symmetric = FALSE,
    diagonal = c("none", "heterogeneous", "homogeneous"),
    weighting = c("marginal", "uniform", "none"),
    se = c("none", "jackknife", "bootstrap"),
    nreplicates = 100, ncpus = getOption("boot.ncpus"),
    family = poisson, weights = NULL,
    start = NA, etastart = NULL, tolerance = 1e-8,
    iterMax = 5000, trace = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{tab}{ a three-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above three will be collapsed.}
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{layer.effect}{ determines the form of the interaction between row-column association and layers.
     See \dQuote{Details} below.}
  \item{symmetric}{ should row and column scores be constrained to be equal? Valid only for square tables.}
  \item{diagonal}{ what type of diagonal-specific parameters to include in the model, if any. This amounts to
    taking quasi-conditional independence, rather than conditional independence, as the baseline model.
    Valid only for square tables.}
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
  This function fits log-multiplicative row-column association models with layer effect, usually called (after
  Wong) RC(M)-L models, typically following the equation:
  \deqn{ log F_{ijk} = \lambda + \lambda^I_i + \lambda^J_j + \lambda^K_k
                       + \lambda^{IK}_{ik} + \lambda^{JK}_{jk}
                       + \sum_{m=1}^M { \phi_{km} \mu_{ikm} \nu_{jkm} } }
  where \eqn{F_{ijk}} is the expected frequency for the cell at the intersection of row i, column j and layer k of
  \code{tab}, and M the number of dimensions. If \code{layer.effect} is set to \sQuote{heterogeneous}, different scores
  will be computed for each level, which is equivalent to fitting separate RC(M) models on the k two-way tables.
  If it is set to \sQuote{homogeneous.scores}, then \eqn{\mu_{ikm} = \mu_{im}} and \eqn{\nu_{ikm} = \nu_{im}} for all
  layers k: only the \eqn{\phi_{km}} are allowed to vary across layers. If it is set to \sQuote{none}, then in addition
  to the previous conditions all \eqn{\phi_{km}} are forced to be equal for all layers k, which amounts to a stability
  of the association across layers. See references for detailed information about the variants of the model, the degrees
  of freedom and the identification constraints applied to the scores.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife or bootstrap standard errors, and plotting.
  The default starting values for association parameters are computed using a singular/eigen value decomposition
  from the results of the model without association component (\dQuote{base model}). In some complex cases, using
  \code{start = NULL} to start with random values can be more efficient, but it is also less stable and can converge
  to non-optimal solutions.
}
\value{
  A \code{rcL} object, with all the components of a \code{\link{gnm}} object, plus an
    \code{assoc} component holding the most relevant association information:
  \item{phi }{The intrisic association parameters, one per dimension and per layer.}
  \item{row }{Row scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{weighting }{The name of the weighting method used, reflected by \code{row.weights}
    and \code{col.weights}.}
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
  Wong, R.S-K. (2010). Association models. SAGE: Quantitative Applications in the Social Sciences.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.rcL}}, \code{\link{gnm}}
}
\examples{
  ## Becker & Clogg (1989), Table 5 (p. 145)
  # See also ?rc for more details
  \dontrun{
  data(color)

  # "Uniform weights" in the authors' terms mean "no weighting" for us,
  # and "average marginals" means "marginal" with rcL
  # See ?rc for "marginals"
  unweighted <- rcL(color, nd=2, weighting="none",
                    layer.effect="heterogeneous", se="jackknife")
  marginal <- rcL(color, nd=2, weighting="marginal",
                  layer.effect="heterogeneous", se="jackknife")
  unweighted
  marginal

  # (our standard errors are much smaller for the marginal-weighted case)
  summary(unweighted)
  summary(marginal)

  opar <- par(mfrow=c(1, 2))
  plot(marginal, layer="Caithness", conf.ellipses=0.95)
  plot(marginal, layer="Aberdeen", conf.ellipses=0.95)
  par(opar)
  }


  ## Wong (2010), Table 4.7 (p. 103), model 9
  \dontrun{
  data(gss7590)

  model <- rcL(gss7590, nd=2, weighting="none", se="jackknife")

  model
  summary(model) # Jackknife standard errors are slightly different
                 # from their asymptotic counterparts

  # See ?plot.rcL for plotting
  }
}

\keyword{ models }
\keyword{ nonlinear }

\name{rcL.trans}
\alias{rcL.trans}
\title{Fitting Row-Column Association Models With Transitional Layer Effect}
\description{
  Fit log-multiplicative row-column association models with transitional layer effect,
  which are related to the RC(M)-L model, with one or several dimensions. Supported
  variants include (for square tables) symmetric (homogeneous) row and column scores,
  possibly combined with separate diagonal parameters.
}
\usage{
rcL.trans(tab, nd = 1,
          symmetric = FALSE,
          diagonal = c("none", "heterogeneous", "homogeneous"),
          weighting = c("marginal", "uniform", "none"),
          se = c("none", "jackknife", "bootstrap"),
          nreplicates = 100, ncpus = getOption("boot.ncpus"),
          family = poisson, weights = NULL,
          start = NULL, etastart = NULL, tolerance = 1e-8,
          iterMax = 5000, trace = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{tab}{ a three-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above three will be collapsed.}
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{symmetric}{ should row and column scores be constrained to be equal? Valid only for square tables.}
  \item{diagonal}{ what type of diagonal-specific parameters to include in the model, if any. This amounts to
    taking quasi-conditional independence, rather than conditional independence, as the baseline model.
    Valid only for square tables.}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{se}{ which method to use to compute standard errors for parameters.}
  \item{nreplicates}{ the number of bootstrap replicates, if enabled.}
  \item{ncpus}{ the number of processes to use for jackknife or bootstrap parallel computing. Defaults to
    the number of cores (see \code{\link[parallel]{detectCores}}), with a maximum of 5, but falls back to 1
    (no parallelization) if package \code{parallel} is not available.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{weights}{ an optional vector of weights to be used in the fitting process.}
  \item{start}{either \code{NA} to use optimal starting values, \code{NULL} to use
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
  \item{\dots}{ more arguments to be passed to \code{\link[gnm]{gnm}}}
}
\details{
  This function fits log-multiplicative row-column association models with regression-type layer effect which are
  \strong{experimental} models combining the principles behind RC(M)-L (Wong, 2010; see \code{\link{rcL}}) and regression-type
  models (Goodman & Hout, 1998). More specifically, like RC(M)-L models, row and column scores are allowed to vary across
  a layer variable, and the pattern of this variation follows the regression-type inspiration: for each dimension, a set of
  scores describes the first layer, another set describes the total variation of these scores need to describe the
  association observed for the last layer, and one parameter per layer describes the position of the layer between the
  first and the last layer. Compared with the RC(M)-L model with homogeneous scores across layers, this models allows
  for a finer description of changes since the ordering and distances of categories on a dimension are allowed to vary,
  and not only the general strength of the association. It is designed to describe transitions from one state to another,
  and is best suited for ordered layer variables like time (though the model is not sensitive to reordering of the layers).

  The general equation of the model is:
  \deqn{ log F_{ijk} = \lambda + \lambda^I_i + \lambda^J_j + \lambda^K_k
               + \lambda^{IK}_{ik} + \lambda^{JK}_{jk}
               + \sum_{m=1}^M { \phi_{mk} (\mu^S_{im} + \psi_{mk} \mu^V_{im}) (\nu^S_{jm} + \psi_{mk} \nu^V_{jm}) }}
  where \eqn{F_{ijk}} is the expected frequency for the cell at the intersection of row i, column j and layer k of
  \code{tab}, and M the number of dimensions. The \eqn{\psi_{mk}} parameter is constrained to be positive, equal to 0
  for the first layer (\eqn{m = 1}), and equal to 1 for the last layer.

  This model should not be confused with another combination of RC(M) models with the regression-type approach,
  presented by Goodman & Hout (1998:180), in which two separate RC(M) associations are used to describe respectively
  the stable and the varying components. In the present model, row and column scores for both components are summed
  \emph{before} entering the multiplicative interaction, which means only one RC(M) association exists.

  The returned object is a generic \code{rcL} association model describing the fitted scores for each layer. To analyze
  more specifically the variation of each (normalized) score from the first to the last layer, use:
  \code{model$assoc$row[,,dim(model$assoc$row)[3]] - model$assoc$row[,,1]}
  (and similarly for column scores).

  Actual model fitting is performed using \code{\link[gnm]{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife or bootstrap standard errors, and plotting.
  The default starting values are taken from a model with a stable RC(M) association (\dQuote{base model}). In some
  complex cases, using \code{start = NULL} to get random starting values can be more efficient, but it is also
  less stable and can converge to non-optimal solutions.
}
\value{
  A \code{rcL} object, with all the components of a \code{\link[gnm]{gnm}} object, plus an
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
  Goodman, L.A., and Hout, M. (1998). Statistical Methods and Graphical Displays for Analyzing How the Association
      Between Two Qualitative Variables Differs Among Countries, Among Groups, Or Over Time: A Modified Regression-Type
      Approach. \emph{Sociological Methodology} 28(1), 175-230.
  Wong, R.S-K. (2010). Association models. SAGE: Quantitative Applications in the Social Sciences.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.rcL}}, \code{\link[gnm]{gnm}}
}

\keyword{ models }
\keyword{ nonlinear }

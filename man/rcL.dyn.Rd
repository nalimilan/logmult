\name{rcL.dyn}
\alias{rcL.dyn}
\title{Fitting Row-Column Association Models With Regression-Type Layer Effect}
\description{
  \code{rcL.dyn} fits log-multiplicative row-column association models with regression-type,
  layer effect with one or several dimensions. Supported variants include (for square tables)
  symmetric (homogeneous) row and column scores, and separate diagonal parameters.
}
\usage{
  rcL.dyn(tab, nd = 1,
          symmetric = FALSE,
          diagonal = c("none", "heterogeneous", "homogeneous"),
          weighting = c("marginal", "uniform", "none"),
          std.err = c("none", "jackknife"), family = poisson,
          start = NA, tolerance = 1e-12, iterMax = 5000, trace = TRUE, ...)
}
\arguments{
  \item{tab}{ a three-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above three will be collapsed as appropriate.}
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{symmetric}{ should row and column scores be constrained to be equal? Valid only for square tables.}
  \item{diagonal}{ what type of diagonal-specific parameters to include in the model, if any. This amounts to
    taking quasi-conditional independence, rather than conditional independence, as the baseline model.
    Valid only for square tables.}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{std.err}{ whether to compute standard errors for scores. \dQuote{jackknife} is the only currently
     supported method, and is disabled by default because it increases computation time.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{start}{either \code{NA} (the default) to use reasonable starting values, \code{NULL} to use
     random starting values, or a vector of starting values for the parameters in the model;
     if a starting value is \code{NA}, the default starting value will be used.}
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
  This function fits log-multiplicative row-column association models with regression-type layer effect, which are
  and experimental model combining the principles behind RC(M)-L (Wong, 2010; see \code{\link{rcL}}) and regression-type
  models (Goodman & Hout, 1983). More specifically, like RC(M)-L models, row and column scores are allowed to vary accross
  a layer variable, and the pattern of this variation follows the regression-type inspiration: for each dimension, a set of
  scores describes the first layer, another set describes the total variation of these scores need to describe the
  association observed for the last layer, and one parameter per layer describes the position of the layer between the
  first and the last layer. Compared with the RC(M)-L model with homogeneous scores accross layers, this models allows
  for a finer description of changes since the ordering and distances of categories on a dimension are allowed to vary,
  and not only the general strength of the association. It is designed to describe transitions from one state to another,
  and is best suited for ordered layer variables like time (though the model is not sensitive to reordering of the layers).

  The general equation of the model is:
  \deqn{ log F_{ijk} = \lambda + \lambda^I_i + \lambda^J_j + \lambda^K_k
               + \sum_{m=1}^M { \phi_{km} (\mu^S_{im} + \psi_km \mu^V_{im}) (\nu^S_{jm} + \psi_km \nu^V{jm}) }}
  where \eqn{F_ijk} is the frequency observed in the cell at the intersection of row i, column j and layer k of
  \code{tab}, and M the number of dimensions. The \eqn{\psi_km} parameter is constrained to be positive, equal to 0
  for the first layer (\eqn{m = 1}), and equal to 1 for the last layer.

  This model should not be confused with another combination of RC(M) models with the regression-type approach,
  presented by Goodman & Hout (1998:180), in which two separate RC(M) associations are used to describe respectively
  the stable and the varying components. In the present model, row and column scores for both components are summed
  \emph{before} entering the multiplicative interaction, which means only one RC(M) association exists.

  The returned object is a generic rcL association describing the fitted scores for each layer. To analyze more
  specifically the variation of each (normalized) score from the first to the last layer, use:
  \code{model$assoc$row[,,dim(model$assoc$row)[3]] - model$assoc$row[,,1]}
  (and similarly for column scores).

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife standard errors, and plotting.
  The default starting values correspond to the parameters computed for the model without association parameters
  (\dQuote{base model}); association parameters start with random starting values. In some complex cases,
  using \code{start = NA} to get completely random starting values can be more efficient, but it is also
  less stable.
}
\value{
  A \code{rcL} object, with all the components of a \code{\link{gnm}} object, plus an
    \code{assoc.rcL} component holding the most relevant association information:
  \item{phi }{The intrisic association parameters, one per dimension and per layer.}
  \item{row }{Row scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{weighting }{The name of the weighting method used, reflected by \code{row.weights}
    and \code{col.weights}.}
  \item{row.weights }{The row weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{col.weights }{The column weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{covmat }{The variance-covariance matrix of the parameter's estimators, if \code{std.err}
    was set. It contains first phi coefficients, then normalized row and column scores, then
    sqrt(phi) * row/column scores.}
  \item{covtype }{The method used to compute the variance-covariance matrix. Currently, only
    \dQuote{jackknife} is supported.}
}
\note{
  Default \code{tolerance} value is very low to ensure the model converges to the right point and that
  coefficients are stable even in most complex cases. This comes at a cost, in particular with high-dimensional
  tables and/or when computing jackknife standard errors. If you know what you are doing, raising \code{tolerance}
  to 1e-6 or an intermediate value can lead to substantial computing time reductions; in that case, be sure to run the
  model several times to check for numerical instability of the scores. (Very large jackknife standard errors are
  a good sign that proper convergence was not attained at the chosen tolerance level.)
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
  \code{\link{plot.rcL}}, \code{\link{gnm}}
}
\examples{

}

\keyword{ models }
\keyword{ nonlinear }

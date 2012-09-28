\name{hmskewL}
\alias{hmskewL}
\alias{print.hmskewL}
\title{Fitting van der Heijen & Mooijaart Skew-Symmetric Association Model With Layer Effect}
\description{
  \code{hmskewL} fits an extension of the skew-symmetric association model proposed in
  van der Heijden & Mooijaart (1995) to describe asymmetry of square tables.
  This model introduces a layer effect by which the strength of skew-symmetric association,
  and optionnally scores, can vary over the levels of the third dimension of the table.
  It can be combined with symmetric association models like a quasi-symmetry
  (the default) or symmetric (homogeneous) RC(M)-L models, with or without layer effect.
}
\usage{
hmskewL(tab, nd.symm = NA,
        layer.effect.skew = c("homogeneous.scores", "heterogeneous", "none"),
        layer.effect.symm = c("uniform", "homogeneous.scores", "heterogeneous", "none"),
        diagonal = c("none", "heterogeneous", "homogeneous"),
        weighting = c("marginal", "uniform", "none"),
        se = c("none", "jackknife", "bootstrap"),
        nreplicates = 100, ncpus = getOption("boot.ncpus"),
        family = poisson, weights = NULL,
        start = NA, etastart = NULL, tolerance = 1e-6,
        iterMax = 5000, trace = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{tab}{ a three-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above three will be collapsed as appropriate. First two dimensions must be symmetric
    (i.e. of the same length).}
  \item{nd.symm}{ the number of dimensions to include in the \emph{symmetric} RC(M) association. Cannot exceed
    \code{2 * min(nrow(tab) - 1, ncol(tab) - 1)} (the model is already saturated at that point).
    If \code{NA}, a quasi-symmetric association is used instead of a RC(M) model.}
  \item{layer.effect.skew}{ determines the form of the interaction between skew-symmetric association and layers.
    See \dQuote{Details} below.}
  \item{layer.effect.symm}{ determines the form of the interaction between symmetric row-column association,
     or quasi-symmetric association (if \code{nd.symm = NA}) and layers. See \dQuote{Details} below.}
  \item{diagonal}{ what type of diagonal-specific parameters to include in the model, if any. Only makes sense
     when \code{nd.symm} is not \code{NA} (else, diagonal parameters are already included).}
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
  \item{verbose}{ a logical value indicating whether progress indicators should be printed,
     including a diagnostic error message if the algorithm restarts.}
  \item{\dots}{ more arguments to be passed to \code{\link{gnm}}}
}
\details{
  This models follows a slightly modified equation similar to that of the original model presented by
  van der Heijden & Mooijaart (1995) (see \code{\link{hmskew}}).

  These models follow the equation:
  \deqn{ log F_{ijk} = q_{ijk} + \phi_k (\nu_{ik} \mu_{jk} - \mu_{jk} \nu_{ik}) }
  where \eqn{F_ijk} is the frequency observed in the cell at the intersection of row i, column j and layer k of
  \code{tab}, and \eqn{q_{ijk}} a quasi-symmetric or a RC(M)-L association.

  If \code{layer.effect.skew} is set to \sQuote{heterogeneous}, different scores will be computed for each level,
  which is equivalent to fitting separate models using \code{\link{hmskew}} on the k two-way tables.
  If it is set to \sQuote{homogeneous.scores}, then \eqn{\mu_{ik} = \mu_i} and \eqn{\nu_{ik} = \nu_i} for all
  layers k: only the \eqn{\phi_k} are allowed to vary accross layers. If it is set to \sQuote{none}, then in addition
  to the previous conditions all \eqn{\phi_{km}} are forced to be equal for all layers k, which amounts to a stability
  of the association accross layers.

  When \code{nd.symm} is different from \code{NA}, the symmetric association works exactly like a call to \code{rcL},
  with parameters \code{nd.symm} and \code{layer.effect.symm} translated respectively to \code{nd} and
  \code{layer.effect}. When \code{nd.symm == NA}, symmetric association parameters are either stable accross layers,
  are multiplied by a layer coefficient (cf. UNIDIFF model), or are different for each layer, when
  \code{layer.effect.symm} is respectively \code{none}, \code{uniform} and \code{heterogeneous}.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife or bootstrap standard errors, and plotting.
  The default starting values correspond to the parameters computed for the model without association parameters
  (\dQuote{base model}); association parameters start with random starting values. In some complex cases,
  using \code{start = NA} to get completely random starting values can be more efficient, but it is also
  less stable.
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
  van der Heijden, P.G.M., and Mooijaart, A. (1995). Some new log bilinear models for the analysis
    of asymmetry in a square contingency table. \emph{Sociol. Methods and Research} 24, 7-29.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.hmskewL}}, \code{\link{hmskew}}, \code{\link{gnm}}
}
\keyword{ models }
\keyword{ nonlinear }

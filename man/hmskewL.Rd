\name{hmskewL}
\alias{hmskewL}
\alias{print.hmskewL}
\title{Fitting van der Heijen & Mooijaart Skew-Symmetric Association Model With Layer Effect}
\description{
  Fits an extension of the skew-symmetric association model proposed in van der Heijden & Mooijaart (1995)
  to describe asymmetry of square tables. This model introduces a layer effect by which the strength of
  skew-symmetric association, and optionnally scores, can vary over the levels of the third dimension of the table.
  Skew-symmetric association can be combined with quasi-symmetry (the default), quasi-independence, or
  symmetric (homogeneous) RC(M) associations, with or without layer effect.
}
\usage{
hmskewL(tab, nd.symm = NA,
        layer.effect.skew = c("homogeneous.scores", "heterogeneous",
                              "none"),
        layer.effect.symm = c("heterogeneous", "uniform",
                              "homogeneous.scores", "none"),
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
    if present, dimensions above three will be collapsed. First two dimensions must be symmetric
    (i.e. of the same length).}
  \item{nd.symm}{ the number of dimensions to include in the \emph{symmetric} RC(M) association. Cannot exceed
    \code{2 * min(nrow(tab) - 1, ncol(tab) - 1)} (the model is already saturated at that point).
    If \code{NA} (the default), a full quasi-symmetric association is used instead of a RC(M) model; if \code{0},
    quasi-independence is used.}
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
  \item{\dots}{ more arguments to be passed to \code{\link{gnm}}}
}
\details{
  This model follows an equation inspired from that presented by van der Heijden & Mooijaart (1995) for two-way tables
  (see \code{\link{hmskew}}):
  \deqn{ log F_{ijk} = q_{ijk} + \phi_k (\nu_{ik} \mu_{jk} - \mu_{jk} \nu_{ik}) }
  where \eqn{F_{ijk}} is the expected frequency for the cell at the intersection of row i, column j and layer k of
  \code{tab}, and \eqn{q_{ij}} a quasi-symmetric specification, with either full interaction parameters, or
  a RC(M) association. See reference for detailed information about the degrees of freedom and the identification
  constraints applied to the scores.

  If \code{layer.effect.skew} is set to \sQuote{heterogeneous}, different scores will be computed for each level,
  which is equivalent to fitting separate models using \code{\link{hmskew}} on the k two-way tables.
  If it is set to \sQuote{homogeneous.scores}, then \eqn{\mu_{ik} = \mu_i} and \eqn{\nu_{ik} = \nu_i} for all
  layers k: only the \eqn{\phi_k} are allowed to vary across layers. If it is set to \sQuote{none}, then in addition
  to the previous conditions all \eqn{\phi_{km}} are forced to be equal for all layers k, which amounts to a stability
  of the association across layers.

  When \code{nd.symm} is different from \code{NA}, the symmetric association works exactly like a call to \code{rcL},
  with parameters \code{nd.symm} and \code{layer.effect.symm} translated respectively to \code{nd} and
  \code{layer.effect}. When \code{nd.symm == NA}, symmetric association parameters are either stable across layers,
  are multiplied by a layer coefficient (UNIDIFF model, see \code{\link{unidiff}}), or are different for each layer,
  when \code{layer.effect.symm} is respectively \code{none}, \code{uniform} and \code{heterogeneous}.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife or bootstrap standard errors, and plotting.
  The default starting values for skew association parameters are computed using an eigen value decomposition from the
  results of the model without skew association component (\dQuote{base model}); if \code{nd.symm} is not \code{NA} and
  strictly positive, random starting values are used. In some complex cases, using \code{start = NULL} to start with
  random values can be more efficient, but it is also less stable and can converge to non-optimal solutions.
}
\value{
  A \code{hmskewL} object, which is a subclass of an \code{rcL.symm} object (see \code{\link{rcL}}) if
    \code{nd.symm} is strictly positive. In addition to this class, it contains a \code{assoc.hmskew} component
    holding information about the \emph{skew-symmetric} association:
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

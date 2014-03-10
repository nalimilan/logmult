\name{svyassocmod}
\alias{svyrc}
\alias{svyhmskew}
\alias{svyyrcskew}
\alias{svyrcL}
\alias{svyrcL.trans}
\alias{svyhmskewL}
\title{Fitting Association Models With Complex Survey Data}
\description{
     Fit association models to data from a complex survey
     design, with inverse-probability weighting and (optionally)
     standard errors based on replicate weights.
}
\usage{
svyrc(formula, design, nd = 1,
      symmetric = FALSE, diagonal = FALSE,
      weighting = c("marginal", "uniform", "none"),
      rowsup = NULL, colsup = NULL,
      Ntotal = nrow(design), exclude = c(NA, NaN),
      se = c("none", "replicate"),
      ncpus = getOption("boot.ncpus"),
      family = quasipoisson, weights = NULL,
      start = NULL, etastart = NULL, tolerance = 1e-8,
      iterMax = 5000, trace = FALSE, verbose = TRUE, ...)

svyhmskew(formula, design, nd.symm = NA, diagonal = FALSE,
          weighting = c("marginal", "uniform", "none"),
          rowsup = NULL, colsup = NULL,
          Ntotal = nrow(design), exclude = c(NA, NaN),
          se = c("none", "replicate"),
          ncpus = getOption("boot.ncpus"),
          family = quasipoisson, weights = NULL,
          start = NULL, etastart = NULL, tolerance = 1e-8,
          iterMax = 5000, trace = FALSE, verbose = TRUE, ...)

svyyrcskew(formula, design, nd.symm = NA, nd.skew = 1, diagonal = FALSE,
           weighting = c("marginal", "uniform", "none"),
           rowsup = NULL, colsup = NULL,
           Ntotal = nrow(design), exclude = c(NA, NaN),
           se = c("none", "replicate"),
           ncpus = getOption("boot.ncpus"),
           family = quasipoisson, weights = NULL,
           start = NA, etastart = NULL, tolerance = 1e-8,
           iterMax = 15000, trace = FALSE, verbose = TRUE, ...)

svyrcL(formula, design, nd = 1,
       layer.effect = c("homogeneous.scores",
                        "heterogeneous", "none"),
       symmetric = FALSE,
       diagonal = c("none", "heterogeneous", "homogeneous"),
       weighting = c("marginal", "uniform", "none"),
       Ntotal = nrow(design), exclude = c(NA, NaN),
       se = c("none", "replicate"),
       ncpus = getOption("boot.ncpus"),
       family = quasipoisson, weights = NULL,
       start = NULL, etastart = NULL, tolerance = 1e-8,
       iterMax = 5000, trace = FALSE, verbose = TRUE, ...)

svyrcL.trans(formula, design, nd = 1,
             symmetric = FALSE,
             diagonal = c("none", "heterogeneous", "homogeneous"),
             weighting = c("marginal", "uniform", "none"),
             Ntotal = nrow(design), exclude = c(NA, NaN),
             se = c("none", "replicate"),
             ncpus = getOption("boot.ncpus"),
             family = quasipoisson, weights = NULL,
             start = NULL, etastart = NULL, tolerance = 1e-8,
             iterMax = 5000, trace = FALSE, verbose = TRUE, ...)

svyhmskewL(formula, design, nd.symm = NA,
           layer.effect.skew = c("homogeneous.scores", "heterogeneous",
                                 "none"),
           layer.effect.symm = c("heterogeneous", "uniform",
                                 "homogeneous.scores", "none"),
           diagonal = c("none", "heterogeneous", "homogeneous"),
           weighting = c("marginal", "uniform", "none"),
           Ntotal = nrow(design), exclude = c(NA, NaN),
           se = c("none", "replicate"),
           ncpus = getOption("boot.ncpus"),
           family = quasipoisson, weights = NULL,
           start = NULL, etastart = NULL, tolerance = 1e-8,
           iterMax = 5000, trace = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{formula}{ a formula specifying margins for the table (using \sQuote{+} only)
                  on which the model will be fitted (passed to
                  \code{\link[survey]{svytable}}); dimensions of the
                  resulting table must match the models expectations.}
  \item{design}{ a survey object; if \code{se == "replicate"}, must be of class
                 \code{\link[survey]{svrepdesign}} (see \dQuote{Details} below).}
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{nd.symm}{ the number of dimensions to include in the \emph{symmetric} RC(M) association. Cannot exceed
    \code{2 * min(nrow(tab) - 1, ncol(tab) - 1)} (quasi-symmetry model).
    If \code{NA} (the default), a full quasi-symmetric association is used instead of a RC(M) model; if \code{0},
    quasi-independence is used.}
  \item{nd.skew}{ the number of dimensions to include in the \emph{skew-symmetric} RC(M) association.}
  \item{layer.effect}{ determines the form of the interaction between row-column association and layers.
     See \dQuote{Details} below.}
  \item{layer.effect.skew}{ determines the form of the interaction between skew-symmetric association and layers.
    See \dQuote{Details} below.}
  \item{layer.effect.symm}{ determines the form of the interaction between symmetric row-column association,
     or quasi-symmetric association (if \code{nd.symm = NA}) and layers. See \dQuote{Details} below.}
  \item{symmetric}{ should row and column scores be constrained to be equal? Valid only for square tables.}
  \item{diagonal}{ what type of diagonal-specific parameters to include in the model, if any. Only makes sense
     when \code{nd.symm} is not \code{NA} (else, diagonal parameters are already included).}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{Ntotal}{ sum of counts to normalize the table to (passed to \code{\link[survey]{svytable}}).
      See \dQuote{Details} below..}
  \item{exclude}{ a vector of values to be exclude when building the table, passed to \code{\link{xtabs}}.}
  \item{rowsup}{ if present, a matrix with the same columns as \code{tab} and rows corresponding to the columns
                of \code{colsup}, giving supplementary (passive) rows.}
  \item{colsup}{ if present, a matrix with the same rows as \code{tab} and columns corresponding to the rows of
                \code{colsup}, giving supplementary (passive) columns.}
  \item{se}{whether to compute replicate standard errors or not (only supported for \code{svrepdesign} objects).}
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
     The model is fitted to a table with probabilities estimated by
     \code{\link[survey]{svytable}} and (when \code{Ntotal = nrow(design)})
     with the sample size equal to the observed sample size, treating the
     resulting table as if it came from iid multinomial sampling, as described
     by Rao and Scott. This assumption affects the fit statistics but not parameter
     point estimates.

     Standard errors that do not rely on this assumption can be computed
     by fitting the model using each series of replicate weights. If your
     data does not come with replicate weights, use
     \code{\link[survey]{as.svrepdesign}} to create them first, and pass
     the resulting \code{\link[survey]{svrepdesign}} object via the
     \code{design} argument.
}
\value{
  An \code{assocmod} object whose exact class depends on the function called.
}
\note{
  Note that printed fit statistics and degrees of freedom rely on the
  iid assumption. This is also the case of the variance-covariance matrix
  returned by the \code{\link[gnm]{vcov.gnm}} function.
}
\references{
  Rao, J.N.K., Scott, A.J. (1984). On Chi-squared Tests For Multiway
  Contingency Tables with Proportions Estimated From Survey Data.
  \emph{Annals of Statistics} 12, 46-60.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{rc}}, \code{\link{hmskew}}, \code{\link{yrcskew}}, \code{\link{rcL}},
  \code{\link{rcL.trans}}, \code{\link{hmskewL}}

  \code{\link[survey]{svytable}}, \code{\link[survey]{svyloglin}},
  \code{\link[survey]{svyglm}}, \code{\link[survey]{as.svrepdesign}}
}
\keyword{ models }
\keyword{ nonlinear }

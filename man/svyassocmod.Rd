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
svyrc(formula, design,
      Ntotal = nrow(design),
      exclude = c(NA, NaN),
      se = c("none", "replicate"), ...)

svyhmskew(formula, design,
          Ntotal = nrow(design),
          exclude = c(NA, NaN),
          se = c("none", "replicate"), ...)

svyyrcskew(formula, design,
           Ntotal = nrow(design),
           exclude = c(NA, NaN),
           se = c("none", "replicate"), ...)

svyrcL(formula, design,
       Ntotal = nrow(design),
       exclude = c(NA, NaN),
       se = c("none", "replicate"), ...)

svyrcL.trans(formula, design,
             Ntotal = nrow(design),
             exclude = c(NA, NaN),
             se = c("none", "replicate"), ...)

svyhmskewL(formula, design,
           Ntotal = nrow(design),
           exclude = c(NA, NaN),
           se = c("none", "replicate"), ...)
}
\arguments{
  \item{formula}{ a formula specifying margins for the table (using \sQuote{+} only)
                  on which the model will be fitted (passed to
                  \code{\link[survey]{svytable}}); dimensions of the
                  resulting table must match the models expectations.}
  \item{design}{ a survey object; if \code{se == "replicate"}, must be of class
                 \code{\link[survey]{svrepdesign}} (see \dQuote{Details} below).}
  \item{Ntotal}{ sum of counts to normalize the table to (passed to \code{\link[survey]{svytable}}).
      See \dQuote{Details} below.}
  \item{exclude}{ a vector of values to be exclude when building the table, passed to \code{\link{xtabs}}.}
  \item{se}{set to \dQuote{replicate} to compute standard errors based on replicate weights
      (see \dQuote{Details} below and \code{\link{se.assoc}}).}
  \item{...}{Other arguments passed to the corresponding functions: see \code{\link{rc}},
      \code{\link{hmskew}}, \code{\link{yrcskew}}, \code{\link{rcL}},
      \code{\link{rcL.trans}} or \code{\link{hmskewL}}.}
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

\name{se.assoc}
\alias{se}
\alias{se.assoc}
\alias{se.rc}
\alias{se.hmskew}
\alias{se.yrcskew}
\alias{se.rcL}
\alias{se.hmskewL}
\title{Standard Errors for Association Models}
\description{Get standard errors for log-multiplicative association scores and intrinsic
             association coefficients.}
\usage{
se(x, ...)

\method{se}{assoc}(x, type = c("se", "quasi.se"), ...)

\method{se}{rc}(x, type = c("se", "quasi.se"), ...)

\method{se}{hmskew}(x, type = c("se", "quasi.se"), ...)

\method{se}{yrcskew}(x, type = c("se", "quasi.se"), ...)

\method{se}{rcL}(x, type = c("se", "quasi.se"), ...)
}
\arguments{
  \item{x}{an \code{\link{assoc}} object with a non-null \code{covmat} component (for \code{se.assoc});.
           or a \code{\link{rc}}, \code{\link{hmskew}}, \code{\link{hmskewL}}, \code{\link{yrcskew}},
           \code{\link{rcL}} or \code{\link{rcL.trans}} object fitted with the \code{se} argument different
           from \dQuote{none} (for other functions).}
  \item{type}{the type of standard errors to be computed (see \dQuote{Details} below).}
  \item{\dots}{currently unused.}
}
\details{
  Currently, only jackknife or bootstrap standard errors are supported, depending on the \code{se} argument
  passed when fitting the model. \strong{Some care is needed before using such standard errors and confidence
  intervals}. First one must ensure all model replicates converged to a correct solution, especially for bootstrap;
  second, when relying on normal confidence intervals computed from these standard errors, one must ensure that the
  coefficients estimators follow a normal distribution. Both checks can be performed by calling
  \code{\link[boot]{plot.boot}} on the \code{boot.results} component of the \code{assoc} object of the models (not
  supported for jackknife), with the \code{index} argument identifying the coefficient of interest (call \code{colnames}
  on the \code{t} member of the \code{boot.results} object to find out the index you need).

  If outliers are present, standard errors and confidence intervals will be artificially large; to fix this, the
  \code{tolerance} argument must be set to a smaller value when fitting the models (which may in turn require
  increasing the value of the \code{iterMax} argument if convergence is too slow. Once outliers are removed, if
  coefficient estimates are still not normally distributed, robust bootstrap confidence intervals can be computed
  using \code{\link[boot]{boot.ci}} on the same object, provided a large number of replicates (> 1000) were computed.

  For each replicate, stable scores and intrinsic association coefficients are identified using an orthogonal
  procrustes analysis to suppress meaningless variations due to random reflections, permutations and rotations of
  dimensions (Milan & Whittaker, 1995). For \code{\link{hmskew}} and \code{\link{hmskewL}} models, a rotation of
  each pair of dimensions and a permutation of pairs of dimensions is performed, but no reflection as it would change
  the sign of intrinsic association coefficients.

  Quasi-standard errors are computed using \code{\link[qvcalc]{qvcalc}}. See the help page for this function for
  details and references about them.
}
\value{
  An object of the same form as the \code{assoc} component of the model, but with standard errors rather than
  the corresponding coefficients.
}
\author{
  Milan Bouchet-Valat
}
\references{
  Milan, L., and J. Whittaker (1995). Application of the Parametric Bootstrap to Models that Incorporate a Singular
  Value Decomposition. \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)} 44(1), 31-49.
}
\seealso{\code{\link{assoc}}, \code{\link{rc}}, \code{\link{hmskew}}, \code{\link{hmskewL}}, \code{\link{yrcskew}},
         \code{\link{rcL}}, \code{\link{rcL.trans}}}
\examples{
  # See ?rc about Wong (2010)
}


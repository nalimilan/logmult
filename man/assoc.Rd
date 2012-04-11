\name{assoc}
\alias{assoc}
\alias{assoc.rc}
\alias{assoc.rc.homog}
\alias{assoc.hmskew}
\alias{assoc.yrcskew}
\alias{assoc.yrcskew.homog}
\alias{assoc.rcL}
\alias{assoc.rcL.homog}
\title{Identifying Scores from Association Models}
\description{Identify log-multiplicative association scores from over-parameterized gnm models}
\usage{\method{assoc}{rc}(model, weights = c("marginal", "uniform", "none"), ...)
       \method{assoc}{rc.homog}(model, weights = c("marginal", "uniform", "none"), ...)
       \method{assoc}{hmskew}(model, weights = c("marginal", "uniform", "none"), ...)
       \method{assoc}{yrcskew}(model, weights = c("marginal", "uniform", "none"), ...)
}
\arguments{
  \item{model}{a \code{\link{gnm}} object, usually obtained using \code{\link{rc}}, \code{\link{hmskew}} or
               \code{\link{yrcskew}}, but not necessarily.}
  \item{weights}{the weights to be used when normalizing scores (see \sQuote{Details}).}
  \item{...}{currently unused.}
}
\details{
  These functions extract parameters from \code{\link{gnm}} log-multiplicative models
  and make them identifiable by imposing the required constraints on them. The general pattern
  is that row and column scores are separately centered around 0 and scaled so that they sum to 1,
  and so that their cross-dimensional correlation is null. From this operation result two series of
  scores (rows and columns) plus an intrinsic association coefficient (phi) for each dimension.

  Most users do not need to call these directly, but they are still made public since they may be
  useful for advanced uses, notably when combining log-multiplicative association components with
  other model specifications. \code{assoc} can be used to identify the scores, the rest of the
  coefficients being extracted manually by the caller.
}
\seealso{\code{\link{rc}}, \code{\link{hmskew}}, \code{\link{yrcskew}}}

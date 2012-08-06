\name{assoc}
\alias{assoc}
\alias{assoc.rc}
\alias{assoc.rc.symm}
\alias{assoc.hmskew}
\alias{assoc.yrcskew}
\alias{assoc.yrcskew.symm}
\alias{assoc.rcL}
\alias{assoc.rcL.symm}
\alias{assoc.rcL.trans}
\alias{assoc.rcL.trans.symm}
\title{Identifying Scores from Association Models}
\description{Identify log-multiplicative association scores from over-parameterized gnm models}
\usage{\method{assoc}{rc}(model, weighting = c("marginal", "uniform", "none"), ...)
       \method{assoc}{rc.symm}(model, weighting = c("marginal", "uniform", "none"), ...)
       \method{assoc}{hmskew}(model, weighting = c("marginal", "uniform", "none"), ...)
       \method{assoc}{yrcskew}(model, weighting = c("marginal", "uniform", "none"), ...)
       \method{assoc}{rcL}(model, weighting = c("marginal", "uniform", "none"), ...)
       \method{assoc}{rcL.symm}(model, weighting = c("marginal", "uniform", "none"), ...)
       \method{assoc}{rcL.trans}(model, weighting = c("marginal", "uniform", "none"), ...)
       \method{assoc}{rcL.trans.symm}(model, weighting = c("marginal", "uniform", "none"), ...)
}
\arguments{
  \item{model}{a \code{\link{gnm}} object, usually obtained using \code{\link{rc}}, \code{\link{hmskew}},
               \code{\link{yrcskew}}, \code{\link{rcL}}, or \code{\link{rcL.trans}}, but not necessarily.}
  \item{weighting}{the weights to be used when normalizing scores (see \sQuote{Details}).}
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
\value{
  An \code{assoc} object with the following components:
  \item{phi }{The intrisic association parameters, one per dimension.}
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
  \item{covmat }{The variance-covariance matrix of the parameter's estimators, if \code{se}
    was set. It contains first phi coefficients, then normalized row and column scores, then
    sqrt(phi) * row/column scores.}
  \item{covtype }{The method used to compute the variance-covariance matrix. Currently, only
    \dQuote{jackknife} is supported.}
}
\seealso{\code{\link{rc}}, \code{\link{hmskew}}, \code{\link{yrcskew}}, \code{\link{rcL}},
         \code{\link{rcL.trans}} }

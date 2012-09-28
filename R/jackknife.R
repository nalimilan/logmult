# Vaguely adapted from bootstrap package by Rob Tibshirani (R port by Friedrich Leisch)
# License: BSD
#
# Three main modifications:
# 1) Estimate the model only once for each cell of the table,
# and compute a mean weighted by cell frequencies at the end
# (Wong, Association models, 2010, p. 28-29; Clogg & Shihadeh, 1994, p. 34-38)
# 2) theta must return a vector rather than a single number,
# to allow handling a series of parameters in the same run
# 3) The variance-covariance matrix is returned, instead of standard errors
#
# Formula taken from original function and adapted to work with table cell weights.
#
# Objects needed by theta() on the cluster nodes can be passed as arguments and handled
# in theta(), or simply passed as *named* arguments to jackknife(): they will be exported
# the the nodes' environments. 
jackknife <- function(x, theta, ..., w=rep(1, length(x)), cl=NULL)
{
    call <- match.call()
    stopifnot(length(w) == length(x))
    w <- as.numeric(w)
    stopifnot(all(w >= 0))
    n <- length(x)
    u <- vector("list", n)

    # Run this first to find out caller errors before running parLapply
    thetahat <- as.numeric(theta(x, ...))

    if(!is.null(cl) && require(parallel)) {
        dots <- list(...)
        dotsnames <- names(dots)

        if(length(dotsnames[dotsnames != ""]) > 0)
            parallel::clusterExport(cl, dotsnames[dotsnames != ""], as.environment(dots[dotsnames != ""]))

        u <- parallel::parLapply(cl, 1:n, function(i, x, theta, ...) theta(x[-i], ...), x, theta, ...)
    }
    else if(!is.null(cl) && require(snow)) {
        dots <- list(...)
        dotsnames <- names(dots)

        if(length(dotsnames[dotsnames != ""]) > 0) {
            # clusterExport() in snow 0.3-3 does not take an 'envir' argument
#             envir <- as.environment(dots[dotsnames != ""])
#             environment(clusterExport) <- envir
            snow::clusterExport(cl, dotsnames[dotsnames != ""])
        }

        u <- snow::clusterApply(cl, 1:n, function(i, x, theta, ...) theta(x[-i], ...), x, theta, ...)
    }
    else {
        for(i in 1:n) {
            u[[i]] <- as.numeric(theta(x[-i], ...))
        }
    }

    u <- do.call(rbind, u)

    # Remove replicates with NAs when computing statistics
    # This can be used by theta() to skip a failed replicate
    u2 <- u[rowSums(is.na(u)) == 0,]

    tot <- sum(w)
    mean.u <- colSums(sweep(u2, 1, w, "*"))/tot
    jack.bias <- (tot - 1) * (mean.u - thetahat)
    dev.u <- sweep(u2, 2, mean.u, "-")

    list(dev = dev.u, bias = jack.bias, values = u)
}


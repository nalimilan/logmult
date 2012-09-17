# Vaguely adapted from bootstrap package, with three modifications:
# 1) Estimate the model only once for each cell of the table,
# and compute a mean weighted by cell frequencies at the end
# (Wong, Association models, 2010, p. 28-29; Clogg & Shihadeh, 1994, p. 34-38)
# 2) theta must return a vector rather than a single number,
# to allow handling a series of parameters in the same run
# 3) The variance-covariance matrix is returned, instead of standard errors
#
# Formula taken from original function and from e.g. Steinhauer & Wuergler, 2010,
# "Leverage and Covariance Matrix Estimation in Finite-Sample IV Regressions", p. 11.
# (adapted to work with table cell weights).
#
# Objects needed by theta() on the cluster nodes can be passed as arguments and handled
# in theta(), or simply passed as *named* arguments to jackknife(): they will be exported
# the the nodes' environments. 
jackknife <- function(x, theta, ..., w=rep(1, length(x)), ncpus=1)
{
    call <- match.call()
    stopifnot(length(w) == length(x))
    w <- as.numeric(w)
    stopifnot(all(w >= 0))
    n <- length(x)
    u <- vector("list", n)

    # Run this first to find out caller errors before running parLapply
    thetahat <- as.numeric(theta(x, ...))

    if(ncpus > 1 && require(parallel)) {
        cl <- makeCluster(ncpus)
        on.exit(stopCluster(cl))

        dots <- list(...)
        dotsnames <- names(dots)

        if(length(dotsnames[dotsnames != ""]) > 0)
            clusterExport(cl, dotsnames[dotsnames != ""], as.environment(dots[dotsnames != ""]))

        u <- parLapply(cl, 1:n, function(i, x, theta, ...) theta(x[-i], ...), x, theta, ...)
    }
    else if(ncpus > 1 && require(snow)) {
        cl <- snow::makeCluster(rep("localhost", ncpus), type="SOCK")
        on.exit(snow::stopCluster(cl))

        dots <- list(...)
        dotsnames <- names(dots)

        if(length(dotsnames[dotsnames != ""]) > 0) {
            # clusterExport() in snow 0.3-3 does not take an 'envir' argument
#             envir <- as.environment(dots[dotsnames != ""])
#             environment(clusterExport) <- envir
            clusterExport(cl, dotsnames[dotsnames != ""])
        }

        u <- snow::clusterApply(cl, 1:n, function(i, x, theta, ...) theta(x[-i], ...), x, theta, ...)
    }
    else {
        for(i in 1:n) {
            u[[i]] <- as.numeric(theta(x[-i], ...))
        }
    }
    u <- do.call(cbind, u)
    # Remove replicates with NAs when computing statistics
    # This can be used by theta() to skip a failed replicate
    u2 <- u[, colSums(is.na(u)) == 0]
    tot <- sum(w)
    mean.u <- rowSums(sweep(u2, 2, w, "*"))/tot
    jack.bias <- (tot - 1) * (mean.u - thetahat)
    dev.u <- sweep(u2, 1, mean.u, "-")
    jack.vcov <- (tot - 1)/tot * sweep(dev.u, 2, w, "*") %*% t(dev.u)
    return(list(jack.vcov = jack.vcov, jack.bias = jack.bias, jack.values = u, 
                call = call))
}

# Additional arguments are needed so that update() finds them even when using parLapply
jackknife.assoc <- function(x, model, repl.verbose=FALSE, ...) {
  tab <- model$data

  if(repl.verbose) {
      iter <- which(!1:length(tab) %in% x)
      if(length(iter) == 1)
          cat(sprintf("Iteration for cell %i of %i\n", iter, length(tab)))
      else
          cat("Initial iteration\n")
  }

  if(sum(tab[-x]) > 0) {
      mat <- tab
      mat[] <- -1
      mat[x] <- 0
      tab <- tab + mat
  }

  replicate.assoc(model, tab, repl.verbose=repl.verbose, ...)
}

boot.assoc <- function(data, indices, args) {
  tab <- args$model$data

  # Create a table not including the observations identified by indices
  cs <- cumsum(tab)
  for(i in setdiff(1:sum(tab), indices)) {
      index <- min(which(i <= cs))
      tab[index] <- tab[index] - 1
  }


  # We need to pass all arguments through "args" to prevent them
  # from being catched by boot(), especially "weights"
  args$tab <- tab
  do.call(replicate.assoc, args)
}

# Replicate model with new data, and combine assoc components into a vector
replicate.assoc <- function(model.orig, tab, assoc1, assoc2, weighting, ...,
                            base=NULL, repl.verbose=FALSE) {
  library(assoc)

  # Remove warnings because we handle this below
  suppressWarnings(model <- update(model.orig, tab=tab,
                                   start=parameters(model.orig), etastart=as.numeric(predict(model.orig)),
                                   verbose=repl.verbose, trace=repl.verbose, se="none"))

  if(!model$converged) {
      cat("Model replicate did not converge.\nData was:\n")
      print(tab)
      cat(sprintf("Trying again with different starting values...\n", model$iterMax))

      suppressWarnings(model <- update(model, start=NA, etastart=NULL,
                                       verbose=TRUE, trace=TRUE, se="none"))

  }

  if(!model$converged) {
      cat(sprintf("Trying once again with random starting values...\n", model$iterMax))

  ass1 <- assoc1(model, weighting=weighting)
  ret <- c(t(ass1$phi))

  if(dim(ass1$row)[3] == 1) {
      # Even if the scores are the same for all layers, we replicate them for simplicity's sake
      for(i in 1:nrow(ass1$phi))
           ret <- c(ret, ass1$row[,,1], ass1$col[,,1],
                    sweep(ass1$row[,,1, drop=FALSE], 2, sqrt(abs(ass1$phi[i,])) * sign(ass1$phi[i,]), "*"),
                    sweep(ass1$col[,,1, drop=FALSE], 2, sqrt(abs(ass1$phi[i,])) * sign(ass1$phi[i,]), "*"))
  }
  else {
      for(i in 1:dim(ass1$row)[3])
           ret <- c(ret, ass1$row[,,i], ass1$col[,,i],
                    sweep(ass1$row[,,i, drop=FALSE], 2, sqrt(abs(ass1$phi[i,])) * sign(ass1$phi[i,]), "*"),
                    sweep(ass1$col[,,i, drop=FALSE], 2, sqrt(abs(ass1$phi[i,])) * sign(ass1$phi[i,]), "*"))
  }

  # For double association models like some hmskew and yrcskew variants
  if(!is.null(assoc2)) {
      ass2 <- assoc2(model, weighting=weighting)
      ret <- c(ret, t(ass2$phi))

      if(dim(ass2$row)[3] == 1) {
          # Even if the scores are the same for all layers, we replicate them for simplicity's sake
          for(i in 1:nrow(ass2$phi))
               ret <- c(ret, ass2$row[,,1], ass2$col[,,1],
                        sweep(ass2$row[,,1, drop=FALSE], 2, sqrt(abs(ass2$phi[i,])) * sign(ass2$phi[i,]), "*"),
                        sweep(ass2$col[,,1, drop=FALSE], 2, sqrt(abs(ass2$phi[i,])) * sign(ass2$phi[i,]), "*"))
      }
      else {
          for(i in 1:dim(ass2$row)[3])
               ret <- c(ret, ass2$row[,,i], ass2$col[,,i],
                        sweep(ass2$row[,,i, drop=FALSE], 2, sqrt(abs(ass2$phi[i,])) * sign(ass2$phi[i,]), "*"),
                        sweep(ass2$col[,,i, drop=FALSE], 2, sqrt(abs(ass2$phi[i,])) * sign(ass2$phi[i,]), "*"))
      }

      # Without the quote(NULL), update.gnm() does call$start <- NULL, which removes it,
      # and eventually restores the default value (NA)
      suppressWarnings(model2 <- update(model, start=quote(NULL), etastart=NULL,
                                        verbose=TRUE, trace=TRUE, se="none"))

      # Random starting values can fail, and we still need to model to know the length of the result
      if(!is.null(model2))
          model <- model2

  if(!model$converged) {
      warning("Model failed to converge three times: ignoring the results of this replicate. Standard errors may not be completely accurate. Consider raising the value of iterMax.")
  }

  if(!model$converged)
      ret[] <- NA

  ret
}


se <- function(x, ...) UseMethod("se", x)

se.default <- function(x, ...) gnm::se(x)

se.rc <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "rc")) 
      stop("x must be a rc object")

  if(length(x$assoc) == 0)
      stop("x must have an association component")

  se.assoc(x$assoc)
}

se.rcL <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "rcL")) 
      stop("x must be a rcL object")

  if(length(x$assoc) == 0)
      stop("x must have an association component")

  se.assoc(x$assoc)
}

se.hmskew <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "hmskew"))
      stop("x must be a hmskew object")

  if(length(x[["assoc"]]) > 0 && length(x$assoc.hmskew) > 0)
      return(list(assoc=se.assoc(x$assoc), assoc.hmskew=se.assoc(x$assoc.hmskew)))
  if(length(x[["assoc"]]) > 0)
      return(se.assoc(x$assoc))
  else if(length(x$assoc.hmskew) > 0)
      return(se.assoc(x$assoc.hmskew))
  else
      stop("x must have an association or a skew-association component")
}

se.yrcskew <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "yrcskew"))
      stop("x must be a yrcskew object")

  if(length(x[["assoc"]]) > 0 && length(x$assoc.yrcskew) > 0)
      return(list(assoc=se.assoc(x$assoc), assoc.yrcskew=se.assoc(x$assoc.yrcskew)))
  if(length(x[["assoc"]]) > 0)
      return(se.assoc(x$assoc))
  else if(length(x$assoc.yrcskew) > 0)
      return(se.assoc(x$assoc.yrcskew))
  else
      stop("x must have an association or a skew-association component")
}

se.assoc <- function(x, type=c("se", "quasi.se"), ...) {
  type <- match.arg(type)

  if(!inherits(x, "assoc"))
      stop("x must be an assoc object")

  if(x$covtype == "none" || length(x$covmat) == 0)
      stop("No covariance matrix found: use the 'se' argument when fitting model")

  if(!(ncol(x$row) == ncol(x$col) &&
       ncol(x$phi) == ncol(x$row)))
      stop("Invalid component length")

  nd <- ncol(x$phi)
  nl <- nrow(x$phi)

  if(nrow(x$covmat) != ncol(x$covmat) ||
     nrow(x$covmat) != nl * nd + 2 * nl * nd * (nrow(x$row) + nrow(x$col)))
      stop("Covariance matrix dimensions do not match association structure")

  std.errs <- list()

  if(type == "quasi.se") {
      get.se <- function(covmat) qvcalc::qvcalc(covmat)$qvframe$quasiSE
  }
  else {
      get.se <- function(covmat) sqrt(diag(covmat))
  }

  covmat <- x$covmat

  std.errs$phi <- x$phi
  for(i in 1:nl)
      std.errs$phi[i,] <- get.se(covmat[seq((i - 1) * nd + 1, i * nd),
                                        seq((i - 1) * nd + 1, i * nd), drop=FALSE])

  std.errs$row <- x$row
  n <- nrow(std.errs$row)
  for(i in 1:dim(std.errs$row)[3]) {
      for(j in 1:nd) {
          int <- seq(nl * nd + n * (j - 1) + 1, nl * nd + n * j)
          std.errs$row[,j,i] <- get.se(covmat[int, int, drop=FALSE])
      }
  }

  std.errs$col <- x$col
  if(inherits(x, "assoc.symm")) {
      std.errs$col <- std.errs$row
  }
  else {
      n <- nrow(std.errs$col)
      start <- nl * nd + nd * nrow(std.errs$row)
      for(i in 1:dim(std.errs$col)[3]) {
          for(j in 1:nd) {
              int <- seq(start + n * (j - 1) + 1, start + n * j)
              std.errs$col[,j,i] <- get.se(covmat[int, int, drop=FALSE])
          }
      }
  }

  std.errs
}

# Vaguely adapted from bootstrap package, with three modifications:
# 1) Estimate the model only once for each cell of the table,
# and compute a mean weighted by cell frequencies at the end
# (Wong, Association models, 2010, p. 28-29)
# 2) theta must return a vector rather than a single number,
# to allow handling a series of parameters in the same run
# 3) The variance-covariance matrix is returned, instead of standard errors
#
# Formula taken from original function and from e.g. Steinhauer & Wuergler, 2010,
# "Leverage and Covariance Matrix Estimation in Finite-Sample IV Regressions", p. 11.
# (adapted to work with table cell weighting).
#
# Objects needed by theta() on the cluster nodes can be passed as arguments and handled
# in theta(), or simply passed as *named* arguments to jackknife(): they will be exported
# the the nodes' environments. 
jackknife <- function(x, theta, ..., w=rep(1, length(x)),
                      ncpus=if(require(parallel)) getOption("cl.cores", detectCores()) else getOption("cl.cores", 2))
{
    call <- match.call()
    stopifnot(length(w) == length(w))
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
    tot <- sum(w)
    mean.u <- rowSums(sweep(u, 2, w, "*"))/tot
    jack.bias <- (tot - 1) * (mean.u - thetahat)
    dev.u <- sweep(u, 1, mean.u, "-")
    jack.vcov <- (tot - 1)/tot * sweep(dev.u, 2, w, "*") %*% t(dev.u)
    return(list(jack.vcov = jack.vcov, jack.bias = jack.bias, jack.values = u, 
                call = call))
}

# Additional arguments are needed so that update() finds them even when using parLapply
theta.assoc <- function(x, model, assoc1, assoc2, family, weighting, ..., base=NULL, verbose=FALSE) {
  data <- model$data

  if(verbose) {
      iter <- which(!1:length(data) %in% x)
      if(length(iter) == 1)
          cat("Iteration for cell", iter, "\n")
      else
          cat("Ending iteration\n")
  }

  if(sum(data[-x]) > 0) {
      tab <- matrix(-1, nrow(data), ncol(data))
      tab[x] <- 0

      model <- update(model, data=data+tab, start=coef(model), verbose=FALSE, trace=TRUE)

      if(!model$converged && !is.null(base)) {
          cat("Model for cell ", which(!1:length(data) %in% x),
              " did not converge, starting again with random values...\n")
          # If we don't specify start, old values are used, which can give very bad initial fits
          base <- update(base, data=data, start=rep(NA, length(coef(base))))
          coefs <- coef(base)
          coefs[base$constrain] <- base$constrainTo
          model <- update(model, iterMax=5 * model$iterMax,
                          start=c(coefs, rep(NA, length(coef(model)) - length(coefs))),
                          verbose=TRUE, trace=TRUE)
      }

      if(!model$converged)
          stop("Model for cell ", which(!1:length(data) %in% x), " did not converge!")
  }

  ass1 <- assoc1(model, weighting=weighting)
  ret <- c(ass1$phi, ass1$row, ass1$col,
           sweep(ass1$row, 2, sqrt(ass1$phi), "*"),
           sweep(ass1$col, 2, sqrt(ass1$phi), "*"))

  # For double association models like some hmskew and yrcskew variants
  if(!is.null(assoc2)) {
      ass2 <- assoc2(model, weighting=weighting)
      ret <- c(ret, ass2$phi, ass2$row, ass2$col,
                    sweep(ass2$row, 2, sqrt(ass2$phi), "*"),
                    sweep(ass2$col, 2, sqrt(ass2$phi), "*"))
  }

  ret
}

# theta.yrcskew <- function(x, model, assoc1, assoc2, family, weighting, ..., base=NULL, verbose=FALSE) {
#   require(gnm)
# 
#   data <- model$data
# 
#   if(verbose) {
#       iter <- which(!1:nrow(data) %in% x)
#       if(length(iter) == 1)
#           cat("Iteration for cell", iter, "\n")
#       else
#           cat("Ending iteration\n")
#   }
# 
#   if(sum(data[-x,"Freq"]) > 0) {
#       data[-x,"Freq"] <- data[-x,"Freq"] - 1
# 
#       model <- update(model, data=data, start=coef(model), verbose=FALSE, trace=TRUE)
# 
#       if(!model$converged && !is.null(base)) {
#           cat("Model for cell ", which(!1:nrow(data) %in% x),
#               " did not converge, starting again with random values...\n")
#           # If we don't specify start, old values are used, which can give very bad initial fits
#           base <- update(base, data=data, start=rep(NA, length(coef(base))))
#           coefs <- coef(base)
#           coefs[base$constrain] <- base$constrainTo
#           model <- update(model, iterMax=5 * model$iterMax,
#                           start=c(coefs, rep(NA, length(coef(model)) - length(coefs))),
#                           verbose=TRUE, trace=TRUE)
#       }
# 
#       if(!model$converged)
#           stop("Model for cell ", which(!1:nrow(data) %in% x), " did not converge!")
#   }
# 
#   ass1 <- assoc1(model, weighting=weighting)
#   ret <- c(ass1$phi, ass1$row, ass1$col,
#            sweep(ass1$row, 2, sqrt(ass1$phi), "*"),
#            sweep(ass1$col, 2, sqrt(ass1$phi), "*"))
# 
#   # For double association models like some hmskew and yrcskew variants
#   if(!is.null(assoc2)) {
#       ass2 <- assoc2(model, weighting=weighting)
#       ret <- c(ret, ass2$phi, ass2$row, ass2$col,
#                     sweep(ass2$row, 2, sqrt(ass2$phi), "*"),
#                     sweep(ass2$col, 2, sqrt(ass2$phi), "*"))
#   }
# 
#   ret
# }

se.rc <- function(model, type=c("se", "quasi.se")) {
  if(!inherits(model, "rc")) 
      stop("model must be a rc object")

  if(length(model$assoc) == 0)
      stop("model must have an association component")

  se.assoc(model$assoc)
}

se.hmskew <- function(model, type=c("se", "quasi.se")) {
  if(!inherits(model, "hmskew")) 
      stop("model must be a hmskew object")

  if(length(model$assoc) > 0 && length(model$assoc.hmskew) > 0)
      return(list(assoc=se.assoc(model$assoc), assoc.hmskew=se.assoc(model$assoc.hmskew)))
  if(length(model$assoc) > 0)
      return(se.assoc(model$assoc))
  else if(length(model$assoc.hmskew) > 0)
      return(se.assoc(model$assoc.hmskew))
  else
      stop("model must have an association or a skew-association component")
}

se.assoc <- function(ass, type=c("se", "quasi.se")) {
  type <- match.arg(type)

  if(!inherits(ass, "assoc"))
      stop("ass must be an assoc object")

  if(ass$covtype == "none" || length(ass$covmat) == 0)
      stop("No covariance matrix found: use the 'std.err' argument when fitting model")

  if(!(ncol(ass$row) == ncol(ass$col) &&
       length(ass$phi) == ncol(ass$row)))
      stop("Invalid component length")

  nd <- ncol(ass$row)

  if(nrow(ass$covmat) != ncol(ass$covmat) ||
     nrow(ass$covmat) != nd + 2 * nd * nrow(ass$row) + 2 * nd * nrow(ass$col))
      stop("Covariance matrix dimensions do not match association structure")

  std.errs <- list()

  if(type == "quasi.se") {
      get.se <- function(covmat) qvcalc::qvcalc(covmat)$qvframe$quasiSE
  }
  else {
      get.se <- function(covmat) sqrt(diag(covmat))
  }

  covmat <- ass$covmat

  std.errs$phise <- ass$phi
  std.errs$phise[] <- get.se(covmat[1:nd, 1:nd, drop=FALSE])

  std.errs$row <- ass$row
  n <- nrow(std.errs$row)
  for(i in 1:nd) {
      int <- seq(nd + n * (i-1) + 1, nd + n * i)
      std.errs$row[,i] <- get.se(covmat[int, int, drop=FALSE])
  }

  std.errs$col <- ass$col
  if(inherits(ass, "assoc.rc.homog")) {
      std.errs$col <- std.errs$row
  }
  else {
      n <- nrow(std.errs$col)
      start <- nd + nd * nrow(std.errs$row)
      for(i in 1:nd) {
          int <- seq(start + n * (i-1) + 1, start + n * i)
          std.errs$col[,i] <- get.se(covmat[int, int, drop=FALSE])
      }
  }

  std.errs
}

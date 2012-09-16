## Skew-symmetric association model (van der Heijden & Mooijaart, 1995)

HMSkew <- function(row, col, inst=NULL) {
  list(predictors = list(substitute(row), substitute(col), substitute(row), substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("%s * %s - %s * %s",
                   predLabels[3], predLabels[2], predLabels[1], predLabels[4])
       },
       call = as.expression(match.call()),
       common = c(1, 1, 2, 2)
       )
   }
class(HMSkew) <- "nonlin"

hmskew <- function(tab, nd.symm=NA, diagonal=FALSE,
                   weighting=c("marginal", "uniform", "none"), se=c("none", "jackknife", "bootstrap"),
                   nreplicates=50, ncpus=getOption("boot.ncpus", if(require(parallel)) min(parallel::detectCores(), 4) else 1),
                   family=poisson, weights=NULL, start=NA, etastart=NULL, tolerance=1e-6, iterMax=15000,
                   trace=TRUE, verbose=TRUE, ...) {
  weighting <- match.arg(weighting)
  se <- match.arg(se)
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(nrow(tab) != ncol(tab))
      stop("tab must be a square table for asymmetry models")

  if(!all(rownames(tab) == colnames(tab)))
      stop("tab must have identical row and column names for asymmetry models")

  if(!is.na(nd.symm) && nd.symm <= 0)
      stop("nd.symm must be NA or strictly positive")

  if(!is.na(nd.symm) && nd.symm/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric association cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(length(dim(tab)) > 2)
      tab <- margin.table(tab, 1:2)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  if(diagonal && !is.na(nd.symm))
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""

  if(is.na(nd.symm))
      basef <- sprintf("Freq ~ %s + %s %s+ Symm(%s, %s)",
                       vars[1], vars[2], diagstr, vars[1], vars[2])
  else
      basef <- sprintf("Freq ~ %s + %s %s+ instances(MultHomog(%s, %s), %i)",
                       vars[1], vars[2], diagstr, vars[1], vars[2], nd.symm)


  if(!is.null(start) && is.na(start)) {
      # Without good starting values, estimation can fail when running start-up iterations
      # with large tables
      cat("Running base model to find starting values...\n")

      # Setting tolerance to a value below 1e-6 can lead to convergence issues with large tables
      # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
      args <- list(formula=as.formula(basef),
                   data=tab, family=family,
                   tolerance=1e-3, iterMax=iterMax, trace=trace)
      base <- do.call("gnm", c(args, list(...)))

      start <- c(parameters(base), rep(NA, 2 * nrow(tab)))

      if(is.null(etastart))
          etastart <- as.numeric(predict(base))

      cat("Running real model...\n")
  }

  f <- sprintf("%s + HMSkew(%s, %s)", basef, vars[1], vars[2])

  # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
  args <- list(formula=as.formula(f), data=tab,
               family=family, start=start, etastart=etastart,
               tolerance=tolerance, iterMax=iterMax, trace=trace)

  model <- do.call("gnm", c(args, list(...)))

  if(is.null(model))
      return(NULL)

  if(!is.na(nd.symm)) {
      model$assoc <- assoc.rc.symm(model, weighting=weighting)
      class(model) <- c("hmskew", "rc.symm", "rc", "assocmod", class(model))
  }
  else {
      model$assoc <- list()
      class(model) <- c("hmskew", "assocmod", class(model))
  }

  model$assoc.hmskew <- assoc.hmskew(model, weighting=weighting)

  model$call <- match.call()


  if(se %in% c("jackknife", "bootstrap")) {
      cat("Computing", se, "standard errors...\n")

      if(is.null(ncpus))
          ncpus <- if(require(parallel)) min(parallel::detectCores(), 5) else 1

      assoc1 <- if(is.na(nd.symm)) assoc.hmskew else assoc.rc.symm
      assoc2 <- if(is.na(nd.symm)) NULL else assoc.hmskew

      if(se == "jackknife") {
          covmat <- jackknife(1:length(tab), jackknife.assoc, w=tab, ncpus=ncpus,
                              model=model, assoc1=assoc1, assoc2=assoc2,
                              weighting=weighting, family=family, ...,
                              base=base, verbose=FALSE)$jack.vcov

          boot.results <- numeric(0)
      }
      else {
          if(!is.null(weights))
              boot.weights <- rep.int(weights, tab)
          else
              boot.weights <- NULL

          boot.results <- boot::boot(1:sum(tab), boot.assoc,
                                     R=nreplicates, ncpus=ncpus, parallel="snow", weights=boot.weights,
                                     args=list(model=model, assoc1=assoc1, assoc2=assoc2,
                                               weighting=weighting, family=family, ...,
                                               weights=weights, base=base))

          covmat <- cov(boot.results$t)
      }

      if(!is.na(nd.symm)) {
          lim <- nd.symm + nd.symm * nrow(tab) + nd.symm * ncol(tab)
          model$assoc$boot.results <- boot.results
          model$assoc$covmat <- covmat[1:lim, 1:lim]
          model$assoc$covtype <- se
      }
      else {
          lim <- 0
      }

      model$assoc.hmskew$covmat <- covmat[seq(lim + 1, nrow(covmat)), seq(lim + 1, ncol(covmat))]
      model$assoc.hmskew$covtype <- se
  }
  else {
      model$assoc$boot.results <- numeric(0)
      model$assoc$covmat <- numeric(0)
      model$assoc$covtype <- "none"
  }

  model
}

assoc.hmskew <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data))])

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, et Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weighting == "uniform")
      p <- rep(1/nrow(tab), nrow(tab))
  else
      p <- rep(1, nrow(tab))

  mu <- parameters(model)[pickCoef(model, sprintf("HMSkew.*(\\Q%s\\E)",
                                            paste(rownames(tab), collapse="\\E|\\Q")))]
  mu1 <- mu[1:nrow(tab)]
  mu2 <- mu[-(1:nrow(tab))]

  if(!(length(mu1) == nrow(tab) && length(mu2) == nrow(tab)))
      stop("No dimensions found. Are you sure this is a van der Heijden & Mooijaart skewness model?")


  sc <- cbind(mu1, mu2)

  if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], dim(tab)[3], nrow(tab))
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
  else
      dg <- numeric(0)


  ## Normalize, cf. Clogg & Shihadeh (1994), eq. 5.3 et 5.4 (p. 83)
  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")
  # Technique proposed in van der Heijden & Mooijaart (1995), Appendix C
  # Weighted SVD taken from Goodman (1991), Appendix 4
  lambda <- sc[,2] %o% sc[,1] - sc[,1] %o% sc[,2]
  lambda0 <- sqrt(p %o% p) * lambda # Eq. A.4.3
  sv <- svd(lambda0)
  sc[,1:2] <- diag(1/sqrt(p)) %*% sv$u[,1:2] # Eq. A.4.7
  phi <- sv$d[1:2]

  # Since both dimensions share the same singular value, we cannot distinguish them
  # and their order is random. Use the convention that we want a positive skew association for cell (1, 2)
  # (implicit in the original article)
  if(sc[1,2] * sc[2,1] - sc[1,1] * sc[2,2] < 0)
      sc <- sc[,2:1]

  # Since rotation is also random, align first category to position 0 on the vertical axis,
  # and on the positive side of the horizontal axis (like original article does with *last* category)
  if(abs(sc[1,2]) > .Machine$double.eps)
      angle <- acos(sc[1,1]/sqrt(sum(sc[1,]^2)))
  else if (sc[1,1] < 0)
      angle <- pi
  else
      angle <- 0

  if(sc[1,2] < -.Machine$double.eps)
      angle <- -angle

  sc <- sc %*% matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)

  ## Prepare objects
  phi <- rbind(c(phi))
  dim(sc)[3] <- 1
  colnames(sc) <- colnames(phi) <- paste("Dim", 1:2, sep="")
  rownames(sc) <- rownames(tab)

  if(length(dg) > 0) {
      # Diag() sorts coefficients alphabetically!
      dg[,order(rownames(tab))] <- dg

      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")

      if(nrow(dg) > 1)
          rownames(dg) <- dimnames(tab)[[3]]
      else
          rownames(dg) <- "All levels"
  }

  obj <- list(phi = phi, row = sc, col = sc, diagonal = dg,
              weighting = weighting, row.weights = p, col.weights = p)

  class(obj) <- c("assoc.hmskew", "assoc.symm", "assoc")
  obj
}

assoc <- function(model, weighting, ...) UseMethod("assoc", model)

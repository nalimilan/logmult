hmskewL <- function(tab, nd.symm=NA, layer.effect.skew=c("homogeneous.scores", "heterogeneous", "none"),
                    layer.effect.symm=c("uniform", "homogeneous.scores", "heterogeneous", "none"),
                    diagonal=c("none", "heterogeneous", "homogeneous"),
                    weighting=c("marginal", "uniform", "none"), se=c("none", "jackknife"),
                    family=poisson, start=NA, etastart=NULL, tolerance=1e-6, iterMax=5000,
                    trace=TRUE, verbose=TRUE, ...) {
  layer.effect.skew <- match.arg(layer.effect.skew)
  layer.effect.symm <- match.arg(layer.effect.symm)
  diagonal <- match.arg(diagonal)
  weighting <- match.arg(weighting)
  se <- match.arg(se)
  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

  if(!is.na(nd.symm) && nd.symm <= 0)
      stop("nd must be strictly positive")

  if(nrow(tab) != ncol(tab))
      stop("tab must be a square table for asymmetry models")

  if(!all(rownames(tab) == colnames(tab)))
      stop("tab must have identical row and column names for asymmetry models")

  if(!is.na(nd.symm) && nd.symm <= 0)
      stop("nd.symm must be NA or strictly positive")

  if(!is.na(nd.symm) && nd.symm/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric association cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(is.na(nd.symm) && diagonal != "none")
     stop("diagonal parameters are redundant with nd.symm=NA")

  if(is.na(nd.symm) && layer.effect.symm == "homogeneous.scores")
      stop("layer.effect.symm == \"homogeneous.scores\" is only valid when nd.symm is not NA")
  else if(!is.na(nd.symm) && layer.effect.symm == "uniform")
      stop("layer.effect.symm == \"uniform\" is only valid when nd.symm is NA")

  if(length(dim(tab)) > 3)
      tab <- margin.table(tab, 1:3)


  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  if(diagonal == "heterogeneous")
      diagstr <- sprintf("+ %s:Diag(%s, %s) ", vars[3], vars[1], vars[2])
  else if(diagonal == "homogeneous")
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""

  f1 <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s",
                vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3])
  # FIXME: we should be able to eliminate 3:Symm(1, 2) in other cases, but this triggers a "numerically singular system" error in the last model
  #if(is.na(nd.symm) && layer.effect.symm == "none")
  #   eliminate <- eval(parse(text=sprintf("quote(Symm(%s, %s))", vars[1], vars[2])))
  #else
      eliminate <- eval(parse(text=sprintf("quote(%s:%s)", vars[1], vars[3])))

  base <- NULL

  nastart <- length(start) == 1 && is.na(start)


  if(nastart) {
      cat("Running base model to find starting values...\n")

      # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
      args <- list(formula=as.formula(paste(f1, diagstr)), data=tab,
                   family=family, eliminate=eliminate,
                   tolerance=1e-6, iterMax=iterMax)
      dots <- as.list(substitute(list(...)))[-1]
      args <- c(args, dots)

      base <- do.call("gnm", args)

      start <- parameters(base)
  }

  if(is.na(nd.symm)) {
      if(layer.effect.symm == "homogeneous.scores") {
          # Handled at the top of the function
          stop()
      }
      else if(layer.effect.symm == "heterogeneous") {
          f2 <- sprintf("+ %s:Symm(%s, %s)",
                        vars[3], vars[1], vars[2])

          if(nastart)
              start <- c(parameters(base), rep(NA, dim(tab)[3] * ((nrow(tab)^2 + nrow(tab))/2 - 1)))
      }
      else if(layer.effect.symm == "uniform") {
          f2 <- sprintf("+ Mult(Exp(%s), Symm(%s, %s))", 
                        vars[3], vars[1], vars[2])

          if(nastart)
              start <- c(parameters(base), rep(NA, dim(tab)[3] + (nrow(tab)^2 + nrow(tab))/2))
      }
      else {
          f2 <- sprintf("+ Symm(%s, %s)", 
                        vars[1], vars[2])

          if(nastart)
              start <- c(parameters(base), rep(NA, (nrow(tab)^2 + nrow(tab))/2 - 1))
      }
  }
  else {
      if(layer.effect.symm == "uniform") {
          # Handled at the top of the function
          stop()
      }
      else if(layer.effect.symm == "homogeneous.scores") {
          f2 <- sprintf("+ MultHomog(%s, %s:%s))",
                        vars[3], vars[1], vars[2])

          if(nastart)
              start <- c(parameters(base), rep(NA, nd.symm * (nrow(tab) + dim(tab)[3] + 2 * dim(tab)[3] * nrow(tab))))
      }
      else if(layer.effect.symm == "heterogeneous") {
          f2 <- sprintf("+ instances(MultHomog(%s:%s, %s:%s), %i)",
                        vars[3], vars[1], vars[3], vars[2], nd.symm)

          if(nastart)
              start <- c(parameters(base), rep(NA, nd.symm * dim(tab)[3] * nrow(tab) + 2 * dim(tab)[3] * nrow(tab)))
      }
      else {
          f2 <- sprintf("+ instances(MultHomog(%s, %s), %i)",
                        vars[1], vars[2], nd.symm)

          if(nastart)
              start <- c(parameters(base), rep(NA, nd.symm * nrow(tab) + 2 * dim(tab)[3] * nrow(tab)))
      }
  }

  if(layer.effect.skew == "heterogeneous") {
      f2.skew <- sprintf("+ HMSkew(%s:%s, %s:%s)",
                         vars[1], vars[3], vars[2], vars[3])

      if(nastart)
          start <- c(start, rep(NA, dim(tab)[3] * 2 * nrow(tab)))
  }
  else if(layer.effect.skew == "homogeneous.scores") {
      f2.skew <- sprintf("+ Mult(Exp(%s), HMSkew(%s, %s))",
                         vars[3], vars[1], vars[2])

      if(nastart)
          start <- c(start, rep(NA, dim(tab)[3] + 2 * nrow(tab)))
  }
  else {
      f2.skew <- sprintf("+ HMSkew(%s, %s)", 
                         vars[1], vars[2])

      if(nastart)
          start <- c(start, rep(NA, 2 * nrow(tab)))
  }

  # Heterogeneous diagonal parameters can make the convergence really slow unless
  # correct starting values are used
  if(!is.null(base)) {
      cat("Running second base model to find starting values...\n")

      # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
      args <- list(formula=as.formula(paste(f1, diagstr, f2, f2.skew)),
                   data=tab, family=family, start=start,
                   eliminate=eliminate, constrain=seq(1, length(parameters(base))), constrainTo=parameters(base),
                   tolerance=1e-3, iterMax=iterMax, trace=trace)
      dots <- as.list(substitute(list(...)))[-1]
      args <- c(args, dots)

      base2 <- do.call("gnm", args)

      # If model fails (can always happen), do not fail completely but start with random values
      if(is.null(base2))
          start <- NULL
      else
          start <- parameters(base2)

      if(is.null(etastart))
          etastart <- as.numeric(predict(base2))
  }

  if(!is.null(base) && is.null(etastart))
      etastart <- as.numeric(predict(base))

  if(!is.null(base))
      cat("Running real model...\n")

  # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
  args <- list(formula=as.formula(paste(f1, diagstr, f2, f2.skew)), data=tab,
               family=family, start=start, etastart=etastart,
               eliminate=eliminate,
               tolerance=tolerance, iterMax=iterMax, trace=trace)
  dots <- as.list(substitute(list(...)))[-1]
  args <- c(args, dots)

  model <- do.call("gnm", args)

  if(is.null(model))
      return(NULL)

  class(model) <- c("hmskewL", class(model))

  model$call <- match.call()

  if(is.na(nd.symm))
      model$assoc <- NULL
  else if(layer.effect.symm == "none")
      model$assoc <- assoc.rc.symm(model, weighting=weighting)
  else
      model$assoc <- assoc.rcL.symm(model, weighting=weighting)

  if(!is.null(model$assoc))
      class(model$assoc) <- c("assoc.rcL", "assoc")


  if(layer.effect.skew == "none")
      model$assoc.hmskew <- assoc.hmskew(model, weighting=weighting)
  else
      model$assoc.hmskew <- assoc.hmskewL(model, weighting=weighting)

  class(model$assoc.hmskew) <-  c("assoc.hmskewL", "assoc")


  if(se == "jackknife") {
      cat("Computing jackknife standard errors...\n")
      assoc1 <- if(is.na(nd.symm)) assoc.hmskewL
                else if(!is.na(nd.symm) && layer.effect.symm == "none") assoc.rc.symm
                else assoc.rcL.symm
      assoc2 <- if(is.na(nd.symm)) NULL else assoc.hmskewL

      covmat <- jackknife(1:length(tab), w=tab, theta.assoc, model, assoc1, assoc2,
                          family, weighting, HMSkew=HMSkew,
                          base=if(!is.null(base2)) base2 else base, verbose=verbose)$jack.vcov


      if(!is.na(nd.symm)) {
          if(layer.effect.symm == "heterogeneous")
              lim <- nd.symm * dim(tab)[3] + nd.symm * nrow(tab) * dim(tab)[3] + nd.symm * ncol(tab) * dim(tab)[3]
          else if(layer.effect.symm == "homogeneous.scores")
              lim <- nd.symm * dim(tab)[3] + nd.symm * nrow(tab) + nd.symm * ncol(tab)
          else if(layer.effect.symm == "none")
              lim <- nd.symm + nd.symm * nrow(tab) + nd.symm * ncol(tab)
          else # "uniform", handled at the top of the function
              stop()

          model$assoc$covmat <- covmat[1:lim, 1:lim]
          model$assoc$covtype <- "jackknife"
      }
      else {
          lim <- 0
      }

      model$assoc.hmskew$covmat <- covmat[seq(lim + 1, nrow(covmat)), seq(lim + 1, ncol(covmat))]
      model$assoc.hmskew$covtype <- "jackknife"
  }
  else {
      if(!is.na(nd.symm)) {
          model$assoc$covmat <- numeric(0)
          model$assoc$covtype <- "none"
      }

      model$assoc.hmskew$covmat <- numeric(0)
      model$assoc.hmskew$covtype <- "none"
  }

  model
}


assoc.hmskewL <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data)),
                             !is.na(dimnames(model$data)[3])])

  nr <- nrow(tab)
  nc <- ncol(tab)
  nl <- dim(tab)[3]

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weighting == "uniform")
      p <- rep(1/nr, nr)
  else
      p <- rep(1, nr)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  homogeneous <- TRUE

  mu1 <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*\\)HMSkew.*\\)[.:]\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                   vars[3], vars[1], vars[2],
                                                   paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]
  mu2 <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*\\)HMSkew.*\\)[.:]\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)\\.1$",
                                                   vars[3], vars[1], vars[2],
                                                   paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]
  phi <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*HMSkew\\(.*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                             vars[3],
                                             paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]


  if(length(phi) == nl &&
     length(mu1) == length(mu2) &&
     length(mu1) == nr) {
      homogeneous <- TRUE
      sc <- array(NA, dim=c(nr, 2, 1))
      sc[,1,1] <- mu1
      sc[,2,1] <- mu2
  }
  else if(length(phi) == nl &&
          length(mu1) == length(mu2) &&
          length(mu1) == nr * nl) {
          homogeneous <- FALSE
          sc <- array(NA, dim=c(nr, nl, 1))
          sc[,1,] <- mu1
          sc[,2,] <- mu2
  }
  else {
      stop("No dimensions found. Are you sure this is a van der Heijden & Mooijaart skewness model with layer effect?")
  }

  layer <- matrix(phi, nl, 1)


  if(length(pickCoef(model, "Diag\\(")) > nr)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], nl, nr)
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nr)
  else
      dg <- numeric(0)



  ## Normalize, cf. Clogg & Shihadeh (1994), eq. 5.3 et 5.4 (p. 83)
  # Center
  sc <- sweep(sc, 2:3, margin.table(sweep(sc, 1, p/sum(p), "*"), 2:3), "-")

  for(l in 1:dim(sc)[3]) {
      # Technique proposed in van der Heijden & Mooijaart (1995), Appendix C
      # Weighted SVD taken from Goodman (1991), Appendix 4
      lambda <- sc[,2,l] %o% sc[,1,l] - sc[,1,l] %o% sc[,2,l]
      lambda0 <- sqrt(p %o% p) * lambda # Eq. A.4.3
      sv <- svd(lambda0)
      sc[,1:2,l] <- diag(1/sqrt(p)) %*% sv$u[,1:2] # Eq. A.4.7
      phi <- sv$d[1]

      # Integrate phi to layer coefficients
      layer[,l] <- sweep(layer[,l, drop=FALSE], 2, phi, "*")
  }

  for(l in 1:dim(sc)[3]) {
      # Since both dimensions share the same singular value, we cannot distinguish them
      # and their order is random. Use the convention that we want a positive skew association for cell (1, 2)
      # (implicit in the original article)
      if(sc[1,2,l] * sc[2,1,l] - sc[1,1,l] * sc[2,2,l] < 0)
          sc[,,l] <- sc[,2:1,l]

      # Since rotation is also random, align first category to position 0 on the vertical axis,
      # and on the positive side of the horizontal axis (like original article does with *last* category)
      if(abs(sc[1,2,l]) > .Machine$double.eps)
          angle <- acos(sc[1,1,l]/sqrt(sum(sc[1,,l]^2)))
      else if (sc[1,1,l] < 0)
          angle <- pi
      else
          angle <- 0

      if(sc[1,2,l] < -.Machine$double.eps)
          angle <- -angle

      sc[,,l] <- sc[,,l] %*% matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)
  }

  ## Prepare objects
  layer <- cbind(layer, layer)

  rownames(sc) <- rownames(tab)
  colnames(sc) <- colnames(layer) <- paste("Dim", 1:2, sep="")
  rownames(layer) <- dimnames(tab)[[3]]

  if(!homogeneous)
      dimnames(sc)[[3]] <- dimnames(tab)[[3]]
  else
      dimnames(sc)[[3]] <- "All levels"

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

  obj <- list(phi = layer, row = sc, col = sc, diagonal = dg,
              weighting = weighting, row.weights = p, col.weights = p)

  class(obj) <- c("assoc.hmskewL", "assoc")
  obj
}

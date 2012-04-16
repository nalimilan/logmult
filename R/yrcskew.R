## RC(M) models with skew-symmetric association (Yamaguchi, 1990)

YRCSkew <- function(mincat, maxcat, skew, inst=NULL) {
  list(predictors = list(substitute(mincat), substitute(maxcat), skew=substitute(skew)),
       term = function(predLabels, varLabels) {
           sprintf("%s * %s * (%s - %s)",
                   predLabels[3], predLabels[1], predLabels[2], predLabels[1])
       },
       call = as.expression(match.call()),
       common = c(1, 1, 2)
       )
  }
class(YRCSkew) <- "nonlin"

yrcskew <- function(tab, nd.symm=NA, nd.skew=1, diagonal=FALSE,
                    weights=c("marginal", "uniform", "none"), std.err=c("none", "jackknife"),
                    family=poisson, start=NULL, tolerance=1e-12, iterMax=15000, trace=TRUE, ...) {
  weights <- match.arg(weights)
  std.err <- match.arg(std.err)
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(nrow(tab) != ncol(tab) || !all(rownames(tab) == colnames(tab)))
      stop("tab must be a square table for asymmetry models")

  if(!is.na(nd.symm) && nd.symm <= 0)
      stop("nd.symm must be NA or strictly positive")

  if(is.na(nd.skew))
      stop("nd.skew must be strictly positive")

  if(nd.skew > 1)
      warning("nd.skew > 1 is untested. You are on your own!")

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

  # Shut up an R CMD check NOTE
  base <- NULL

  if(is.null(start)) {
      # Without good starting values, estimation can fail when running start-up iterations
      # with large tables
      cat("Running base model to find starting values...\n")

      eval(parse(text=sprintf("base <- gnm(%s, data=tab, family=family, iterMax=%i, trace=%s, ...)",
                              basef, iterMax, if(trace) "TRUE" else "FALSE")))

      coefs <- coef(base)
      coefs[base$constrain] <- base$constrainTo
      start <- c(coefs, rep(NA, nd.skew * (nrow(tab) + 1)))

      cat("Running real model...\n")
  }

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  # We integrate computations in the formula rather than doing them separately
  # because gnm() does not seem to allow using objects outside of the data argument
  # in formulas called from functions. This seems to be a problem with how the formula's
  # environment is handled, and also happens when calling gnm() directly without eval().
  f <- sprintf("%s + instances(YRCSkew(%s[ifelse(ordered(%s) < ordered(%s), %s, %s)], %s[ifelse(ordered(%s) > ordered(%s), %s, %s)], sign(as.numeric(%s) - as.numeric(%s))), %s)",
                                       basef,
                                       vars[1], vars[1], vars[2], vars[1], vars[2],
                                       vars[1], vars[1], vars[2], vars[1], vars[2],
                                       vars[2], vars[1], nd.skew)
  eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, start=start, tolerance=%e, iterMax=%i, trace=%s, ...)",
                           f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))

  if(is.null(model))
      return(NULL)

  class(model) <- c("yrcskew", "rc.symm", "rc", class(model))

  if(!is.na(nd.symm))
      model$assoc <- assoc.rc.symm(model, weights=weights)
  else
      model$assoc <- list()

  model$assoc.yrcskew <- assoc.yrcskew(model, weights=weights)


  if(std.err == "jackknife") {
      assoc1 <- if(is.na(nd.symm)) assoc.yrcskew else assoc.rc.symm
      assoc2 <- if(is.na(nd.symm)) NULL else assoc.yrcskew

      cat("Computing jackknife standard errors...\n")
      covmat <- jackknife(1:length(tab), w=tab, theta.assoc, model, assoc1, assoc2,
                          family, weights, YRCSkew=YRCSkew, base=base)$jack.vcov

      if(!is.na(nd.symm)) {
          lim <- nd.symm + nd.symm * nrow(tab) + nd.symm * ncol(tab)
          model$assoc$covmat <- covmat[1:lim, 1:lim]
          model$assoc$covtype <- "jackknife"
      }
      else {
          lim <- 0
      }

      model$assoc.yrcskew$covmat <- covmat[seq(lim + 1, nrow(covmat)), seq(lim + 1, ncol(covmat))]
      model$assoc.yrcskew$covtype <- "jackknife"
  }
  else if(length(model$assoc) > 0) {
      model$assoc$covmat <- numeric(0)
      model$assoc$covtype <- "none"
  }

  model
}

## Model with separate skew-symmetric association (RC_SK)
assoc.yrcskew <- function(model, weights=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

#   # Reproduce original 2x2 table if it was modified (e.g. by yrcskew())
#   # gnm doesn't include coefficients for NA row/columns, so get rid of them too
#   tab <- as.data.frame(model$data[!is.na(rownames(model$data)),
#                                   !colnames(model$data) %in% c("mincat", "maxcat", "skew")])
#   tab <- xtabs(sprintf("Freq ~ %s + %s", colnames(tab)[1], colnames(tab)[2]),
#                data=tab)
#   tab <- as.table(tab)
  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data))])

  weights <- match.arg(weights)

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, et Becker & Clogg (1989), p. 144.
  if(weights == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weights == "uniform")
      p <- rep(1/nrow(tab), nrow(tab))
  else
      p <- rep(1, nrow(tab))

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")


  ## Get skew association coefficients
  sc <- matrix(NA, nrow(tab), 0)

  nd <- 0
  while(TRUE) {
      musk <- coef(model)[pickCoef(model, sprintf("YRCSkew.*inst = %s\\)(%s)", nd+1, paste(rownames(tab), collapse="|")))]

      if(length(musk) != nrow(tab))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      sc <- cbind(sc, musk)
  }

  if(nd <= 0) {
      musk <- coef(model)[pickCoef(model, sprintf("YRCSkew.*\\)\\.%s.*[^\\.].$", vars[1]))]

      if(length(musk) == nrow(tab)) {
          nd <- 1
          sc <- cbind(sc, musk)
      }
      else {
          stop("No skew dimensions found. Are you sure this is a Yamaguchi RC_SK model?")
      }
  }

  skew <- coef(model)[pickCoef(model, "YRCSkew.*\\)skew")]
  if(length(skew) != nd)
      stop("Skew coefficients not found. Are you sure this is a Yamaguchi RC_SK model?")

  dg <- coef(model)[pickCoef(model, "Diag\\(")]
  dg <- dg[match(order(names(dg)), order(rownames(tab)))]

  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd) {
    lambda <- lambda + abs(skew[i]) * sc[,i] %o% sc[,i]
  }
  lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
  sv <- svd(lambda0)
  sc[] <- diag(1/sqrt(p)) %*% sv$u[,1:nd] # Eq. A.4.7
  # Preserve the sign of skew
  phisk <- sign(skew) * sv$d[1:nd]


  ## Prepare objects
  phi <- rbind(c(phisk))
  dim(sc)[3] <- 1
  names(phisk) <- NULL
  colnames(sc) <- colnames(phisk) <- paste("Dim", 1:nd, sep="")
  rownames(sc) <- rownames(tab)
  if(length(dg) > 0) {
      dg <- rbind(c(dg))
      colnames(dg) <- rownames(tab)
  }

  obj <- list(phi = phisk, row = sc, col = sc,  diagonal = dg,
              weighting = weights, row.weights = p, col.weights = p)

  class(obj) <- c("assoc.yrcskew", "assoc.symm", "assoc")
  obj
}


## Model with skew-symmetric association homogeneous to symmetric association (HM_(S+SK))
assoc.yrcskew.homog <- function(model, weights=c("marginal", "unit"), ...) {
  if(!"gnm" %in% class(model)) stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)), !is.na(colnames(model$data))])

  weights <- match.arg(weights)

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, et Becker & Clogg (1989), p. 144.
  if(weights == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else
      p <- rep(1/nrow(tab), nrow(tab))

  ## Get symmetric association coefficients
  sc <- matrix(NA, nrow(tab), 0)

  nd <- 0
  while(TRUE) {
      mu <- coef(model)[pickCoef(model, sprintf("MultHomog\\(.*inst = %s\\)", nd+1))]

      if(!(length(mu) == nrow(tab)))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      sc <- cbind(sc, mu)
  }

  if(nd <= 0) {
      mu <- coef(model)[pickCoef(model, "MultHomog\\(")]

      if(length(mu) == nrow(tab)) {
          nd <- 1
          sc <- cbind(sc, mu)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with symmetric row and column scores?")
      }
  }

  ## Get skew association coefficients
  scsk <- matrix(NA, nrow(tab), 0)

  ndsk <- 0
  while(TRUE) {
      musk <- coef(model)[pickCoef(model, sprintf("MultHomog.*skew.*inst = %s\\)\\)(%s)", nd+1, paste(rownames(tab), collapse="|")))]

      if(!(length(musk) == nrow(tab)))
          break

      # This is a correct dimension, add it
      ndsk <- ndsk + 1

      scsk <- cbind(scsk, musk)
  }

  if(ndsk <= 0) {
      musk <- coef(model)[pickCoef(model, sprintf("MultHomog.*skew.*\\)(%s)$", paste(rownames(tab), collapse="|")))]

      if(length(musk) == nrow(tab)) {
          ndsk <- 1
          scsk <- cbind(scsk, musk)
      }
      else {
          stop("No skew dimensions found. Are you sure this is a row-column association model with symmetric row and column scores plus skewness?")
      }
  }

  skew <- coef(model)[pickCoef(model, "MultHomog.*skew\\)$")]
  if(length(skew) != ndsk)
      stop("skew coefficients not found. Are you sure this is a row-column association model with symmetric row and column scores plus skewness?")

  dg <- coef(model)[pickCoef(model, "Diag\\(")]
  dg <- dg[match(order(names(dg)), order(rownames(tab)))]


  ## Normalize, cf. Clogg & Shihadeh (1994), eq. 5.3 et 5.4 (p. 83)
  ## Symmetric part
  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4
  # We use eigen() here rather than svd() because matrix is square and symmetric
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd)
      lambda <- lambda + sc[,i] %o% sc[,i]
  lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
  eigen <- eigen(lambda0, symmetric=TRUE)
  sc[,1:nd] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
  phi <- eigen$values[1:nd]

  ## Skew part
  # Center
  scsk <- sweep(scsk, 2, colSums(sweep(scsk, 1, p, "*")), "-")

  # Integrate skew coefficient, changing sign if needed
  if(skew < 0) {
      skew <- -skew
      scsk <- -scsk
  }
  scsk <- sweep(scsk, 2, sqrt(skew), "*")

  # Technique proposed in Goodman (1991), Appendix 4
  # We use eigen() here rather than svd() because matrix is square and symmetric
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:ndsk)
      lambda <- lambda + scsk[,i] %o% scsk[,i]
  lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
  eigen <- eigen(lambda0, symmetric=TRUE)
  scsk[,1:ndsk] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:ndsk] # Eq. A.4.7
  phisk <- eigen$values[1:ndsk]


  ## Prepare objects
  phi <- rbind(c(phi))
  dim(sc)[3] <- dim(scsk)[3] <- 1
  colnames(sc) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
  colnames(scsk) <- paste("Dim", 1:ndsk, sep="")
  rownames(sc) <- rownames(scsk) <- rownames(tab)
  if(length(dg) > 0) {
      dg <- rbind(dg)
      colnames(dg) <- rownames(tab)
  }

  obj <- list(phi = phi, row = sc, col = sc, phisk = phi, rowsk = scsk, colsk = scsk,
              diagonal = dg, weighting = weights, row.weights = p, col.weights = p)

  class(obj) <- c("assoc.yrcskew.homog", "assoc.symm", "assoc")
  obj
}

assoc <- function(model, weights, ...) UseMethod("assoc", model)

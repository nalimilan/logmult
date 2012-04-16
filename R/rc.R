## RC(M) model

rc <- function(tab, nd=1, homogeneous=FALSE, diagonal=FALSE,
               weights=c("marginal", "uniform", "none"), std.err=c("none", "jackknife"),
               family=poisson, start=NULL, tolerance=1e-12, iterMax=5000, trace=TRUE, ...) {
  weights <- match.arg(weights)
  std.err <- match.arg(std.err)
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(is.na(nd) || nd <= 0)
      stop("nd must be strictly positive")

  if(homogeneous && nd/2 > min(nrow(tab), ncol(tab)) - 1)
     stop("Number of dimensions of homogeneous model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!homogeneous && nd > min(nrow(tab), ncol(tab)) - 1)
     stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  if(diagonal)
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""

  if(homogeneous) {
      f <- sprintf("Freq ~ %s + %s %s+ instances(MultHomog(%s, %s), %i)",
                   vars[1], vars[2], diagstr, vars[1], vars[2], nd)

      if(is.null(start))
          eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, tolerance=%e, iterMax=%i, trace=%s, ...)",
                                  f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))
      else
          eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, start=start, tolerance=%e, iterMax=%i, trace=%s, ...)",
                                  f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))
  }
  else {
      if(is.null(start)) {
          base <- gnm(as.formula(sprintf("Freq ~ %s + %s %s", vars[1], vars[2], diagstr)),
                      family=family, data=tab)

          # residSVD evaluates the variable names in parent.frame(), which uses any object
          # called "vars" in the global environment if not handled like this
          res <- eval(parse(text=sprintf("residSVD(base, %s, %s, %i)", vars[1], vars[2], nd)))
          start <- c(rep(NA, length(coef(base))), res)
      }

      f <- sprintf("Freq ~ %s + %s %s+ instances(Mult(%s, %s), %i)",
                   vars[1], vars[2], diagstr, vars[1], vars[2], nd)
      eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, start=start, tolerance=%e, iterMax=%i, trace=%s, ...)",
                               f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))
  }

  if(is.null(model))
      return(NULL)

  newclasses <- if(homogeneous) c("rc.homog", "rc") else "rc"
  class(model) <- c(newclasses, class(model))

  model$assoc <- assoc(model, weights=weights)

  if(std.err == "jackknife") {
      cat("Computing jackknife standard errors...\n")
      model$assoc$covmat <- jackknife(1:length(tab), w=tab, theta.assoc, model,
                                      getS3method("assoc", class(model)), NULL,
                                      family, weights)$jack.vcov
      scnames <- c(t(outer(paste("D", 1:nd, " ", vars[1], ".", sep=""), rownames(tab), paste, sep="")),
                   t(outer(paste("D", 1:nd, " ", vars[2], ".", sep=""), colnames(tab), paste, sep="")))
      rownames(model$assoc$covmat) <- colnames(model$assoc$covmat) <-
          c(paste("Dim", 1:nd, sep=""), scnames, paste(scnames, "*", sep=""))

      model$assoc$covtype <- "jackknife"
  }
  else {
      model$assoc$covmat <- numeric(0)
      model$assoc$covtype <- "none"
  }

  model
}

assoc.rc <- function(model, weights=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data))])

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weights <- match.arg(weights)
  if(weights == "marginal") {
      rp <- prop.table(margin.table(tab, 1))
      cp <- prop.table(margin.table(tab, 2))
  }
  else if(weights == "uniform") {
      rp <- rep(1/nrow(tab), nrow(tab))
      cp <- rep(1/ncol(tab), ncol(tab))
  }
  else {
      rp <- rep(1, nrow(tab))
      cp <- rep(1, ncol(tab))
  }

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  # Prepare matrices before filling them
  row <- matrix(NA, nrow(tab), 0)
  col <- matrix(NA, ncol(tab), 0)

  nd <- 0
  while(TRUE) {
      mu <- coef(model)[pickCoef(model, sprintf("Mult.*inst = %s\\)\\.%s", nd+1, vars[1]))]
      nu <- coef(model)[pickCoef(model, sprintf("Mult.*inst = %s\\)\\.%s", nd+1, vars[2]))]

      if(!(length(mu) == nrow(tab) && length(nu) == ncol(tab)))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      row <- cbind(row, mu)
      col <- cbind(col, nu)
  }

  if(nd <= 0) {
      mu <- coef(model)[pickCoef(model, sprintf("Mult.*\\)\\.%s.*[^\\.].$", vars[1]))]
      nu <- coef(model)[pickCoef(model, sprintf("Mult.*\\)\\.%s.*[^\\.].$", vars[2]))]

      if(length(mu) == nrow(tab) && length(nu) == ncol(tab)) {
          nd <- 1

          row <- cbind(row, mu)
          col <- cbind(col, nu)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model?")
      }
  }

  dg <- coef(model)[pickCoef(model, "Diag\\(")]

  # Center
  row <- sweep(row, 2, colSums(sweep(row, 1, rp/sum(rp), "*")), "-")
  col <- sweep(col, 2, colSums(sweep(col, 1, cp/sum(cp), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd) {
      lambda <- lambda + row[,i] %o% col[,i]
  }
  lambda0 <- lambda * sqrt(rp %o% cp) # Eq. A.4.3
  sv <- svd(lambda0)
  row[] <- diag(1/sqrt(rp)) %*% sv$u[,1:nd] # Eq. A.4.7
  col[] <- diag(1/sqrt(cp)) %*% sv$v[,1:nd] # Eq. A.4.7
  phi <- sv$d[1:nd]

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first row category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(row[1,i] < 0) {
          row[,i] <- -row[,i]
          col[,i] <- -col[,i]
      }
  }

  ## Prepare objects
  phi <- rbind(c(phi))
  dim(row)[3] <- dim(col)[3] <- 1
  colnames(row) <- colnames(col) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
  rownames(row) <- rownames(tab)
  rownames(col) <- colnames(tab)

  if(length(dg) > 0) {
      dg <- rbind(c(dg))
      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab) else paste(rownames(tab), colnames(tab), sep=":")
  }

  obj <- list(phi = phi, row = row, col = col, diagonal = dg,
              weighting = weights, row.weights = rp, col.weights = cp)

  class(obj) <- c("rc.assoc", "assoc")
  obj
}

## RC(M) model with homogeneous row and column scores
assoc.rc.homog <- function(model, weights=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data))])

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weights <- match.arg(weights)
  if(weights == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weights == "uniform")
      p <- rep(1/nrow(tab), nrow(tab))
  else
      p <- rep(1, nrow(tab))


  sc <- matrix(NA, nrow(tab), 0)

  nd <- 0
  while(TRUE) {
      mu <- coef(model)[pickCoef(model, sprintf("MultHomog.*inst = %s\\)", nd+1))]

      if(!(length(mu) == nrow(tab)))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      sc <- cbind(sc, mu)
  }

  if(nd <= 0) {
      mu <- coef(model)[pickCoef(model, "MultHomog")]

      if(length(mu) == nrow(tab)) {
          nd <- 1
          sc <- cbind(sc, mu)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with homogeneous row and column scores?")
      }
  }

  dg <- coef(model)[pickCoef(model, "Diag\\(")]
  dg <- dg[match(order(names(dg)), order(rownames(tab)))]

  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4, but with eigenvalues decomposition
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd)
      lambda <- lambda + sc[,i] %o% sc[,i]
  lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
  eigen <- eigen(lambda0, symmetric=TRUE)
  sc[,1:nd] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
  phi <- eigen$values[1:nd]

  # FIXME: does not apply to homogeneous scores
  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(sc[1,i] < 0)
          sc[,i] <- -sc[,i]
  }

  ## Prepare objects
  phi <- rbind(c(phi))
  dim(sc)[3] <- 1
  colnames(sc) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
  rownames(sc) <- rownames(tab)
  if(length(dg) > 0) {
      dg <- rbind(c(dg))
      colnames(dg) <- rownames(tab)
  }


  obj <- list(phi = phi, row = sc, col= sc, diagonal = dg,
              weighting = weights, row.weights = p, col.weights = p)

  class(obj) <- c("rc.homog.assoc", "rc.assoc", "assoc")
  obj
}

assoc <- function(model, weights, ...) UseMethod("assoc", model)

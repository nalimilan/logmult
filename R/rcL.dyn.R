## RC(M) with regression-type and transitional variation

RCReg <- function(row, col, layer, inst=NULL) {
  list(predictors = list(R1=substitute(row), C1=substitute(col), substitute(layer), R2=substitute(row), C2=substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("(%s + %s * (%s)^2) * (%s + %s * (%s)^2)",
                   predLabels[1], predLabels[4], predLabels[3],
                   predLabels[2], predLabels[5], predLabels[3])
       },
       call = as.expression(match.call()),
       match = c(1, 2, 3, 1, 2),
       start = function(theta) {
           theta[attr(theta, "assign") == 3] <- sqrt(seq(0, 1, length.out=sum(attr(theta, "assign") == 3)))
           theta
       })
   }
class(RCReg) <- "nonlin"

RCRegHomog <- function(row, col, layer, inst=NULL) {
  list(predictors = list(R1=substitute(row), C1=substitute(col), substitute(layer), R2=substitute(row), C2=substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("(%s + %s * (%s)^2) * (%s + %s * (%s)^2)",
                   predLabels[1], predLabels[4], predLabels[3],
                   predLabels[2], predLabels[5], predLabels[3])
       },
       call = as.expression(match.call()),
       match = c(1, 2, 3, 1, 2),
       common = c(1, 1, 3, 2, 2),
       start = function(theta) {
           theta[attr(theta, "assign") == 3] <- sqrt(seq(0, 1, length.out=sum(attr(theta, "assign") == 3)))
           theta
       })
   }
class(RCRegHomog) <- "nonlin"

RCTrans <- function(row, col, layer, inst=NULL) {
  list(predictors = list(R1=substitute(row), C1=substitute(col), substitute(layer), R2=substitute(row), C2=substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("(%s * %s * (1 - (%s)^2)) + (%s * %s * (%s)^2)",
                   predLabels[1], predLabels[2], predLabels[3],
                   predLabels[4], predLabels[5], predLabels[3])
       },
       call = as.expression(match.call()),
       match = c(1, 2, 3, 1, 2),
       start = function(theta) {
           theta[attr(theta, "assign") == 3] <- seq(0, 1, length.out=sum(attr(theta, "assign") == 3))
           theta
       })
  }
class(RCTrans) <- "nonlin"

RCTransHomog <- function(row, col, layer, inst=NULL) {
  list(predictors = list(R1=substitute(row), C1=substitute(col), substitute(layer), R2=substitute(row), C2=substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("(%s * %s * (1 - (%s)^2)) + (%s * %s * (%s)^2)",
                   predLabels[1], predLabels[2], predLabels[3],
                   predLabels[4], predLabels[5], predLabels[3])
       },
       call = as.expression(match.call()),
       match = c(1, 2, 3, 1, 2),
       common = c(1, 1, 3, 2, 2),
       start = function(theta) {
           theta[attr(theta, "assign") == 3] <- seq(0, 1, length.out=sum(attr(theta, "assign") == 3))
           theta
       })
  }
class(RCTransHomog) <- "nonlin"

rcL.dyn <- function(tab, nd=1, type=c("regression", "transition"),
                    symmetric=FALSE, diagonal=c("none", "homogeneous", "heterogeneous"),
                    weighting=c("marginal", "uniform", "none"), std.err=c("none", "jackknife"),
                    family=poisson, start=NA, tolerance=1e-12, iterMax=50000, trace=TRUE, ...) {
  type <- match.arg(type)
  diagonal <- match.arg(diagonal)
  weighting <- match.arg(weighting)
  std.err <- match.arg(std.err)
  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

  if(is.na(nd) || nd <= 0)
      stop("nd must be strictly positive")

  if(symmetric && nd/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!symmetric && nd > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  if(length(dim(tab)) > 3)
      tab <- margin.table(tab, 1:3)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  if(diagonal == "homogeneous")
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else if(diagonal == "heterogeneous")
      diagstr <- sprintf("+ %s:Diag(%s, %s) ", vars[3], vars[1], vars[2])
  else
      diagstr <- ""

  variation <- if(type == "regression") "RCReg" else "RCTrans"

  if(symmetric)
      variation <- paste(variation, "Homog", sep="")


  if(!is.null(start) && is.na(start)) {
      cat("Running base model to find starting values...\n")

      if(symmetric) {
          base <- gnm(as.formula(sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s + instances(MultHomog(%s, %s), %i)",
                                         vars[1], vars[2], vars[3],
                                         vars[1], vars[3], vars[2], vars[3], diagstr,
                                         vars[1], vars[2], nd)),
                      family=family, data=tab)

          start <- coef(base)[seq(1, length(coef(base)) - nrow(tab) * nd)]
          for(i in 1:nd)
              start <- c(start, coef(base)[pickCoef(base, paste("Mult.*inst =", i))],
                         sqrt(seq(0, 1, length.out=dim(tab)[3])), rep(NA, nrow(tab)))
      }
      else {
          base <- gnm(as.formula(sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s + instances(Mult(%s, %s), %i)",
                                         vars[1], vars[2], vars[3],
                                         vars[1], vars[3], vars[2], vars[3], diagstr,
                                         vars[1], vars[2], nd)),
                      family=family, data=tab)

          start <- coef(base)[seq(1, length(coef(base)) - (nrow(tab) + ncol(tab)) * nd)]
          for(i in 1:nd)
              start <- c(start, coef(base)[pickCoef(base, paste("Mult.*inst =", i))],
                         sqrt(seq(0, 1, length.out=dim(tab)[3])), rep(NA, nrow(tab) + ncol(tab)))

      }

      cat("Running real model...\n")
    }

  f <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s+ instances(%s(%s, %s, %s), %i)",
               vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3], diagstr,
               variation, vars[1], vars[2], vars[3], nd)

  args <- list(formula=eval(as.formula(f)), data=substitute(tab),
               family=substitute(family),
               # For RCReg, both constraints are really needed: the computed scores are wrong without them
               # (rows are ordered along an oblique axis, and columns get weird values)
               constrain=sprintf("(RCReg|RCTrans).*\\).\\Q%s\\E(\\Q%s\\E|\\Q%s\\E)$",
                                 vars[3], head(dimnames(tab)[[3]], 1), tail(dimnames(tab)[[3]], 1), nd),
               constrainTo=rep(0:1, nd),
               start=start,
               tolerance=tolerance, iterMax=iterMax, trace=trace)
  dots <- as.list(substitute(list(...)))[-1]
  args <- c(args, dots)

  model <- do.call("gnm", args)


  if(is.null(model))
      return(NULL)

  newclasses <- if(symmetric) c("rcL.dyn.symm", "rcL.dyn", "rcL") else c("rcL.dyn", "rcL")
  class(model) <- c(newclasses, class(model))

  model$assoc <- assoc(model, weighting=weighting)

  if(std.err == "jackknife") {
      cat("Computing jackknife standard errors...\n")
      model$assoc$covmat <- jackknife(1:length(tab), w=tab, theta.assoc, model,
                                      getS3method("assoc", class(model)), NULL,
                                      family, weighting)$jack.vcov
      scnames <- t(outer(paste(vars[3], ".", dimnames(tab)[[3]], sep=""),
                         c(t(outer(paste("D", 1:nd, " ", vars[1], ".", sep=""), rownames(tab), paste, sep="")),
                           t(outer(paste("D", 1:nd, " ", vars[2], ".", sep=""), colnames(tab), paste, sep=""))),
                         paste))
      rownames(model$assoc$covmat) <- colnames(model$assoc$covmat) <-
          c(t(outer(paste(vars[3], dimnames(tab)[[3]], sep=""), paste("Dim", 1:nd, sep=""), paste)),
              scnames, paste(scnames, "*", sep=""))

      model$assoc$covtype <- "jackknife"
  }
  else {
      model$assoc$covmat <- numeric(0)
      model$assoc$covtype <- "none"
  }

  model
}

# Number of constraints applied:
#   - two layer coefficients (normally already present in the model)
#   - two for each dimension of the start scores
#   - two for each dimension of the end scores
assoc.rcL.dyn <- function(model, weighting=c("marginal", "uniform", "none")) {
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
  if(weighting == "marginal") {
      rp <- prop.table(margin.table(tab, 1))
      cp <- prop.table(margin.table(tab, 2))
  }
  else if(weighting == "uniform") {
      rp <- rep(1/nr, nr)
      cp <- rep(1/nc, nc)
  }
  else {
      rp <- rep(1, nr)
      cp <- rep(1, nc)
  }
  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")


  # Find out the number of dimensions
  nd <- 0
  while(length(pickCoef(model, sprintf("(RCReg|RCTrans).*inst = %s\\)\\.[RC][12]%s", nd+1, vars[1]))) > 0)
      nd <- nd + 1

  # One dimension, or none
  if(nd <= 0) {
      mu <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*\\)\\.R1\\Q%s\\E(\\Q%s\\E)$", vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*\\)\\.C1\\Q%s\\E(\\Q%s\\E)$", vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]

      mu1 <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*\\)\\.R2\\Q%s\\E(\\Q%s\\E)",vars[1],
                                                 paste(rownames(tab), collapse="\\E|\\Q")))]
      nu1 <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*\\)\\.C2\\Q%s\\E(\\Q%s\\E)", vars[2],
                                                 paste(colnames(tab), collapse="\\E|\\Q")))]

      phi <- coef(model)[pickCoef(model, sprintf("RCReg.*\\)\\.\\Q%s\\E(\\Q%s\\E)", vars[3],
                                                paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      if(length(mu) == nr && length(nu) == nc
         && length(mu1) == nr && length(nu1) == nc
         && length(phi) == nl) {
          nd <- 1

          row <- matrix(mu, nr, 1)
          col <- matrix(nu, nc, 1)
          row1 <- matrix(mu1, nr, 1)
          col1 <- matrix(nu1, nc, 1)
          layer <- matrix(phi, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with regression-type layer effect?")
      }
  }
  else {
      # Several dimensions: prepare matrices before filling them
      row <- matrix(NA, nr, nd)
      col <- matrix(NA, nc, nd)
      row1 <- matrix(NA, nr, nd)
      col1 <- matrix(NA, nc, nd)
      layer <- matrix(NA, nl, nd)

      for(i in 1:nd) {
          mu <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*inst = %s\\)\\.R1\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]
          nu <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*inst = %s\\)\\.C1\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[2],
                                                    paste(colnames(tab), collapse="\\E|\\Q")))]

          mu1 <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*inst = %s\\)\\.R2\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]
          nu1 <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*inst = %s\\)\\.C2\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[2],
                                                    paste(colnames(tab), collapse="\\E|\\Q")))]

          phi <- coef(model)[pickCoef(model, sprintf("(RCReg|RCTrans).*inst = %s\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[3],
                                                    paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

          if(length(mu) == nr && length(nu) == nc
               && length(mu1) == nr && length(nu1) == nc
               && length(phi) == nl) {
              row[,i] <- mu
              col[,i] <- nu
              row1[,i] <- mu1
              col1[,i] <- nu1
              layer[,i] <- phi
          }
          else {
              stop("Invalid dimensions found. Are you sure this is a row-column association model with regression-type layer effect?")
          }
      }
  }

  if(length(pickCoef(model, "Diag\\(") > 0)) {
      dg <- matrix(NA, nl, nr)
      dg[] <- coef(model)[pickCoef(model, "Diag\\(")]
  }
  else {
      dg <- numeric(0)
  }


  # Layer coefficients are squared internally by RCReg
  layer <- layer^2

  # Replace constrained coefficients
  layer[1,] <- 0
  layer[nrow(layer),] <- 1

  # Center
  row <- sweep(row, 2, colSums(sweep(row, 1, rp/sum(rp), "*")), "-")
  col <- sweep(col, 2, colSums(sweep(col, 1, cp/sum(cp), "*")), "-")
  row1 <- sweep(row1, 2, colSums(sweep(row1, 1, rp/sum(rp), "*")), "-")
  col1 <- sweep(col1, 2, colSums(sweep(col1, 1, cp/sum(cp), "*")), "-")


  ## Prepare objects
  if(any(grepl("RCTrans", names(coef(model))))) {
      # Scale
      phir <- sqrt(colSums(sweep(row^2, 1, rp, "*")))
      phic <- sqrt(colSums(sweep(col^2, 1, cp, "*")))
      row <- sweep(row, 2, phir, "/")
      col <- sweep(col, 2, phic, "/")
      phi <- phir * phic

      phi1r <- sqrt(colSums(sweep(row1^2, 1, rp, "*")))
      phi1c <- sqrt(colSums(sweep(col1^2, 1, cp, "*")))
      row1 <- sweep(row1, 2, phi1r, "/")
      col1 <- sweep(col1, 2, phi1c, "/")
      phi1 <- phi1r * phi1c

      layer <- cbind(sweep(1 - layer, 2, phi, "*"), sweep(layer, 2, phi1, "*"))
      colnames(layer) <- c(paste("Dim", 1:nd, " S", sep=""), paste("Dim", 1:nd, " E", sep=""))
      rownames(layer) <- dimnames(tab)[[3]]

      row <- array(cbind(row, row1), dim=c(nrow(row), 2 * nd, 1), dimnames=list(rownames(tab), colnames(layer)))
      col <- array(cbind(col, col1), dim=c(nrow(col), 2 * nd, 1), dimnames=list(colnames(tab), colnames(layer)))
  }
  else if(any(grepl("RCReg", names(coef(model))))) {
      row <- array(row, dim=c(nrow(row), nd, nl)) +
                 sweep(array(row1, dim=c(nrow(row1), nd, nl)), 3:2, layer, "*")
      col <- array(col, dim=c(nrow(col), nd, nl)) +
                 sweep(array(col1, dim=c(nrow(col1), nd, nl)), 3:2, layer, "*")

      phir <- sqrt(margin.table(sweep(row^2, 1, rp, "*"), 2:3))
      phic <- sqrt(margin.table(sweep(col^2, 1, cp, "*"), 2:3))
      row <- sweep(row, 2:3, phir, "/")
      col <- sweep(col, 2:3, phic, "/")
      layer <- t(phir * phic)

      # Order dimensions according to phi on first layer category
      ord <- order(abs(layer[1,]), decreasing=TRUE)
      layer <- layer[,ord, drop=FALSE]
      row <- row[,ord,, drop=FALSE]
      col <- col[,ord,, drop=FALSE]

      colnames(layer) <- colnames(row) <- colnames(col) <- paste("Dim", 1:nd, sep="")
      rownames(row) <- rownames(tab)
      rownames(col) <- colnames(tab)
      rownames(layer) <- dimnames(row)[[3]] <- dimnames(col)[[3]] <- dimnames(tab)[[3]]
  }
  else {
      stop("Invalid model")
  }


  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first row category: this ensures the results are stable when jackknifing.
  for(i in 1:dim(row)[3]) {
      for(j in 1:nd) {
          if(row[1,j,i] < 0) {
              row[,j,i] <- -row[,j,i]
              col[,j,i] <- -col[,j,i]
          }
      }
  }

  # Scores can switch sides from one layer to another
  # Reverse axes for each layer state so that the scores are positively correlated to the first layer scores for rows
  for(i in seq_len(dim(row)[3] - 1)) {
      chg <- ifelse(diag(as.matrix(cor(row[,,i+1], row[,,1]))) >= 0, 1, -1)
      row[,,i+1] <- sweep(row[,,i+1, drop=FALSE], 2, chg, "*")
      col[,,i+1] <- sweep(col[,,i+1, drop=FALSE], 2, chg, "*")
  }

  if(length(dg) > 0) {
      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")
      rownames(dg) <- dimnames(tab)[[3]]
  }

  obj <- list(phi = layer, row = row, col = col, diagonal = dg,
              weighting = weighting, row.weights = rp, col.weights = cp)

  class(obj) <- c("assoc.rcL.dyn", "assoc.rcL", "assoc")
  obj
}


assoc.rcL.dyn.symm <- function(model, weighting=c("marginal", "uniform", "none")) {
    if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data)),
                             !is.na(dimnames(model$data)[3])])

  nr <- nrow(tab)
  nc <- ncol(tab)
  nl <- dim(tab)[3]

  stopifnot(nr == nc)

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weighting == "uniform")
      rp <- rep(1/nr, nr)
  else
      rp <- rep(1, nr)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")


  # Find out the number of dimensions
  nd <- 0
  while(length(pickCoef(model, sprintf("(RCRegHomog|RCTransHomog).*inst = %s\\)[RC][12]\\.\\Q%s\\E\\|\\Q%s\\E",
                                       nd+1, vars[1], vars[2]))) > 0)
      nd <- nd + 1

  # One dimension, or none
  if(nd <= 0) {
      mu <- coef(model)[pickCoef(model, sprintf("(RCRegHomog|RCTransHomog).*\\)[RC]1\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]

      mu1 <- coef(model)[pickCoef(model, sprintf("(RCRegHomog|RCTransHomog).*\\)[RC]2\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                 vars[1], vars[2],
                                                 paste(rownames(tab), collapse="\\E|\\Q")))]

      phi <- coef(model)[pickCoef(model, sprintf("(RCRegHomog|RCTransHomog).*\\)\\.\\Q%s\\E(\\Q%s\\E)$", vars[3],
                                                paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      if(length(mu) == nr && length(mu1) == nr && length(phi) == nl) {
          nd <- 1

          sc <- matrix(mu, nr, 1)
          sc1 <- matrix(mu1, nr, 1)
          layer <- matrix(phi, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with regression-type layer effect?")
      }
  }
  else {
      # Several dimensions: prepare matrices before filling them
      sc <- matrix(NA, nr, nd)
      sc1 <- matrix(NA, nr, nd)
      layer <- matrix(NA, nl, nd)

      for(i in 1:nd) {
          mu <- coef(model)[pickCoef(model, sprintf("(RCRegHomog|RCTransHomog).*inst = %i\\)[RC]1\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1], vars[2],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]

          mu1 <- coef(model)[pickCoef(model, sprintf("(RCRegHomog|RCTransHomog).*inst = %i\\)[RC]2\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1], vars[2],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]

          phi <- coef(model)[pickCoef(model, sprintf("(RCRegHomog|RCTransHomog).*inst = %i\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[3],
                                                    paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

          if(length(mu) == nr && length(mu1) == nr && length(phi) == nl) {
              sc[,i] <- mu
              sc1[,i] <- mu1
              layer[,i] <- phi
          }
          else {
              stop("Invalid dimensions found. Are you sure this is a row-column association model with regression-type layer effect?")
          }
      }
  }

  if(length(pickCoef(model, "Diag\\(") > 0)) {
      dg <- matrix(NA, nl, nr)
      dg[] <- coef(model)[pickCoef(model, "Diag\\(")]
  }
  else {
      dg <- numeric(0)
  }


  # Layer coefficients are squared internally by RCReg
  layer <- layer^2

  # Replace constrained coefficients
  layer[1,] <- 0
  layer[nrow(layer),] <- 1

  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")
  sc1 <- sweep(sc1, 2, colSums(sweep(sc1, 1, p/sum(p), "*")), "-")


  ## Prepare objects
  if(any(grepl("RCTrans", names(coef(model))))) {
      # Scale
      phi <- sqrt(colSums(sweep(sc^2, 1, p, "*")))
      sc <- sweep(sc, 2, phi, "/")

      phi1 <- sqrt(colSums(sweep(sc1^2, 1, p, "*")))
      sc1 <- sweep(sc1, 2, phi1, "/")

      layer <- cbind(sweep(1 - layer, 2, phi, "*"), sweep(layer, 2, phi1, "*"))
      colnames(layer) <- c(paste("Dim", 1:nd, " S", sep=""), paste("Dim", 1:nd, " E", sep=""))
      rownames(layer) <- dimnames(tab)[[3]]

      sc <- array(cbind(sc, sc1), dim=c(nrow(sc), 2 * nd, 1), dimnames=list(rownames(tab), colnames(layer)))
  }
  else if(any(grepl("RCReg", names(coef(model))))) {
      sc <- array(sc, dim=c(nrow(sc), nd, nl)) +
                sweep(array(sc1, dim=c(nrow(sc1), nd, nl)), 3:2, layer, "*")

      phi <- sqrt(margin.table(sweep(sc^2, 1, p, "*"), 2:3))
      sc <- sweep(sc, 2:3, phi, "/")
      layer <- t(phi)

      # Order dimensions according to phi on first layer category
      ord <- order(abs(layer[1,]), decreasing=TRUE)
      layer <- layer[,ord, drop=FALSE]
      sc <- sc[,ord,, drop=FALSE]

      colnames(layer) <- colnames(sc) <- paste("Dim", 1:nd, sep="")
      rownames(sc) <- rownames(tab)
      rownames(layer) <- dimnames(sc)[[3]] <- dimnames(tab)[[3]]
  }
  else {
      stop("Invalid model")
  }


  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:dim(sc)[3]) {
      for(j in 1:nd) {
          if(sc[1,j,i] < 0)
              sc[,j,i] <- -sc[,j,i]
      }
  }

  # Scores can switch sides from one layer to another
  # Reverse axes for each layer state so that the scores are positively correlated to the first layer scores for rows
  for(i in seq_len(dim(sc)[3] - 1)) {
      chg <- ifelse(diag(as.matrix(cor(sc[,,i+1], sc[,,1]))) >= 0, 1, -1)
      sc[,,i+1] <- sweep(sc[,,i+1, drop=FALSE], 2, chg, "*")
  }

  if(length(dg) > 0) {
      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")
      rownames(dg) <- dimnames(tab)[[3]]
  }

  obj <- list(phi = layer, row = sc, col = sc, diagonal = dg,
              weighting = weighting, row.weights = p, col.weights = p)

  class(obj) <- c("assoc.rcL.dyn", "assoc.rcL", "assoc.symm", "assoc")
  obj
}

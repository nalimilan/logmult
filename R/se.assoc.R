
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

se.hmskewL <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "hmskew"))
      stop("x must be a hmskew object")

  if(length(x[["assoc"]]) > 0 && length(x$assoc.hmskew) > 0)
      return(list(assoc=se.assoc(x$assoc), assoc.hmskew=se.assoc(x$assoc.hmskewL)))
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


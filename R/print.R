print.rc <- function(x, ...) {
  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(round(ass$phi[1,], d=2))
  cat("\nNormalized row scores:\n")
  print(round(ass$row[,,1], d=2))
  cat("\nNormalized column scores:\n")
  print(round(ass$col[,,1], d=2))

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(round(ass$diag[1,], d=2))
  }

  cat("\nNormalization weights:", ass$weighting, "\n")
  cat("Deviance:", round(x$deviance, d=2), "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(extractAIC(x, k=log(sum(x$data)))[2]), " - ",
      "AIC:", round(extractAIC(x)[2]), "\n")
}

print.rc.homog <- function(x, ...) {
  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(round(ass$phi[1,], d=2))
  cat("\nNormalized scores:\n")
  print(round(ass$row[,,1], d=2))

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(round(ass$diag[1,], d=2))
  }

  cat("\nNormalization weights:", ass$weighting, "\n")
  cat("Deviance:", round(x$deviance, d=2), "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(extractAIC(x, k=log(sum(x$data)))[2]), " - ",
      "AIC:", round(extractAIC(x)[2]), "\n")
}

print.hmskew <- function(x, ...) {
  ass <- x$assoc

  if(length(ass) > 0) {
      cat("Intrinsic symmetric association coefficients:\n")
      print(round(ass$phi[1,], d=2))
      cat("\nNormalized symmetric association scores:\n")
      print(round(ass$row[,,1], d=2))
  }

  ass <- x$assoc.hmskew

  cat("Intrinsic skew association coefficients:\n")
  print(round(ass$phi[1,], d=2))
  cat("\nNormalized skew association scores:\n")
  print(round(ass$row[,,1], d=2))

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(round(ass$diag[1,], d=2))
  }

  cat("\nNormalization weights:", ass$weighting, "\n")
  cat("Deviance:", round(x$deviance, d=2), "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(extractAIC(x, k=log(sum(x$data)))[2]), " - ",
      "AIC:", round(extractAIC(x)[2]), "\n")
}

print.yrcskew <- function(x, ...) {
  ass <- x$assoc

  if(length(ass) > 0) {
      cat("Intrinsic symmetric association coefficients:\n")
      print(round(ass$phi[1,], d=2))
      cat("\nNormalized symmetric association scores:\n")
      print(round(ass$row[,,1], d=2))
  }

  ass <- x$assoc.yrcskew

  cat("Intrinsic skew association coefficients:\n")
  print(round(ass$phi[1,], d=2))
  cat("\nNormalized skew association scores:\n")
  print(round(ass$row[,,1], d=2))

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(round(ass$diag[1,], d=2))
  }

  cat("\nNormalization weights:", ass$weighting, "\n")
  cat("Deviance:", round(x$deviance, d=2), " - ",
      "Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(extractAIC(x, k=log(sum(x$data)))[2]), " - ",
      "AIC:", round(extractAIC(x)[2]), "\n")
}

print.rcL <- function(x, ...) {
  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(round(ass$phi, d=2))

  if(dim(ass$row)[3] == 1) {
      cat("\nNormalized row scores for all layers:\n")
      print(round(ass$row[,,1], d=2))
  }
  else {
      cat("\nNormalized row scores:\n")
      print(round(ass$row, d=2))
  }

  if(dim(ass$col)[3] == 1) {
      cat("\nNormalized column scores for all layers:\n")
      print(round(ass$col[,,1], d=2))
  }
  else {
      cat("\nNormalized column scores:\n")
      print(round(ass$col, d=2))
  }

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(round(ass$diag, d=2))
  }

  cat("\nNormalization weights:", ass$weighting, "\n")
  cat("Deviance:", round(x$deviance, d=2), "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(extractAIC(x, k=log(sum(x$data)))[2]), " - ",
      "AIC:", round(extractAIC(x)[2]), "\n")
}

print.rcL.homog <- function(x, ...) {
  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(round(ass$phi, d=2))

  if(dim(ass$row)[3] == 1) {
      cat("\nNormalized scores for all layers:\n")
      print(round(ass$row[,,1], d=2))
  }
  else {
      cat("\nNormalized scores:\n")
      print(round(ass$row, d=2))
  }

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(round(ass$diag, d=2))
  }

  cat("\nNormalization weights:", ass$weighting, "\n")
  cat("Deviance:", round(x$deviance, d=2), "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(extractAIC(x, k=log(sum(x$data)))[2]), " - ",
      "AIC:", round(extractAIC(x)[2]), "\n")
}

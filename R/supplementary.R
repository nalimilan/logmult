
sup.scores.rc <- function(model, tab, ass, rowsup, colsup,
                          symmetry=c("asymmetric", "symmetric", "skew-symmetric"), str="Mult") {
  symmetry <- match.arg(symmetry)
  stopifnot(!is.null(rowsup) || !is.null(colsup))
  stopifnot(symmetry =="asymmetric" || nrow(rowsup) == ncol(colsup))

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  phi <- ass$phi
  row <- ass$row
  col <- ass$col
  rp <- ass$row.weights
  cp <- ass$col.weights
  weighting <- ass$weighting

  # Model with association parameters constrained to their ajusted value
  args <- list()
  args$formula <- model$formula
  args$data <- tab
  args$family <- model$family
  args$tolerance <- model$tolerance
  args$verbose <- FALSE
  args$constrain <- which(grepl(str, names(parameters(model))))

  if(symmetry != "asymmetric")
      args$constrainTo <- sweep(row[,, 1], 2, sqrt(phi[1,]), "*")
  else
      args$constrainTo <- sweep(rbind(row[,, 1], col[,, 1]), 2, sqrt(phi[1,]), "*")

  model2 <- do.call("gnm", args)

  # Diagonal does not make sense in this model since table is not square,
  # or has NA diagonal if both rowsup and colsup have been supplied
  if(grepl("Diag", deparse(args$formula)))
      args$formula <- as.formula(paste(deparse(args$formula),
                                       sprintf("- Diag(%s, %s)", vars[1], vars[2])))

      if(is.null(colsup)) {
          args$data <- as.table(rowsup)
      }
      else if(is.null(rowsup)) {
          args$data <- as.table(colsup)
      }
      else {
          # Block matrix that is symmetric if rowsup and colsup are transposes of each other
          args$data <- as.table(cbind(rbind(matrix(NA, nrow(colsup), ncol(rowsup)), rowsup),
                                      rbind(colsup, matrix(NA, nrow(rowsup), ncol(colsup)))))
          dimnames(args$data) <- list(c(rownames(model$data), rownames(rowsup)),
                                      c(colnames(model$data), colnames(colsup)))
      }

      names(dimnames(args$data)) <- names(dimnames(model$data))

      if(symmetry == "skew-symmetric") {
          hmterm <- sprintf("+ HMSkew(%s, %s)", vars[1], vars[2])

          args2 <- args
          args2$constrain <- NULL
          args2$constrainTo <- NULL
          args2$formula <- as.formula(sub(paste("\\Q", hmterm, "\\E", sep=""),
                                          "", deparse(args$formula)))

          base <- do.call("gnm", args2)
          args$start <- c(rep(NA, length(parameters(base))), residEVD(base, 1, skew=TRUE))

          hmnames <- names(gnm:::gnmTools(gnm:::gnmTerms(as.formula(paste("Freq ~ -1", hmterm)),
                                                         data=args$data),
                                          as.data.frame(args$data), FALSE)$start)

          cnames <- c(names(parameters(base)), hmnames)
      }
      else {
          cnames <- names(gnm:::gnmTools(gnm:::gnmTerms(args$formula, data=args$data),
                                         as.data.frame(args$data), FALSE)$start)
      }

      args$constrain <- which(cnames %in% names(parameters(model2)))
      args$constrainTo <- parameters(model2)[match(cnames, names(parameters(model2)), nomatch=0)]
      msup <- do.call("gnm", args)

      sup <- parameters(msup)[setdiff(pickCoef(msup, str), args$constrain)]

      if(symmetry != "asymmetric")
          dim(sup) <- c(sum(nrow(rowsup)), ncol(phi))
      else
          dim(sup) <- c(sum(nrow(rowsup), ncol(colsup)), ncol(phi))

      sup <- sweep(sup, 2, sqrt(phi[1,]), "/")


  if(length(rowsup) > 0) {
      if(weighting == "none")
          rpsup <- rep(1, nrow(rowsup))
      else if(weighting == "uniform")
          rpsup <- rep(1/nrow(rowsup), nrow(rowsup))
      else
          rpsup <- prop.table(apply(rowsup, 1, sum, na.rm=TRUE))

      rp <- c(rp, rpsup)

      row2 <- rbind(cbind(row[,,1]), sup[seq(nrow(rowsup)), , drop=FALSE])
      dim(row2)[3] <- 1
      names(dimnames(row2)) <- names(dimnames(row))
      rownames(row2) <- c(rownames(row), rownames(rowsup))
      colnames(row2) <- colnames(row)
      row <- row2
  }

  if(length(colsup) > 0) {
      if(weighting == "none")
          cpsup <- rep(1, ncol(colsup))
      else if(weighting == "uniform")
          cpsup <- rep(1/ncol(colsup), ncol(colsup))
      else
          cpsup <- prop.table(apply(colsup, 2, sum, na.rm=TRUE))

      cp <- c(cp, cpsup)

      if(symmetry != "asymmetric") {
          col <- row
          rp <- cp <- (rp + cp)/2
      }
      else {
          col2 <- rbind(cbind(col[,,1]), sup[-seq(nrow(rowsup)), , drop=FALSE])
          dim(col2)[3] <- 1
          names(dimnames(col2)) <- names(dimnames(col))
          rownames(col2) <- c(rownames(col), colnames(colsup))
          colnames(col2) <- colnames(col)
          col <- col2
      }
  }

  list(row=row, col=col, row.weights=rp, col.weights=cp)
}

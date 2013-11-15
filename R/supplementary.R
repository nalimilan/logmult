
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
  weighting <- ass$weighting

  args <- list()
  args$formula <- model$formula
  args$family <- model$family
  args$tolerance <- model$tolerance
  args$verbose <- FALSE


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

  args$constrain <- which(cnames %in% names(parameters(model)) & grepl(str, cnames))

  if(symmetry != "asymmetric")
      args$constrainTo <- sweep(cbind(row[,, 1]), 2, sqrt(phi[1,]), "*")
  else if(length(rowsup) > 0 && length(colsup) > 0)
      args$constrainTo <- sweep(rbind(cbind(row[,, 1]), cbind(col[,, 1])), 2, sqrt(phi[1,]), "*")
  else if(length(rowsup) > 0)
      args$constrainTo <- sweep(cbind(col[,, 1]), 2, sqrt(phi[1,]), "*")
  else
      args$constrainTo <- sweep(cbind(row[,, 1]), 2, sqrt(phi[1,]), "*")

  msup <- do.call("gnm", args)

  sup <- parameters(msup)[setdiff(pickCoef(msup, str), args$constrain)]

  if(symmetry != "asymmetric") {
      dim(sup) <- c(sum(nrow(rowsup)), ncol(phi))

      if(weighting == "none")
          rpsup <- rep(1, nrow(rowsup))
      else if(weighting == "uniform")
          rpsup <- rep(1/nrow(rowsup), nrow(rowsup))
      else
          rpsup <- prop.table(apply(rowsup, 1, sum, na.rm=TRUE) + apply(colsup, 2, sum, na.rm=TRUE))

      rsup <- sup[seq(nrow(rowsup)), , drop=FALSE]

      rsup <- sweep(rsup, 2, colSums(sweep(rsup, 1, rpsup/sum(rpsup), "*")), "-")
      rsup <- sweep(rsup, 2, sqrt(phi[1,]), "/")

      row2 <- rbind(cbind(row[,,1]), rsup)
      dim(row2)[3] <- 1
      names(dimnames(row2)) <- names(dimnames(row))
      rownames(row2) <- c(rownames(row), rownames(rowsup))
      colnames(row2) <- colnames(row)
      row <- col <- row2
  }
  else {
      dim(sup) <- c(sum(nrow(rowsup), ncol(colsup)), ncol(phi))

      if(length(rowsup) > 0) {
          if(weighting == "none")
              rpsup <- rep(1, nrow(rowsup))
          else if(weighting == "uniform")
              rpsup <- rep(1/nrow(rowsup), nrow(rowsup))
          else
              rpsup <- prop.table(apply(rowsup, 1, sum, na.rm=TRUE))

          rsup <- sup[seq(nrow(rowsup)), , drop=FALSE]

          rsup <- sweep(rsup, 2, colSums(sweep(rsup, 1, rpsup/sum(rpsup), "*")), "-")
          rsup <- sweep(rsup, 2, sqrt(phi[1,]), "/")

          row2 <- rbind(cbind(row[,,1]), rsup)
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

          csup <- sup[seq(NROW(rowsup) + 1, nrow(sup)), , drop=FALSE]

          csup <- sweep(csup, 2, colSums(sweep(csup, 1, cpsup/sum(cpsup), "*")), "-")
          csup <- sweep(csup, 2, sqrt(phi[1,]), "/")

          col2 <- rbind(cbind(col[,,1]), csup)
          dim(col2)[3] <- 1
          names(dimnames(col2)) <- names(dimnames(col))
          rownames(col2) <- c(rownames(col), colnames(colsup))
          colnames(col2) <- colnames(col)
          col <- col2
      }
  }

  if(length(dim(tab)) == 3) {
      row.weights.sup <- apply(rowsup, c(1, 3), sum, na.rm=TRUE)
      col.weights.sup <- apply(colsup, c(2, 3), sum, na.rm=TRUE)
  }
  else {
      row.weights.sup <- as.matrix(apply(rowsup, 1, sum, na.rm=TRUE))
      col.weights.sup <- as.matrix(apply(colsup, 2, sum, na.rm=TRUE))
  }

  list(row=row, col=col,
       row.weights=rbind(ass$row.weights, row.weights.sup),
       col.weights=rbind(ass$col.weights, col.weights.sup))
}

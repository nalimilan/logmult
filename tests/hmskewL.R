# Artificial example to check the hmskewL model

library(logmult)
data(ocg1973)

tab <- array(ocg1973, dim=c(nrow(ocg1973), ncol(ocg1973), 2))

model <- hmskewL(tab[5:1, 5:1,], weighting="uniform")
ass <- model$assoc

# First score for Farmers is slightly different from the original article
stopifnot(isTRUE(all.equal(round(ass$row[,,1] * sqrt(ass$phi[1,1]), d=2)[5:1,],
                                 matrix(c(-0.08, -0.2, -0.23, -0.11, 0.61,
                                           0.34,  0.3, -0.13, -0.51, 0), 5, 2),
                           check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(round(ass$row[,,1] * sqrt(ass$phi[2,1]), d=2)[5:1,],
                                 matrix(c(-0.08, -0.2, -0.23, -0.11, 0.61,
                                           0.34,  0.3, -0.13, -0.51, 0), 5, 2),
                           check.attributes=FALSE)))

model2 <- hmskewL(tab[5:1, 5:1,], weighting="uniform", layer.effect.skew="heterogeneous")
stopifnot(isTRUE(all.equal(model$assoc$phi, model2$assoc$phi)))
stopifnot(isTRUE(all.equal(model$assoc$row[,,1], model2$assoc$row[,,1])))

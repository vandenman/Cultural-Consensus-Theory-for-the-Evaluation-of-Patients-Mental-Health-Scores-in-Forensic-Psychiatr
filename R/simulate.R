getProbsOrderedR <- function(delta, location, scale, cdf = plogis) {

  # replicate of what Stan does
  nc <- length(delta) + 1L
  prob <- numeric(nc)

  vals <- cdf(location - delta, 0, scale)
  prob[1] <- 1 - vals[1]
  for (c in 2:(nc - 1)) {
    prob[c] <- vals[c - 1L] - vals[c]
  }
  prob[nc] <- vals[nc - 1]

  return(prob)
}

getProbsOrderedR2 <- function(delta, mu = 0, s = 1, cdf = plogis) {

  # As written in AB.
  nc <- length(delta) + 1L
  probs <- numeric(nc)

  probs[1] <- cdf(delta[1], mu, s)
  for (c in 2:(nc - 1)) {
    probs[c] <- cdf(delta[c], mu, s) - sum(probs[1:(c - 1)])
  }
  probs[nc] <- 1 - sum(probs[1:(nc - 1)])
  return(probs)
}

dOrdered <- function(x, delta, location = 0, scale = 1, log = TRUE, cdf = plogis) {
  prob <- getProbsOrderedR(delta, location, scale, cdf = cdf)
  if (log)
    return(log(prob[x]))
  else return(prob[x])
}

simulate <- function(nr, ni, nc, np = 1, nl = 1,
                     ncatR = 1, ncatP = 1) {

  # constants
  sdLt <- 0.4
  sdBeta <- 1
  sdEta <- 1
  sdReg <- .4
  muKappa <- .25
  sdKappa <- 0.5
  muAlpha <- 1
  sdAlpha <- 0.25
  muZeta  <- 1
  sdZeta  <- 0.25

  # rater influence
  designR <- sort(rep_len(1:ncatR, length.out = nr))
  muBeta <- rnorm(ncatR, 0, .3)
  muBeta <- muBeta - mean(muBeta)

  # patient latent parametersinfluence
  designP <- sort(rep_len(1:ncatP, length.out = np))
  muEta <- matrix(rnorm(nl * ncatP), nl, ncatP)

  # same loadings for all patients
  reg <- abs(rnorm(ni, 0, sdReg))
  # which item loads on which latent variable?
  latentMembership <- sort(rep_len(1:nl, length.out = ni))

  eta <- matrix(rnorm(nl*np, muEta[, designP], sdEta), nl, np)

  for (l in 1:nl) {
    muEta[l, ] <- muEta[l, ] - mean(muEta[l, ])
    eta[l, ] <- eta[l, ] - mean(eta[l, ])
    if (np > 1L)
      eta[l, ] <- eta[l, ] / sd(eta[l, ])
  }

  # put the item with the highest loading first and force it to be positive
  # reg <- unlist(tapply(reg, latentMembership, function(x) x[order(abs(x), decreasing = TRUE)]))
  # idxPos <- match(1:nl, latentMembership)
  # reg[idxPos] <- abs(reg[idxPos])

  lt <- kappa <- matrix(NA, ni, np)
  for (p in seq_len(np)) {
    muLt <- reg * (eta[, p] + muEta[, designP[p]])[latentMembership]
    lt[, p] <- rnorm(ni, muLt, sdLt)
    kappa[, p] <- abs(rnorm(ni, muKappa, sdKappa))
    kappa[, p] <- kappa[, p] / mean(kappa[, p])
  }
  lt <- lt

  # betas: fix mean and variance to observed
  beta <- rnorm(nr, muBeta[designR], sdBeta)
  # obs <- do.call(rbind, tapply(beta, designR, function(x) c(mean(x), sd(x))))
  # tru <- cbind(muBeta, sdBeta)
  # plot(beta, col = designR, pch = 16)
  # abline(h = muBeta, col = unique(designR))
  # abline(h = obs[, 1], col = unique(designR), lty = 2)
  # legend("topleft", c("observed", "true"), col = 1, lty = 1:2, bty = "n")

  alpha <- abs(rnorm(nr, muAlpha, sdAlpha))
  zeta  <- abs(rnorm(nr, muZeta,  sdZeta))

  # identification constraint: mean = 0 (unbounded parameters) or mean = 1 (positive parameters)
  # beta  <- beta  - mean(beta)
  alpha <- alpha / mean(alpha)
  zeta  <- zeta  / mean(zeta)

  gamma <- (1:(nc - 1))
  gamma <- log(gamma / (nc - gamma))

  x <- array(NA_integer_, c(nr, ni, np))
  for (r in seq_len(nr)) {

    delta <- alpha[r] * gamma + beta[r]
    for (i in seq_len(ni)) {
      for (p in seq_len(np)) {

        location <- lt[i, p]
        scale    <- kappa[i, p] / zeta[r]

        prob <- getProbsOrderedR(delta, location, scale)
        # browser()
        # prob == dOrderedStan(location / scale, delta / scale)
        x[r, i, p] <- sample(1:nc, 1L, prob = prob)

      }
    }
  }

  output <- list(
    x = x, alpha = alpha, beta = beta,
    kappa = kappa, zeta = zeta,
    eta = eta, reg = reg, latentMembership = latentMembership,
    muBeta = muBeta, sdBeta = sdBeta,
    muEta = muEta, sdEta = sdEta,
    designR = designR, designP = designP,
    lt = lt, nr = nr, ni = ni, nc = nc, np = np, nl = nl, ncatR = ncatR, ncatP = ncatP
  )
  output <- lapply(output, drop)
  class(output) <- "CCT_sim"
  return(output)

}

reshapeSim <- function(x) {
  nr <- x$nr
  ni <- x$ni
  np <- x$np
  values <- x$x
  dimnames <- list()
  nms <- list(
    "Rater"   = seq_len(nr),
    "Item"    = seq_len(ni),
    "Patient" = seq_len(np)
  )
  nms <- nms[lengths(nms) > 1L]
  dimnames(values) <- nms
  df <- as.data.frame.table(values, responseName = "Outcome")
  if (x$ncatR > 1L)
    df$ncatR <- x$designR[df$Rater]
  if (x$ncatP > 1L)
    df$ncatP <- x$designP[df$Patient]

  return(df)
}

obsList2datList <- function(obsList) {

  if (obsList[["np"]] == 1L) { # LTRM
    return(obsList[c("x", "nr", "ni", "nc")])
  } else { # Extended LTRM
    datList <- obsList[c("x", "nr", "ni", "np", "nc", "nl", "ncatR", "ncatP", "latentMembership",
                         "designR", "designP")]
    datList[["firstItemOfLv"]] <- match(
      seq_len(datList[["nl"]]),
      obsList[["latentMembership"]]
    )
    return(datList)
  }
}

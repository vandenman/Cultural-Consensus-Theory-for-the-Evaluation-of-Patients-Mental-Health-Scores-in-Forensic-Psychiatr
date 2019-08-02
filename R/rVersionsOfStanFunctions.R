getProbsOrderedcdfR <- function(delta, location, scale, cdf = plogis) {

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

dOrderedLogistic <- function(x, delta, location = 0, scale = 1, log = TRUE) {
  prob <- getProbsOrderedcdfR(delta, location, scale)
  if (log)
    return(log(prob[x]))
  else return(prob[x])
}

dOrderedProbit <- function(x, delta, location = 0, scale = 1, log = TRUE) {
  prob <- getProbsOrderedcdfR(delta, location, scale, pnorm)
  if (log)
    return(log(prob[x]))
  else return(prob[x])
}

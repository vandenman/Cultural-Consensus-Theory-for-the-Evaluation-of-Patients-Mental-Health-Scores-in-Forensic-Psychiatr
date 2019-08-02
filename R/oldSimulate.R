
simulate1 <- function(nr, ni, nc) {
  
  lt <- rbeta(ni, 1, 1)
  lambda <- 1:nc
  lambda2 <- (1 - lambda) / lambda
  lambda2 <- 1 / (1 + lambda2)
  
  # y <- 1.0 / (1.0 + exp(plogis(lt)))
  y <- lt
  
  prob <- matrix(NA, ni, nc)
  
  # sdD <- 1
  x <- matrix(NA, nr, ni)
  for (i in seq_len(ni)) {
    
      prob[i, 1] <- 1 - pnorm(y[i] - lambda2[1]) # 1 - sum(prob[i, 1:(nc - 1)])
			
			for (c in 2:(nc - 1)) {
				# P(x = c)    =  P(x <= c)                        - [P(x = c-1) + P(x = c-2) + ... + P(x = 1)]
				prob[i, c] <- pnorm(y[i] - lambda2[c - 1L]) - pnorm(y[i] - lambda2[c])
			}
      prob[i, nc] <- pnorm(y[i] - lambda2[nc - 1L])
			# prob[i, nc] <- 1 - sum(prob[i, -nc]) # pnorm(y[i] - lambda2[nc])
			
			x[, i] <- sample(1:nc, nr, prob = prob[i, ], replace = TRUE)
  }
  
  
  return(list(
		x = x, prob = prob,
		lt = lt, y = y,
		# d = d, a = a, b = b,
		# delta = delta, lambda = lambda,
		# beta = beta,
		# design = design, y_cov = y_cov,
		nr = nr, ni = ni, nc = nc
	))
  
}

simulate2 <- function(nr, ni, nc) {

  lt <- rbeta(ni, 1, 1)
  alpha <- rgamma(nr, 1, 1)
  beta <- rgamma(nr, 1, 1)
  lambda <- (1:nc) / nc
  lambda2 <- (1 - lambda) / lambda
  # lambda2 <- 1 / (1 + lambda2)

  y <- ilogit(lt)

  prob <- matrix(NA, ni, nc)

  sdD <- 1
  x <- matrix(NA, nr, ni)
  for (i in seq_len(ni)) {
    for (r in seq_len(nr)) {

      delta <- 1 / (1 + lambda2^alpha[r] / beta[r])
      
      prob[i, 1] <- 1 - pnorm(y[i] - delta[1], sd = sdD) # 1 - sum(prob[i, 1:(nc - 1)])
      
      for (c in 2:(nc - 1)) {
        # P(x = c)  = P(x <= c_{i-1})                         - P(x <= c_{i})
        prob[i, c] <- pnorm(y[i] - delta[c - 1L], sd = sdD) - pnorm(y[i] - delta[c], sd = sdD)
      }
      prob[i, nc] <- pnorm(y[i] - delta[nc - 1L], sd = sdD)
      if (anyNA(prob[i, ]))
        browser()
      x[r, i] <- sample(1:nc, 1L, prob = prob[i, ])
    }
  }

  return(list(
		x = x, prob = prob, alpha = alpha, beta = beta,
		lt = lt, y = y,
		# a = a, b = b,
		# delta = delta, lambda = lambda,
		# beta = beta,
		# design = design, y_cov = y_cov,
		nr = nr, ni = ni, nc = nc
	))
  
}

simulate3 <- function(nr, ni, nc) {

  lt <- rnorm(ni, 0, 1)
  lt <- lt - mean(lt)
  alpha <- rgamma(nr, 1, 1)
  beta <- rnorm(nr, 0, 1)
  lambda <- (1:(nc - 1)) / nc
  lambda <- log(lambda / (1 - lambda))
  print(lambda)
  
  y <- lt

  prob <- matrix(NA, ni, nc)

  sdD <- 1
  x <- matrix(NA, nr, ni)
  for (r in seq_len(nr)) {
    
    delta <- alpha[r] * lambda + beta[r]
    for (i in seq_len(ni)) {
      prob[i, 1] <- 1 - plogis(y[i] - delta[1], scale = sdD) # 1 - sum(prob[i, 1:(nc - 1)])
      
      for (c in 2:(nc - 1)) {
        # P(x = c)  = P(x <= c_{i-1})                         - P(x <= c_{i})
        prob[i, c] <- plogis(y[i] - delta[c - 1L], scale = sdD) - plogis(y[i] - delta[c], scale = sdD)
      }
      prob[i, nc] <- plogis(y[i] - delta[nc - 1L], scale = sdD)
      x[r, i] <- sample(1:nc, 1L, prob = prob[i, ])
    }
  }
  
  if (any(abs(rowSums(prob) - 1) > sqrt(.Machine$double.eps)))
    browser()

  return(list(
		x = x, prob = prob, alpha = alpha, beta = beta,
		lt = lt,
		# a = a, b = b,
		# delta = delta, lambda = lambda,
		# beta = beta,
		# design = design, y_cov = y_cov,
		nr = nr, ni = ni, nc = nc
	))
  
}

simulate5 <- function(nr, ni, nc, ncat) {

  lt <- rnorm(ni, 0, 1)
  lt <- lt - mean(lt)
  alpha <- rgamma(nr, 1, 1)
  beta <- rnorm(nr, 0, 1)
  lambda <- (1:(nc - 1)) / nc
  lambda <- log(lambda / (1 - lambda))
  
  design <- sample(ncat, nr, TRUE)
  catEffect <- rnorm(ncat)
  
  error <- rnorm(nr, sd = 1)
  y <- matrix(error, nrow = ni, ncol = nr, TRUE)
  y <- sweep(y, 1, lt, `+`)

  prob <- matrix(NA, ni, nc)

  sdD <- 1
  x <- matrix(NA, nr, ni)
  for (r in seq_len(nr)) {
    
    delta <- alpha[r] * lambda + beta[r]
    for (i in seq_len(ni)) {
      prob[i, 1] <- 1 - plogis(y[i, r] - delta[1] + catEffect[design[r]], scale = sdD) # 1 - sum(prob[i, 1:(nc - 1)])
      
      for (c in 2:(nc - 1)) {
        # P(x = c)  = P(x <= c_{i-1})                         - P(x <= c_{i})
        prob[i, c] <- plogis(y[i, r] - delta[c - 1L] + catEffect[design[r]], scale = sdD) - 
          plogis(y[i, r] - delta[c] + catEffect[design[r]], scale = sdD)
      }
      prob[i, nc] <- plogis(y[i, r] - delta[nc - 1L] + catEffect[design[r]], scale = sdD)
      x[r, i] <- sample(1:nc, 1L, prob = prob[i, ])
    }
  }
  
  if (any(abs(rowSums(prob) - 1) > sqrt(.Machine$double.eps)))
    browser()

  return(list(
		x = x, prob = prob, alpha = alpha, beta = beta,
		# kappa = kappa, zeta = zeta,
		lt = lt, y = y, error = error, 
		design = design, catEffect = catEffect,
		nr = nr, ni = ni, nc = nc, ncat = ncat
	))
  
}

simulate6 <- function(nr, ni, nc, ncat, np) {
  
  lt <- matrix(rnorm(ni*np, 0, 1), ni, np)
  for (p in seq_along(np))
    lt[, p] <- lt[, p] - mean(lt[, p])
  # lt <- lt - mean(lt)
  alpha <- rgamma(nr, 1, 1)
  beta <- rnorm(nr, 0, 1)
  lambda <- (1:(nc - 1)) / nc
  lambda <- log(lambda / (1 - lambda))
  
  design <- sample(ncat, nr, TRUE)
  catEffect <- rnorm(ncat)
  
  sdError <- abs(rnorm(nr))
  
  y <- error <- array(NA_real_, c(ni, nr, np))
  for (p in 1:np) {
    y[, , p] <- error[, , p] <- matrix(rnorm(ni*nr, sd = rep(sdError, each=ni)), ni, nr)
    y[, , p] <- sweep(y[, , p], 1, lt[, p], `+`)
  }
  # browser()
  # plot(
  #   matrixStats::colSds(error[, , p]),
  #   sdError
  # )

  prob <- matrix(NA, ni, nc)
  
  sdD <- 1
  x <- array(NA_integer_, c(nr, ni, np))
  for (p in seq_len(np)) {
    for (r in seq_len(nr)) {
      
      delta <- alpha[r] * lambda + beta[r]
      for (i in seq_len(ni)) {
        prob[i, 1] <- 1 - plogis(y[i, r, p] - delta[1] + catEffect[design[r]], scale = sdD) # 1 - sum(prob[i, 1:(nc - 1)])
        
        for (c in 2:(nc - 1)) {
          # P(x = c)  = P(x <= c_{i-1})                         - P(x <= c_{i})
          prob[i, c] <- plogis(y[i, r, p] - delta[c - 1L] + catEffect[design[r]], scale = sdD) - 
            plogis(y[i, r, p] - delta[c] + catEffect[design[r]], scale = sdD)
        }
        prob[i, nc] <- plogis(y[i, r, p] - delta[nc - 1L] + catEffect[design[r]], scale = sdD)
        x[r, i, p] <- sample(1:nc, 1, TRUE, prob = prob[i, ])
      }
    }
  }
  
  if (any(abs(rowSums(prob) - 1) > sqrt(.Machine$double.eps)))
    browser()
  
  return(list(
    x = x, prob = prob, alpha = alpha, beta = beta,
    lt = lt, y = y, error = error, sdError = sdError,
    design = design, catEffect = catEffect,
    nr = nr, ni = ni, nc = nc, ncat = ncat
  ))
  
}


simulateModel1 <- function(nr, ni, nc) {

  lt <- rnorm(ni, 0, 1)
  beta <- rnorm(nr, 0, 1)

  alpha <- abs(rnorm(nr, 1, .25))
  kappa <- abs(rnorm(ni, 1, .25))
  zeta  <- abs(rnorm(nr,  1, .25))
  
  alpha <- alpha / mean(alpha)
  kappa <- kappa / mean(kappa)
  zeta  <- zeta  / mean(zeta)
  
  gamma <- (1:(nc - 1))
  gamma <- log(gamma / (nc - gamma))
  
  x <- matrix(NA, nr, ni)
  for (r in seq_len(nr)) {

    delta <- alpha[r] * gamma + beta[r]
    for (i in seq_len(ni)) {
      
      location <- lt[i]
      scale    <- kappa[i]# / zeta[r]

      prob <- getProbsOrderedLogisticR(delta, location, scale)
      x[r, i] <- sample(1:nc, 1L, prob = prob)

    }
  }
  
  return(list(
		x = x, prob = prob, alpha = alpha, beta = beta,
		kappa = kappa, zeta = zeta,
		lt = lt, nr = nr, ni = ni, nc = nc
	))
  
}

simulateModel2 <- function(nr, ni, nc, np, nl) {

  eta <- matrix(rnorm(nl*np), nl, np)
  reg <- rnorm(ni, 0)
  # index i says that item i loads on the latent variable latentMembership[i]
  latentMembership <- sort(rep_len(1:nl, 20))

  # put the item with the highest loading first and force it to be positive
  # this avoids label switching problems
  reg <- unlist(tapply(reg, latentMembership, function(x) x[order(abs(x), decreasing = TRUE)]))
  idxPos <- match(1:nl, latentMembership)
  reg[idxPos] <- abs(reg[idxPos])

  # identification constraint: sign of product of loadings is positive -- failed?
  # signs <- tapply(reg, latentMembership, function(x) sign(prod(x)))
  # reg2 <- signs[3] * reg[latentMembership==3]
  # tapply(reg2, latentMembership, function(x) sign(prod(x)))
  lt <- kappa <- matrix(NA, ni, np)
  for (p in seq_len(np)) {
    lt[, p] <- reg * eta[latentMembership, p] + rnorm(ni, 0, .5)
    kappa[, p] <- abs(rnorm(ni, .25, .5))
    kappa[, p] <- kappa[, p] / mean(kappa[, p])
  }

  beta <- rnorm(nr, 0, 1)
  alpha <- abs(rnorm(nr, 1, .25))
  zeta  <- abs(rnorm(nr,  1, .25))

  # identification constraint: mean = 0 (unbounded parameters) or mean = 1 (positive parameters)
  beta  <- beta  - mean(beta)
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
        scale    <- kappa[i, p]# / zeta[r]
        
        prob <- getProbsOrderedLogisticR(delta, location, scale)
        x[r, i, p] <- sample(1:nc, 1L, prob = prob)
        
      }
    }
  }
  
  return(list(
		x = x, prob = prob, alpha = alpha, beta = beta,
		kappa = kappa, zeta = zeta,
		eta = eta, reg = reg, latentMembership = latentMembership,
		lt = lt, nr = nr, ni = ni, nc = nc, np = np
	))
  
}

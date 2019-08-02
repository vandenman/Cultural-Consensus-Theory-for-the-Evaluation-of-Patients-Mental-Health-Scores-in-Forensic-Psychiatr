DIM <- function(x) c(NROW(x), NCOL(x))

genData <- function(nr, ni, nc = 5, design = matrix(0, ni, 1), beta = 0, y_cov = matrix(0, nr, 1)) {
	
	design <- as.matrix(design)
	y_cov  <- as.matrix(y_cov)
	ny     <- NCOL(y_cov)
	stopifnot(
		all.equal(DIM(design), c(ni, ny)),
	  all.equal(DIM(beta),   c(ny, 1L)),
		all.equal(DIM(y_cov),  c(nr, ny))
	)
	
	muT <- .5
	sdT <- 2
	muA <- muB <- 1
	sdA <- sdB <- .1
	sdD <- .1
	
	lt <- rbeta(ni, muT*sdT, (1.0 - muT)*sdT) # latent truth
	# d <- matrix(rnorm(nr*ni, 0, sdD), nr, ni) # discrepancy
	d <- matrix(0, nr, ni) # discrepancy
	y <- matrix(NA, nr, ni)
	logitLt <- qlogis(lt)
	for (i in seq_len(nr)) {
		# y[i, ] <- lt#ilogit(logit(lt) + d[i, ])             # latent appraisal
		# logistic_cdf(d[r, k] + logit(lt[k]), 0, 1);
	  y[i, ] <- plogis(logitLt, -d[i, ], 1)#ilogit(logit(lt) + d[i, ])             # latent appraisal
	}
	
	a <- rgamma(n = nr, shape = (muA / sdA)^2, rate = muA / sdA^2)      # rater scale
	b <- rgamma(n = nr, shape = (muB / sdB)^2, rate = muB / sdB^2)      # rater shift
	lambda <- seq(1 / nc, 1, 1 / nc)                     # unbiased thresholds
	lambda2 <- (1 - lambda) / lambda
	
	delta <- matrix(NA, nr, nc)
	for (i in 1:nr) {
		for (c in 1:nc) {
		  delta[i, c] <- b[i] / (b[i] + lambda2[c]^a[i]) # biased thresholds
		}
	}

	prob <- array(NA, dim = c(nr, ni, nc))
	x <- matrix(NA, nr, ni)
	for (i in 1:nr) {
		for (k in 1:ni) {
			
			mu <- beta %*% (design[k, ] * y_cov[i, ])
			prob[i, k, 1] <- pnorm(delta[i, 1], y[i, k] + mu, sdD)
			# browser()
			for (c in 2:(nc - 1)) {
				# P(x = c)    =  P(x <= c)                        - [P(x = c-1) + P(x = c-2) + ... + P(x = 1)]
				prob[i, k, c] <- pnorm(delta[i, c], y[i, k] + mu, sdD) - sum(prob[i, k, 1:(c - 1)])
			}
			prob[i, k, nc] <- 1 - sum(prob[i, k, 1:(nc - 1)])
			
			# x[i, k] <- rcat(1, prob[i, k, 1:nc])
			x[i, k] <- sample(1:nc, 1, prob = prob[i, k, 1:nc])
		}
	}
	
	return(list(
		x = x, prob = prob,
		lt = lt, d = d, y = y, a = a, b = b,
		delta = delta, lambda = lambda,
		beta = beta,
		design = design, y_cov = y_cov,
		nr = nr, ni = ni, nc = nc
	))
	
}

corplot <- function(x, y) {
	
	cc <- cor(x, y)
	if (!is.na(cc)) {
		main <- sprintf("Cor: %.3f", cor(x, y))
	} else {
		main <- "Cor: NA"
	}
	plot(x, y, main = main, xlab = "")
	abline(0, 1)
}

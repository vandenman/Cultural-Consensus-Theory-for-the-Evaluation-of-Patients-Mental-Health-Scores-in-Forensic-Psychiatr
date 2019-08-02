DIM <- function(x) c(NROW(x), NCOL(x))

genData <- function(nr, ni, nc = 5, design = matrix(0, ni, 1), beta = 0, y_cov = matrix(0, nr, 1)) {
	
	design <- as.matrix(design)
	y_cov <- as.matrix(y_cov)
	ny <- NCOL(y_cov)
	stopifnot(
		all.equal(DIM(design), c(ni, ny)),
	  	all.equal(DIM(beta), c(ny, 1)),
		all.equal(DIM(y_cov), c(nr, ny))
	)
	
	muT <- .5
	sdT <- 2
	muA <- 1
	muB <- 1
	sdA <- sdB <- 1
	sdD <- .1
	
	lt <- rbeta(ni, muT*sdT, (1.0 - muT)*sdT)            # latent truth
	d <- matrix(rnorm(nr*ni, 0, sdD), nr, ni)            # discrepancy
	y <- matrix(NA, nr, ni)
	for (i in seq_len(nr))
		y[i, ] <- lt#ilogit(logit(lt) + d[i, ])             # latent appraisal
	
	a <- rgamma(n = nr, shape = (muA / sdA)^2, rate = muA / sdA^2)      # rater scale
	b <- rgamma(n = nr, shape = (muB / sdB)^2, rate = muB / sdB^2)      # rater shift
	lambda <- seq(1 / nc, 1, 1 / nc)                     # unbiased thresholds
	
	delta <- matrix(NA, nr, nc)
	for (i in 1:nr) {
		for (c in 1:nc) {
			delta[i, c] <- b[i]*lambda[c]^a[i] / (b[i]*lambda[c]^a[i] + (1-lambda[c])^a[i]) # biased thresholds
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
	plot(x, y, main = main)
	abline(0, 1)
}

library(rstan)
## Simple ----
nr = 50
ni = 20
nc = 5

test <- genData(nr, ni, nc)
datList <- list(
	# data
	x = test$x,
	nr = nr,
	ni = ni,
	nc = nc,
	# hyperparameters
	muT = .5,
	sdT = 2,
	muA = 1,
	muB = 1,
	sdA = 1,
	sdB = 1,
	sdD = .1
)
res <- stan(file = "stanTranslation.stan", data = datList, chains = 2)
s <- summary(res)
ss <- s$summary
rnms <- rownames(ss)
aest <- ss[grepl("a[", rnms, fixed = TRUE), "mean"]
best <- ss[grepl("b[", rnms, fixed = TRUE), "mean"]
ltest <- ss[grepl("lt[", rnms, fixed = TRUE), "mean"]

# looks good
par(bty = "n", las = 1)
corplot(test$a, aest)
corplot(test$b, best)
corplot(test$lt, ltest)

## Covariates ----

nr = 50
ni = 20
nc = 5

y_cov <- matrix(rep(0:1, each = floor(nr/2)), nr, 1)
design <- matrix(rep(0:1, each = floor(ni/2)), ni, 1)
beta <- 1

test <- genData(nr, ni, nc, design = design, beta = beta, y_cov = y_cov)
datList <- list(
	# data
	x = test$x,
	nr = nr,
	ni = ni,
	nc = nc,
	
	nycat = 1,
	design = design,
	covariates = y_cov,
	
	# hyperparameters
	muT = .5,
	sdT = 2,
	muA = 1,
	muB = 1,
	sdA = 1,
	sdB = 1,
	sdD = .1
)
res <- stan(file = "withCovariates.stan", data = datList, chains = 2)
s <- summary(res)
ss <- s$summary

rmns <- rownames(ss)
aest <- ss[grepl("a[", rnms, fixed = TRUE), "mean"]
best <- ss[grepl("b[", rnms, fixed = TRUE), "mean"]
ltest <- ss[grepl("lt[", rnms, fixed = TRUE), "mean"]
betaest <- ss[nrow(ss)-1, "mean"]

# looks good
par(bty = "n", las = 1)
corplot(test$a, aest)
corplot(test$b, best)
corplot(test$lt, ltest[-21])
corplot(test$beta, betaest)

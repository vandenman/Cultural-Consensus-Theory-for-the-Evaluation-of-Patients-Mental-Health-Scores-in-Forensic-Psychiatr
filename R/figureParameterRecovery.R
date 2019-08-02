rm(list = ls())

library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
source("R/simulate.R")
source("R/utils.R")


# settings for variational inference
iter           = 3e4L # maximum number of iterations
tol_rel_obj    = 1e-3 # convergence tolerance
output_samples = 1e3L # number of samples drawn from the approximate posterior

# LTRM ----

nr <- 300 # no raters
ni <- 200 # no items
nc <- 5   # no response options

# nr <- 25  # no raters
# ni <- 300 # no items
# nc <- 5   # no response options


set.seed(1)
obsList1 <- simulate(nr, ni, nc)
datList1 <- list(x = obsList1$x, nr = nr, ni = ni, nc = nc)

if (!file.exists("StanFits/fitRecovery1.rds")) {
  LTRM <- stan_model("Stan/LTRM.stan")
  fit1 <- vb(LTRM, datList1, output_samples = output_samples, iter = iter)
  saveRDS(list(fit = fit1, obsList = obsList1), file = "StanFits/fitRecovery1.rds")
} else {
  obj <- readRDS("StanFits/fitRecovery1.rds")
  fit1 <- obj$fit
  obsList1 <- obj$obsList
}

graph1 <- makeGraph(fit1, obsList1, c("alpha", "beta", "zeta", "lt", "kappa"), nrow = 2)
saveGraph("parameterRecoveryModel1.pdf", graph1)

# Extended LTRM ----
nr <- 20
ni <- 30
nc <- 5
np <- 50
nl <- 3
ncatR <- 2
ncatP <- 5

set.seed(42)
obsList2 <- simulate(nr, ni, nc, np, nl, ncatR, ncatP)

datList2 <- list(
  nr = nr, ni = ni, nc = nc, np = np, nl = nl,
  ncatP = ncatP, ncatR = ncatR,
  x                = obsList2$x,
  latentMembership = obsList2$latentMembership, # assumed to be known
  firstItemOfLv    = match(1:nl, obsList2$latentMembership),
  designR          = obsList2$designR,
  designP          = obsList2$designP
)

if (!file.exists("StanFits/fitRecovery2.rds")) {
  extendedLTRM <- stan_model("Stan/extendedLTRM.stan")
  fit2 <- vb(extendedLTRM, datList2, output_samples = output_samples)
  saveRDS(list(fit = fit2, obsList = obsList2), file = "StanFits/fitRecovery2.rds")
} else {
  obj <- readRDS("StanFits/fitRecovery2.rds")
  fit2 <- obj$fit
  obsList2 <- obj$obsList
}

graph2 <- makeGraph(fit2, obsList2, c("lt", "eta", "kappa", "reg", "alpha", "beta", "zeta"))
saveGraph("parameterRecoveryModel2.pdf", graph2)

# mean of beta matches observed means closely
estsBeta <- getEsts(c("beta", "muBeta"), as.matrix(fit2))

df <- data.frame(
  true = rep(tapply(obsList$beta, obsList$designR, mean), 2),
  est  = c(tapply(estsBeta$beta, obsList$designR, mean), estsBeta$muBeta),
  what = rep(c("Mean of posterior mean for individual Betas", "Estimated population mean"), each = 2)
)
ggplot(data = df, aes(x = true, y = est, group = what, color = what)) +
  geom_abline() +
  geom_point()
plot(tapply(obsList$beta, obsList$designR, mean), tapply(estsBeta$beta, obsList$designR, mean),
     pch = unique(obsList$designR))
points(tapply(obsList$beta, obsList$designR, mean), estsBeta$muBeta, col = 2, pch = unique(obsList$designR))
abline(0, 1)

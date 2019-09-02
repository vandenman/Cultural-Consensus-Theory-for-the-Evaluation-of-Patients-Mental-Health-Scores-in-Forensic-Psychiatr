rm(list = ls())

library(rstan)
rstan_options(auto_write = TRUE)
library(bayesplot)
color_scheme_set("viridis")
source("R/simulate.R")


model1 <- stan_model("Stan/simpleModel1.stan")

nr <- 150
ni <- 20
nc <- 5

obsList <- simulate1(nr, ni, nc)
plot(colSums(obsList$x), obsList$lt)

datList <- list(x = obsList$x, nr = nr, ni = ni, nc = nc, muT = .5, sdT = 2)

fit <- sampling(model1, datList, chains = 2)

samples <- as.matrix(fit)[, 1:ni]

bayesplot::mcmc_trace(samples, pars = paste0("lt[", 1:4, "]"))
bayesplot::mcmc_recover_intervals(x = samples, true = obsList$lt)
qplot(obsList$lt, colMeans(samples)) + geom_abline() + theme_bw()


model2 <- stan_model("Stan/simpleModel2.stan")

nr <- 50
ni <- 50
nc <- 5

obsList <- simulate2(nr, ni, nc)
plot(colSums(obsList$x), obsList$lt)

datList <- list(x = obsList$x, nr = nr, ni = ni, nc = nc, muT = .5, sdT = 2)

fit2 <- sampling(model2, datList, chains = 2)

samples2Arr <- as.array(fit2)
samples2Mat <- as.matrix(fit2)

idxLt <- startsWith(colnames(samples2Mat), "lt[")
estLt <- colMeans(samples2Mat[, idxLt])
idxAlpha <- startsWith(colnames(samples2Mat), "alpha[")
estAlpha <- colMeans(samples2Mat[, idxAlpha])
idxBeta <- startsWith(colnames(samples2Mat), "beta[")
estBeta <- colMeans(samples2Mat[, idxBeta])

bayesplot::mcmc_trace(samples2Arr, pars = paste0("lt[", 1:4, "]"))
bayesplot::mcmc_recover_intervals(x = samples2Mat[, idxLt], true = obsList$lt)
bayesplot::mcmc_recover_intervals(x = samples2Mat[, idxAlpha], true = obsList$alpha)
bayesplot::mcmc_recover_intervals(x = samples2Mat[, idxAlpha], true = obsList$alpha)
qplot(obsList$lt, estLt) + geom_abline() + theme_bw()
qplot(obsList$alpha, estAlpha) + geom_abline() + theme_bw()
qplot(obsList$beta, estBeta) + geom_abline() + theme_bw()



model3 <- stan_model("Stan/simpleModel3.stan")

nr <- 50
ni <- 50
nc <- 5

obsList <- simulate3(nr, ni, nc)
plot(colSums(obsList$x), obsList$lt)

datList <- list(x = obsList$x, nr = nr, ni = ni, nc = nc, muT = .5, sdT = 2)

fit3 <- sampling(model3, datList, chains = 2)

samples3Arr <- as.array(fit3)
samples3Mat <- as.matrix(fit3)

idxLt <- startsWith(colnames(samples3Mat), "lt[")
estLt <- colMeans(samples3Mat[, idxLt])
idxAlpha <- startsWith(colnames(samples3Mat), "alpha[")
estAlpha <- colMeans(samples3Mat[, idxAlpha])
idxBeta <- startsWith(colnames(samples3Mat), "beta[")
estBeta <- colMeans(samples3Mat[, idxBeta])

bayesplot::mcmc_trace(samples3Arr, pars = paste0("lt[", 1:4, "]"))
bayesplot::mcmc_recover_intervals(x = samples3Mat[, idxLt], true = obsList$lt)
bayesplot::mcmc_recover_intervals(x = samples3Mat[, idxAlpha], true = obsList$alpha)
bayesplot::mcmc_recover_intervals(x = samples3Mat[, idxBeta], true = obsList$beta)
qplot(obsList$lt, estLt) + geom_abline() + theme_bw()
qplot(obsList$alpha, estAlpha) + geom_abline() + theme_bw()
qplot(obsList$beta, estBeta) + geom_abline() + theme_bw()


model4 <- stan_model("Stan/simpleModel4.stan")

nr <- 50
ni <- 40
nc <- 5

obsList <- simulate4(nr, ni, nc)
plot(colSums(obsList$x), obsList$lt)

datList <- list(x = obsList$x, nr = nr, ni = ni, nc = nc,
                muLt       = 0,
                muBeta     = 0,
                sdLt       = 1,
                shapeAlpha = .1,
                rateAlpha  = .1,
                sdBeta     = 1,
                sdError    = 1)

fit4 <- sampling(model4, datList, chains = 2)

samples4Arr <- as.array(fit4)
samples4Mat <- as.matrix(fit4)

idxLt <- startsWith(colnames(samples4Mat), "lt[")
estLt <- colMeans(samples4Mat[, idxLt])
idxAlpha <- startsWith(colnames(samples4Mat), "alpha[")
estAlpha <- colMeans(samples4Mat[, idxAlpha])
idxBeta <- startsWith(colnames(samples4Mat), "beta[")
estBeta <- colMeans(samples4Mat[, idxBeta])
idxError <- startsWith(colnames(samples4Mat), "error[")
estError <- colMeans(samples4Mat[, idxError])

bayesplot::mcmc_trace(samples4Arr, pars = c(paste0("lt[",    sample(ni, 2), "]"), 
                                            paste0("alpha[", sample(nr, 2), "]"),
                                            paste0("beta[",  sample(nr, 2), "]"),
                                            paste0("error[", sample(ni, 2), ",", sample(nr, 2), "]")))
bayesplot::mcmc_recover_intervals(x = samples4Mat[, idxLt], true = obsList$lt)
bayesplot::mcmc_recover_intervals(x = samples4Mat[, idxAlpha], true = obsList$alpha)
bayesplot::mcmc_recover_intervals(x = samples4Mat[, idxBeta], true = obsList$beta)
# bayesplot::mcmc_recover_intervals(x = samples4Mat[, idxError], true = c(obsList$error))
qplot(obsList$lt, estLt) + geom_abline() + theme_bw()
qplot(obsList$alpha, estAlpha) + geom_abline() + theme_bw()
qplot(obsList$beta, estBeta) + geom_abline() + theme_bw()
qplot(c(obsList$error), estError) + geom_abline() + theme_bw()



model5 <- stan_model("Stan/simpleModel5.stan")

nr <- 50
ni <- 40
nc <- 5
ncat <- 3

obsList <- simulate5(nr, ni, nc, ncat)
plot(colSums(obsList$x), obsList$lt)

datList <- list(x = obsList$x, nr = nr, ni = ni, nc = nc, ncat = ncat,
                raterCategory = obsList$design,
                muLt       = 0,
                muBeta     = 0,
                sdLt       = 1,
                shapeAlpha = .1,
                rateAlpha  = .1,
                sdBeta     = 1,
                sdError    = 1,
                priorEffRater = 1)

fit5 <- vb(model5, datList)

# fit5 <- sampling(model5, datList, chains = 2)

samples5Arr <- as.array(fit5)
samples5Mat <- as.matrix(fit5)

idxLt <- startsWith(colnames(samples5Mat), "lt[")
estLt <- colMeans(samples5Mat[, idxLt])
idxAlpha <- startsWith(colnames(samples5Mat), "alpha[")
estAlpha <- colMeans(samples5Mat[, idxAlpha])
idxBeta <- startsWith(colnames(samples5Mat), "beta[")
estBeta <- colMeans(samples5Mat[, idxBeta])
idxError <- startsWith(colnames(samples5Mat), "error[")
estError <- colMeans(samples5Mat[, idxError])
idxCatEffect <- startsWith(colnames(samples5Mat), "fixedRater[")
estCatEffect <- colMeans(samples5Mat[, idxCatEffect])

bayesplot::mcmc_trace(samples5Arr, pars = c(paste0("lt[",    sample(ni, 2), "]"), 
                                            paste0("alpha[", sample(nr, 2), "]"),
                                            paste0("beta[",  sample(nr, 2), "]"),
                                            paste0("error[", sample(ni, 2), ",", sample(nr, 2), "]")))
bayesplot::mcmc_recover_intervals(x = samples5Mat[, idxLt], true = obsList$lt)
bayesplot::mcmc_recover_intervals(x = samples5Mat[, idxAlpha], true = obsList$alpha)
bayesplot::mcmc_recover_intervals(x = samples5Mat[, idxBeta], true = obsList$beta)
bayesplot::mcmc_recover_intervals(x = samples5Mat[, idxCatEffect], true = obsList$catEffect)
# bayesplot::mcmc_recover_intervals(x = samples4Mat[, idxError], true = c(obsList$error))
qplot(obsList$lt, estLt) + geom_abline() + theme_bw()
qplot(obsList$alpha, estAlpha) + geom_abline() + theme_bw()
qplot(obsList$beta, estBeta) + geom_abline() + theme_bw()
qplot(c(obsList$error), estError) + geom_abline() + theme_bw()
qplot(c(obsList$catEffect), estCatEffect) + geom_abline() + theme_bw()



# plotScatter <- function()


model6 <- stan_model("Stan/simpleModel6.stan")

nr <- 30
ni <- 20
nc <- 5
np <- 100
ncat <- 7

obsList <- simulate6(nr, ni, nc, ncat, np)
plot(colSums(obsList$x[, , 1]), obsList$lt[, 1])
sapply(seq_len(np), function(p) cor(colSums(obsList$x[, , p]), obsList$lt[, p]))

datList <- list(x = obsList$x, nr = nr, ni = ni, nc = nc, ncat = ncat, np = np,
                raterCategory = obsList$design,
                muLt       = 0,
                muBeta     = 0,
                sdLt       = 1,
                shapeAlpha = .1,
                rateAlpha  = .1,
                sdBeta     = 1,
                priorEffRater = 1)

fit6 <- vb(model6, datList)

# fit5 <- sampling(model5, datList, chains = 2)

samples6Arr <- as.array(fit6)
samples6Mat <- as.matrix(fit6)


idxLt <- startsWith(colnames(samples6Mat), "lt[")
estLt <- colMeans(samples6Mat[, idxLt])
idxAlpha <- startsWith(colnames(samples6Mat), "alpha[")
estAlpha <- colMeans(samples6Mat[, idxAlpha])
idxBeta <- startsWith(colnames(samples6Mat), "beta[")
estBeta <- colMeans(samples6Mat[, idxBeta])
idxError <- startsWith(colnames(samples6Mat), "error[")
estError <- colMeans(samples6Mat[, idxError])
idxCatEffect <- startsWith(colnames(samples6Mat), "fixedRater[")
estCatEffect <- colMeans(samples6Mat[, idxCatEffect])
idxSdError <- startsWith(colnames(samples6Mat), "sdError[")
estSdError <- colMeans(samples6Mat[, idxSdError])


# bayesplot::mcmc_trace(samples6Arr, pars = c(paste0("lt[",    sample(ni, 2), "]"), 
#                                             paste0("alpha[", sample(nr, 2), "]"),
#                                             paste0("beta[",  sample(nr, 2), "]"),
#                                             paste0("error[", sample(ni, 2), ",", sample(nr, 2), "]")))
# bayesplot::mcmc_recover_intervals(x = samples6Mat[, idxLt], true = c(obsList$lt))
# bayesplot::mcmc_recover_intervals(x = samples6Mat[, idxAlpha], true = obsList$alpha)
# bayesplot::mcmc_recover_intervals(x = samples6Mat[, idxBeta], true = obsList$beta)
# bayesplot::mcmc_recover_intervals(x = samples6Mat[, idxCatEffect], true = obsList$catEffect)
# bayesplot::mcmc_recover_intervals(x = samples4Mat[, idxError], true = c(obsList$error))
gList <- list(
  qplot(c(obsList$lt), estLt) + geom_abline() + theme_bw(),
  qplot(obsList$alpha, estAlpha) + geom_abline() + theme_bw(),
  qplot(obsList$beta, estBeta) + geom_abline() + theme_bw(),
  qplot(c(obsList$error), estError) + geom_abline() + theme_bw(),
  qplot(c(obsList$catEffect), estCatEffect) + geom_abline() + theme_bw(),
  qplot(c(obsList$sdError), estSdError) + geom_abline() + theme_bw()
)
gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = gList))






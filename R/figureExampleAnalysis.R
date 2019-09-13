rm(list = ls())

library(rstan)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
library(ggplot2)
source("R/simulate.R")
source("R/utils.R")

nr <- 10
ni <- 20
nc <- 5
np <- 50
nl <- 3
ncatR <- 2
ncatP <- 5
# seed <- 7
seed <- pi

set.seed(seed)
obsList <- simulate(nr, ni, nc, np, nl, ncatR, ncatP)
datList <- obsList2datList(obsList)

# fit LTRM
if (!file.exists("StanFits/fitExampleAnalysis.rds")) {
  extendedLTRM <- stan_model("Stan/extendedLTRM.stan")
  fit <- vb(extendedLTRM, datList, output_samples = 2e4L, seed = seed)
  saveRDS(list(fit = fit, obsList = obsList), file = "StanFits/fitExampleAnalysis.rds")
} else {
  obj <- readRDS("StanFits/fitExampleAnalysis.rds")
  fit <- obj$fit
  obsList <- obj$obsList
}

# parameter recovery looks fine
graph <- makeGraph(fit, obsList, c("lt", "eta", "alpha", "beta", "kappa", "reg", "zeta"))

# names of constructs for this fictious example
constructNms <- c("Aggressiveness", "Anxiety", "Depression")

samples <- as.matrix(fit)
mat <- computeSimilarities(obsList, samples)

# do the sample means observations disagree?
graph00 <- plotPostMeanVersusObsMean(obsList, samples)
graph0  <- plotDiffPostMeanVersusDiffObsMean(mat, constructNms)

# there are some differences
head(mat[mat[, "DiffPosterior"] > 1 & mat[, "DiffObserved"] < 0.05, ], 20)

# find the patient combination where the difference between observed and posterior means is largest,
# and the difference in observed means is small.
ii <- which.max((mat[, "DiffPosterior"] - mat[, "DiffObserved"]) * (mat[, "DiffObserved"] < 0.05))
idx1 <- mat[ii, "Patient 1"]
idx2 <- mat[ii, "Patient 2"]

obs <- rbind(c(obsList$x[, , idx1]), c(obsList$x[, , idx2]))
apply(obs, 1L, function(x) which.max(table(x)))

se <- function(x) sd(x) / sqrt(length(x))
mu1 <- tapply(obs[1L, ], rep(obsList$latentMembership, each = nr), mean)
mu2 <- tapply(obs[2L, ], rep(obsList$latentMembership, each = nr), mean)
se1 <- tapply(obs[1L, ], rep(obsList$latentMembership, each = nr), se)
se2 <- tapply(obs[2L, ], rep(obsList$latentMembership, each = nr), se)
tb2show <- cbind(rbind(mu1, mu2), rowMeans(obs))
dimnames(tb2show) <- list(c("Patient 1", "Patient 2"), c(constructNms, "overall Mean"))

s1 <- getSamplesOfPatient(samples, idx1, obsList)
s2 <- getSamplesOfPatient(samples, idx2, obsList)

graph1 <- plotEtaDensity(s1$eta, s2$eta, constructNms)
obj <- plotEtaDiffDensity(s1$eta, s2$eta, constructNms)
graph2 <- obj$graph
print(tb2show)
round(rbind(se1, se2), 2)
print(round(obj$probs, 2))

saveGraph("corrObsMeanPostMean.pdf",     graph00, height = 10, width = 10)
saveGraph("diffCorrObsMeanPostMean.pdf", graph0,  height = 10, width = 10)
saveGraph("twoPatientsDensity.pdf",      graph1,  height = 7, width = 15)
saveGraph("twoPatientsDiffDensity.pdf",  graph2,  height = 7, width = 15)
xtable::xtable(tb2show)

# let's find out where the difference between patient 1 and 2 lies!
idxLM <- obsList$latentMembership
nmIdxLM <- constructNms[idxLM]
# difference in item difficulty
plot(obsList$kappa[, idx1], obsList$kappa[, idx2], col = idxLM, pch = idxLM, cex = 1.25,
     bty = "n", las = 1)
abline(0, 1)
legend("topright", legend = paste(constructNms), col = 1:3, pch = 1:3, bty = "n")

# means of item difficulty per latent construct
by(obsList$kappa[, c(idx1, idx2)], nmIdxLM, colMeans)

# means of item truth per latent construct
by(obsList$lt[, c(idx1, idx2)], nmIdxLM, colMeans)

# patients comitted different crimes
obsList$designP[c(idx1, idx2)]

# group level means of the latent constructs
t(obsList$muEta[, obsList$designP[c(idx1, idx2)]])

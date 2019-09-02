rm(list = ls())

library(rstan)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
library(bayesplot)
library(tibble)
color_scheme_set("viridis")
source("R/simulate.R")
source("R/utils.R")

model1 <- stan_model("Stan/modelFigure1.stan")
# expose_stan_functions(model4)
nr <- 1
ni <- 50
nc <- 5

obsList <- simulate(nr, ni, nc)

datList <- list(
  x = matrix(obsList$x, nr, ni), nr = nr, ni = ni, nc = nc#,
  # muLt       = 0,
  # muBeta     = 0,
  # sdLt       = 10,
  # sdBeta     = 10,
  # shapeAlpha = .1,
  # rateAlpha  = .1,
  # shapeZeta  = .1,
  # rateZeta   = .1,
  # shapeKappa = .1,
  # rateKappa  = .1
)

fit1 <- vb(model1, datList)
fit1 <- sampling(model1, datList)
samples1 <- as.matrix(fit1)

plotTrueVsEst(samples1, "lt",    obsList)
plotTrueVsEst(samples1, "alpha", obsList)
plotTrueVsEst(samples1, "beta",  obsList)
# plotTrueVsEst(samples1, "zeta",  obsList)
plotTrueVsEst(samples1, "kappa", obsList)



tb <- data.frame()
for (nm in names(ests)) {
  nm2show <- if (nm == "lt") "italic(T)" else nm
  tb2add <- tibble(
    estimated = ests[[nm]],
    true      = obsList[[nm]],
    parameter = nm2show
  )
  tb <- rbind(tb, tb2add)
}


grey  <- "grey80"
grey2 <- "grey83"
ggplot(tb, aes(x = true, y = estimated)) +
  geom_abline() + geom_point() +
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) +
  labs(x = "True value", y = "Posterior Mean") +
  theme(
    axis.ticks       = element_line(size = 1),
    panel.grid       = element_line(colour = grey),
    panel.grid.minor = element_line(colour = grey2),
    strip.background = element_rect(colour = grey, fill = "white"),
    panel.background = element_rect(fill = "transparent", colour = grey),
    text             = element_text(size = 30)
  )


model2 <- stan_model("Stan/modelFigure2.stan")
nr <- 50
ni <- 20
nc <- 5
np <- 30
nl <- 3

obsList <- simulateModel2(nr, ni, nc, np, nl)
plot(colSums(obsList$x[, , 1]), obsList$lt[, 1])
sapply(seq_len(np), function(p) cor(colSums(obsList$x[, , p]), obsList$lt[, p]))

datList <- list(
  x = obsList$x, nr = nr, ni = ni, nc = nc, np = np, nl = nl,
  latentMembership = obsList$latentMembership, # assumed to be known
  firstItemOfLv    = match(1:nl, obsList$latentMembership),
  muLt       = 0,
  muBeta     = 0,
  sdLt       = 10,
  sdBeta     = 10,
  shapeAlpha = .1,
  rateAlpha  = .1,
  shapeZeta  = .1,
  rateZeta   = .1,
  shapeKappa = .1,
  rateKappa  = .1
)

fit2 <- vb(model2, datList)
samples2 <- as.matrix(fit2)

plotTrueVsEst(samples2, "lt",    obsList)
plotTrueVsEst(samples2, "alpha", obsList)
plotTrueVsEst(samples2, "beta",  obsList)
# plotTrueVsEst(samples4Mat, "zeta",  obsList)
plotTrueVsEst(samples2, "kappa", obsList)
plotTrueVsEst(samples2, "reg", obsList)
obsList2 <- obsList
obsList2$eta <- t(obsList2$eta)
plotTrueVsEst(samples2, "eta", obsList2)




model3 <- stan_model("Stan/modelFigure3.stan")
nr <- 50
ni <- 40
nc <- 5
np <- 30
nl <- 3
ncatR <- 3
ncatP <- 5

obsList <- simulateModel3(nr, ni, nc, np, nl, ncatR, ncatP)
plot(scale(colSums(obsList$x[, , 1])), scale(obsList$lt[, 1])); abline(0, 1)
sapply(seq_len(np), function(p) cor(colSums(obsList$x[, , p]), obsList$lt[, p]))
plot(obsList$muBeta, tapply(obsList$beta, obsList$designR, mean))

datList <- list(
  x = obsList$x, nr = nr, ni = ni, nc = nc, np = np, nl = nl,
  latentMembership = obsList$latentMembership, # assumed to be known
  firstItemOfLv    = match(1:nl, obsList$latentMembership),
  designR = obsList$designR, designP = obsList$designP,
  muLt       = 0,
  muBeta     = 0,
  sdLt       = 10,
  sdBeta     = .1,
  shapeAlpha = .1,
  rateAlpha  = .1,
  shapeZeta  = .1,
  rateZeta   = .1,
  shapeKappa = .1,
  rateKappa  = .1
)

fit3 <- vb(model3, datList)
samples3 <- as.matrix(fit3)

nms <- c("lt", "alpha", "beta", "kappa", "eta", "reg", "muBeta", "sdBeta", "muEta", "sdEta", "sdLt")
ests <- getEsts(nms, samples3)

# debugonce(plotTrueVsEst)

plotTrueVsEst(samples3, "lt",     obsList)
plotTrueVsEst(samples3, "alpha",  obsList)
plotTrueVsEst(samples3, "beta",   obsList)
plotTrueVsEst(samples3, "kappa",  obsList)
plotTrueVsEst(samples3, "reg",    obsList)
plotTrueVsEst(samples3, "eta",    obsList)
plotTrueVsEst(samples3, "muBeta", obsList)
plotTrueVsEst(samples3, "muBeta", list(muBeta = tapply(obsList$beta, obsList$designR, mean)))
plot(obsList$muBeta, tapply(obsList$beta, obsList$designR, mean))
plot(ests$muBeta, tapply(obsList$beta, obsList$designR, mean)); abline(0, 1)
plotTrueVsEst(samples3, "sdBeta", obsList)
plotTrueVsEst(samples3, "muEta",  obsList)


# long format ----
library(randomForest)
library(foreach)
model3b <- stan_model("Stan/modelFigure3LongFormat.stan")
nr <- 50
ni <- 40
nc <- 5
np <- 30
nl <- 3
ncatR <- 3
ncatP <- 5

obsList <- simulate(nr, ni, nc, np, nl, ncatR, ncatP)
plot(scale(colSums(obsList$x[, , 1])), scale(obsList$lt[, 1])); abline(0, 1)
sapply(seq_len(np), function(p) cor(colSums(obsList$x[, , p]), obsList$lt[, p]))
plot(obsList$muBeta, tapply(obsList$beta, obsList$designR, mean))

datLongFormat <- reshapeSim(obsList)

nObs <- nrow(datLongFormat)
nHoldOut <- 0.1 * nObs
idxHoldOut <- sample(nObs, nHoldOut, replace = FALSE)

datTrain <- datLongFormat[-idxHoldOut, ]
datTest  <- datLongFormat[idxHoldOut,  ]

datList <- list(
  nObs = nrow(datTrain),
  x                 = datTrain[, "Outcome"],
  idxRater          = as.numeric(datTrain[, "Rater"]),
  idxItem           = as.numeric(datTrain[, "Item"]),
  idxPatient        = as.numeric(datTrain[, "Patient"]),
  nMissing          = nHoldOut,
  idxRaterMissing   = as.numeric(datTest[, "Rater"]),
  idxItemMissing    = as.numeric(datTest[, "Item"]),
  idxPatientMissing = as.numeric(datTest[, "Patient"]),

  nr = nr, ni = ni, nc = nc, np = np, nl = nl,
  latentMembership = obsList$latentMembership, # assumed to be known
  firstItemOfLv    = match(1:nl, obsList$latentMembership),
  designR = obsList$designR, designP = obsList$designP,
  muLt       = 0,
  muBeta     = 0,
  sdLt       = 10,
  sdBeta     = .1,
  shapeAlpha = .1,
  rateAlpha  = .1,
  shapeZeta  = .1,
  rateZeta   = .1,
  shapeKappa = .1,
  rateKappa  = .1
)

fit3b <- vb(model3b, datList)
samples3b <- as.matrix(fit3b)

nms <- c("lt", "alpha", "beta", "kappa", "eta", "reg", "muBeta", "sdBeta", "muEta", "sdEta", "xPred")
ests <- getEsts(nms, samples3b)

# debugonce(plotTrueVsEst)

plotTrueVsEst(samples3, "lt",     obsList)
plotTrueVsEst(samples3, "alpha",  obsList)
plotTrueVsEst(samples3, "beta",   obsList)
plotTrueVsEst(samples3, "kappa",  obsList)
plotTrueVsEst(samples3, "reg",    obsList)
plotTrueVsEst(samples3, "eta",    obsList)
plotTrueVsEst(samples3, "muBeta", obsList)
plotTrueVsEst(samples3, "muBeta", list(muBeta = tapply(obsList$beta, obsList$designR, mean)))
plot(obsList$muBeta, tapply(obsList$beta, obsList$designR, mean))
plot(ests$muBeta, tapply(obsList$beta, obsList$designR, mean)); abline(0, 1)
plotTrueVsEst(samples3, "sdBeta", obsList)
plotTrueVsEst(samples3, "muEta",  obsList)

idxNew <- which(startsWith(colnames(samples3b), "xPred["))

preds <- integer(nHoldOut)
for (i in seq_along(preds)) {
  j <- idxNew[i]
  tb <- table(c(1:5, samples3b[, j])) - rep.int(1, 5)
  preds[i] <- unname(which.max(tb)) # posterior mode
}
tb <- table(preds, datTest[, colnames(datTest) == "Outcome"])
sum(diag(tb)) / sum(tb)

# compare predictions
datTrain2 <- datTrain
datTrain2$Outcome <- factor(datTrain2$Outcome)
rf <- randomForest(Outcome ~ ., data = datTrain2)
plot(rf$err.rate[, 1])
# rf$
preds <- predict(rf, newdata = datTest[, colnames(datTest) != "Outcome"])
confusionTable <- table(preds, datTest[, colnames(datTest) == "Outcome"])
predAccuracy <- sum(diag(confusionTable)) / sum(confusionTable)

rf <- foreach(
  ntree = rep(50, 8),
  .combine = randomForest::combine,
  .multicombine=TRUE,
  .packages='randomForest') %dopar% {
    randomForest(Outcome ~ ., data = datTrain2, ntree = ntree)
  }

install.packages("beepr")
beepr::beep(5)


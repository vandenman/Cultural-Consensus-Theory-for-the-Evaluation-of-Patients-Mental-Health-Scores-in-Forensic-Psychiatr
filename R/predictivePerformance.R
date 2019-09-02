rm(list = ls())

# setup ----
library(rstan)
library(randomForest)
library(gbm)
library(caret)
library(doParallel)
rstan_options(auto_write = TRUE)
source("R/simulate.R")
source("R/utils.R")

fitRandomforest <- function(datTrain, datTest) {

  tStart <- Sys.time()
  datTrain <- as.data.frame(lapply(datTrain, factor))
  datTest  <- as.data.frame(lapply(datTest, factor))

  rfObj <- ranger::ranger(dependent.variable.name = "Outcome", data = datTrain,
                          oob.error = TRUE, verbose = TRUE, num.trees = 1e4)

  preds <- predict(rfObj, data = datTest[, colnames(datTest) != "Outcome"])
  confusionTable <- table(preds$predictions, datTest[, colnames(datTest) == "Outcome"])
  predAccuracy <- sum(diag(confusionTable)) / sum(confusionTable)
  tStop <- Sys.time()
  return(list(
    # trainAccuracy   = sum(diag(rfObj$confusion.matrix)) / sum(rfObj$confusion.matrix),
    predAccuracy    = predAccuracy,
    confusionTable  = confusionTable,
    rfObj           = rfObj,
    time            = tStop - tStart
  ))
}

fitBoosting <- function(datTrain, datTest, n.trees = 500, cv.folds = 5, shrinkage = 1e-3,
                        train = FALSE, unTrained = TRUE) {

  tStart <- Sys.time()
  for (i in 1:6) {
    datTrain[[i]] <- factor(datTrain[[i]])
    datTest[[i]]  <- factor(datTest[[i]])
  }
  # datTrain <- as.data.frame(lapply(datTrain, factor))
  # datTest  <- as.data.frame(lapply(datTest, factor))

  # code adapted from from Kuhn (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software. 5(28)
  # gbmGrid <- expand.grid(
  #   interaction.depth = c(1, 5, 10),
  #   n.trees           = seq(800, 1200, 100),
  #   shrinkage         = 1e-3,
  #   n.minobsinnode    = c(1, 10, 20)
  # )

  if (unTrained) {
    # don't attempt to optimize hyperparameters
    gbmObj <- gbm(
      formula = Outcome ~ .,
      data    = datTrain,
      distribution = "multinomial"
    )

  } else {
    # optimize hyperparameters

    gbmGrid <- expand.grid(
      interaction.depth = 1,#c(1, 5, 10),
      n.trees           = seq(800, 1200, 200),
      shrinkage         = 1e-3,
      n.minobsinnode    = c(20, 50, 100)
    )

    if (train) {
      k <- 3L
      nSeeds <- 9
      seeds <- vector("list", k + 1L)
      for (i in seq_len(k)) seeds[[i]] <- (1+(i-1)*nSeeds):(i*nSeeds)
      seeds[[k + 1L]] <- 0L

      control <- trainControl(
        method = "cv",
        seeds = seeds,
        index = createFolds(datTrain$Outcome, k = k)
      )
      cat("Traning gbm, this may take several hours!\n")
      gbmFit <- train(
        x              = datTrain[colnames(datTrain) != "Outcome"],
        y              = datTrain[["Outcome"]],
        method         = "gbm",
        trControl      = control,
        verbose        = TRUE,
        bag.fraction   = 0.5,
        train.fraction = 0.2,
        # nTrain         = 30,
        tuneGrid       = gbmGrid,
        distribution   = "multinomial"
      )
      cat("Saving trained model.\n")
      saveRDS(gbmFit, file.path("data", "gbmFit.rds"))
    } else {
      cat("Loading trained model.\n")
      file <- file.path("data", "gbmFit.rds")
      if (file.exists(file)) {
        gbmFit <- readRDS(file)
      } else {# fallback of
        cat("No Trained model exists, using falling back to reasonable defaults.\n")
        gbmFit <- list(bestTune = list(
          n.trees           = 1e3L,
          interaction.depth = 1L,
          shrinkage         = 1e-3,
          n.minobsinnode    = 50L
        ))
      }
    }
    n.trees           <- gbmFit$bestTune[, "n.trees"]
    interaction.depth <- gbmFit$bestTune[, "interaction.depth"]
    shrinkage         <- gbmFit$bestTune[, "shrinkage"]
    n.minobsinnode    <- gbmFit$bestTune[, "n.minobsinnode"]

    gbmObj <- gbm(
      formula = Outcome ~ .,
      data    = datTrain,
      distribution = "multinomial",
      n.trees = n.trees,
      shrinkage = shrinkage,
      interaction.depth = interaction.depth,
      n.minobsinnode = n.minobsinnode
    )
  }
  predsProbs <- predict(gbmObj, newdata = datTest[, colnames(datTest) != "Outcome"], n.trees = gbmObj$n.trees,
                    type = "response")
  preds <- apply(predsProbs, 1L, which.max)
  confusionTable <- table(preds, datTest[, colnames(datTest) == "Outcome"])
  predAccuracy <- sum(diag(confusionTable)) / sum(confusionTable)
  tStop <- Sys.time()
  return(list(
    # trainAccuracy   = sum(diag(rfObj$confusion.matrix)) / sum(rfObj$confusion.matrix)
    predAccuracy    = predAccuracy,
    confusionTable  = confusionTable,
    gbmObj          = gbmObj,
    time            = tStop - tStart
  ))
}

fitLTRM <- function(datTrain, datTest, obsList) {

  nr <- obsList$nr
  ni <- obsList$ni
  nc <- obsList$nc
  np <- obsList$np
  nl <- obsList$nl
  ncatR <- obsList$ncatR
  ncatP <- obsList$ncatP

  cat("Compiling Stan model...\n")
  model <- stan_model("Stan/extendedLTRM_longFormat.stan")

  tStart <- Sys.time()
  nHoldOut <- nrow(datTest)
  datList <- list(
    nObs              = nrow(datTrain),
    x                 = datTrain[, "Outcome"],
    idxRater          = as.numeric(datTrain[, "Rater"]),
    idxItem           = as.numeric(datTrain[, "Item"]),
    idxPatient        = as.numeric(datTrain[, "Patient"]),
    nMissing          = nHoldOut,
    idxRaterMissing   = as.numeric(datTest[, "Rater"]),
    idxItemMissing    = as.numeric(datTest[, "Item"]),
    idxPatientMissing = as.numeric(datTest[, "Patient"]),
    nr = nr, ni = ni, nc = nc, np = np, nl = nl, ncatR = ncatR, ncatP = ncatP,
    latentMembership = c(obsList$latentMembership), # assumed to be known
    firstItemOfLv    = if (nl == 1L) 1L else c(match(1:nl, obsList$latentMembership)),
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

  cat("Starting Variational Bayes...\n")
  fit <- vb(model, datList)
  samples <- as.matrix(fit)

  idxNew <- which(startsWith(colnames(samples), "xPred["))

  preds <- integer(nHoldOut)
  for (i in seq_along(preds)) {
    j <- idxNew[i]
    tb <- table(c(1:5, samples[, j])) # add 1:5 so table always goes from 1-5.
    preds[i] <- unname(which.max(tb)) # posterior mode (unaffected by the 1-5)
  }
  confusionTable <- table(preds, datTest[, colnames(datTest) == "Outcome"])
  predAccuracy <- sum(diag(confusionTable)) / sum(confusionTable)
  tStop <- Sys.time()
  return(list(
    predAccuracy    = predAccuracy,
    confusionTable  = confusionTable,
    stanObj         = fit,
    time            = tStop - tStart
  ))
}

fitMean <- function(datTrain, datTest) {

  tStart <- Sys.time()
  obs <- tapply(datTrain$Outcome, datTrain[colnames(datTrain) %in% c("Rater", "Item", "Patient")], mean, na.rm = TRUE)

  predsMean <- numeric(nrow(datTest))
  predsMode <- numeric(nrow(datTest))
  for (o in seq_along(predsMean)) {

    r <- datTest$Rater[o]
    i <- datTest$Item[o]
    p <- datTest$Patient[o]
    datTest

    ests <- c(
      obs[r, i,  ],
      obs[r,  , p],
      obs[,  i, p]
    )

    predsMean[o] <- mean(ests, na.rm = TRUE)
    predsMode[o] <- unname(which.max(table(ests)))

  }

  confusionTableMean <- table(round(predsMean), datTest[, colnames(datTest) == "Outcome"])
  confusionTableMode <- table(predsMode, datTest[, colnames(datTest) == "Outcome"])
  predAccuracyMean <- sum(diag(confusionTableMean)) / sum(confusionTableMean)
  predAccuracyMode <- sum(diag(confusionTableMode)) / sum(confusionTableMode)
  tStop <- Sys.time()
  # return(list(
  #   predAccuracyMean   = predAccuracyMean,
  #   predAccuracyMode   = predAccuracyMode,
  #   confusionTableMean = confusionTableMean,
  #   confusionTableMode = confusionTableMode
  # ))
  #
  return(list(
    predAccuracy   = predAccuracyMode,
    confusionTable = confusionTableMode,
    time           = tStop - tStart
  ))
}

fitAllModels <- function(datTrain, datTest, obsList) {

  # fit the machine learning models
  cat("Random Forest\n")
  resRF   <- fitRandomforest(datTrain, datTest)
  cat("Boosting\n")
  resGBM  <- fitBoosting(datTrain, datTest)

  # fit the mode estimate
  cat("Mode\n")
  resMean <- fitMean(datTrain, datTest)

  # fit the extended LTRM
  cat("Extended LTRM\n")
  resLTRM <- fitLTRM(datTrain, datTest, obsList)

  allRes <- list(
    resRF   = resRF,
    resGBM  = resGBM,
    resMean = resMean,
    resLTRM = resLTRM
  )
  return(allRes)
}

# Dense data: all rater rate all patients ----
nr <- 20
ni <- 30
nc <- 5
np <- 50
nl <- 3
ncatR <- 2
ncatP <- 3

set.seed(42)
obsList <- simulate(nr, ni, nc, np, nl, ncatP = ncatP, ncatR = ncatR)
datLongFormat <- reshapeSim(obsList)

# create a training data set and a test data set
nObs <- nrow(datLongFormat)
nHoldOut <- floor(0.2 * nObs)
idxHoldOut <- sample(nObs, nHoldOut, replace = FALSE)

datTrain <- datLongFormat[-idxHoldOut, ]
datTest  <- datLongFormat[idxHoldOut,  ]

results1 <- fitAllModels(datTrain, datTest, obsList)

# # fit the machine learning models
# resRF   <- fitRandomforest(datTrain, datTest)
# resGBM  <- fitBoosting(datTrain, datTest)
#
# # fit the mode estimate
# resMean <- fitMean(datTrain, datTest)
#
# # fit the extended LTRM
# resLTRM <- fitLTRM(datTrain, datTest, obsList)
#

sapply(results1, `[`, "time")
sapply(results1, `[`, "predAccuracy")
sapply(results1, `[`, "confusionTable")

tb <- do.call(rbind, sapply(results1, `[`, "predAccuracy"))
tb <- cbind.data.frame(c("LTRM", "Random Forest", "Boosting", "Sample Mode"), tb[c(4, 1, 2, 3)])
colnames(tb) <- c("Method", "Prediction Accuracy")
print(xtable::xtable(tb[order(tb[, 2], decreasing = TRUE), ]), include.rownames = FALSE)

# Sparse data: all rater rate 10 patients ----
# the first simulation is somewhat biased because many raters rate almost all patients.
# here we simulate a more sparse design.
nr <- 20
ni <- 30
nc <- 5
np <- 50
nl <- 3
ncatR <- 2
ncatP <- 3

set.seed(31415)
obsList <- simulate(nr, ni, nc, np, nl, ncatP = ncatP, ncatR = ncatR)
datLongFormat <- reshapeSim(obsList)

# create a training data set and a test data set
# each rater only rates 10 patients

idxRated <- matrix(NA_integer_, 10, nr)
for (r in seq_len(nr)) {
  idxRated[, r] <- seq_len(10) + 4L * (r - 1L)
  ii <- idxRated[, r] > np
  idxRated[ii, r] <- (idxRated[ii, r] %% np)
}
# idxRated <- (idxRated %% 50L) + 1L # mod 50 to wrap 51 to 0 and + 1 to make 0 into 1
table(idxRated)

idx <- integer()
for (r in seq_len(nr)) {
  ii <- which(datLongFormat$Rater == r & datLongFormat$Patient %in% idxRated[, r])
  idx <- c(idx, ii)
}
datLongFormatKept <- datLongFormat[idx, ]
nObs <- nrow(datLongFormatKept)

nHoldOut <- floor(0.2 * nObs)
idxHoldOut <- sample(nObs, nHoldOut, replace = FALSE)

datTrain <- datLongFormatKept[-idxHoldOut, ]
datTest  <- datLongFormatKept[idxHoldOut,  ]

results2 <- fitAllModels(datTrain, datTest, obsList)

sapply(results2, `[`, "time")
sapply(results2, `[`, "predAccuracy")
sapply(results2, `[`, "confusionTable")

tb <- do.call(rbind, sapply(results2, `[`, "predAccuracy"))
tb <- cbind.data.frame(c("LTRM", "Random Forest", "Boosting", "Sample Mode"), tb[c(4, 1, 2, 3)])
colnames(tb) <- c("Method", "Prediction Accuracy")
print(xtable::xtable(tb[order(tb[, 2], decreasing = TRUE), ]), include.rownames = FALSE)


# Informed ML: MAP estimates as features for ML algorithms ----
# setup from dense data
nr <- 20
ni <- 30
nc <- 5
np <- 50
nl <- 3
ncatR <- 2
ncatP <- 3

set.seed(2718282)
obsList <- simulate(nr, ni, nc, np, nl, ncatP = ncatP, ncatR = ncatR)
datLongFormat <- reshapeSim(obsList)

# create a training data set and a test data set
nObs <- nrow(datLongFormat)
nHoldOut <- floor(0.2 * nObs)
idxHoldOut <- sample(nObs, nHoldOut, replace = FALSE)

datTrain <- datLongFormat[-idxHoldOut, ]
datTest  <- datLongFormat[idxHoldOut,  ]

resLTRM <- fitLTRM(datTrain, datTest, obsList)

# augment the training dataset with posterior means
samples <- as.matrix(resLTRM$stanObj)
postMeans <- getEsts(c("lt", "kappa", "alpha", "beta", "zeta", "eta", "reg"), samples)
dim(postMeans$lt)    <- c(ni, np)
dim(postMeans$kappa) <- c(ni, np)
dim(postMeans$eta)   <- c(nl, np)

for (nm in c("alpha", "beta", "zeta")) {
  datTrain[[nm]] <- postMeans[[nm]][datTrain$Rater]
  datTest[[nm]]  <- postMeans[[nm]][datTest$Rater]
}

datTrain[["lt"]]    <- postMeans[["lt"]]   [cbind(datTrain$Item, datTrain$Patient)]
datTrain[["kappa"]] <- postMeans[["kappa"]][cbind(datTrain$Item, datTrain$Patient)]
datTrain[["eta"]]   <- postMeans[["eta"]]  [cbind(obsList$latentMembership[datTrain$Item], datTrain$Patient)]
datTrain[["reg"]]   <- postMeans[["reg"]]  [datTrain$Item]

datTest[["lt"]]    <- postMeans[["lt"]]   [cbind(datTest$Item, datTest$Patient)]
datTest[["kappa"]] <- postMeans[["kappa"]][cbind(datTest$Item, datTest$Patient)]
datTest[["eta"]]   <- postMeans[["eta"]]  [cbind(obsList$latentMembership[datTest$Item], datTest$Patient)]
datTest[["reg"]]   <- postMeans[["reg"]]  [datTest$Item]

# fit the machine learning models
idx <- colnames(datTrain)[1:6]
resRF_naive    <- fitRandomforest(datTrain[idx], datTest[idx])
resRF_informed <- fitRandomforest(datTrain,      datTest)

resGBM_naive    <- fitBoosting(datTrain[idx], datTest[idx])
resGBM_informed <- fitBoosting(datTrain,      datTest)

tb <- matrix(c(
  resRF_naive$predAccuracy,  resRF_informed$predAccuracy,
  resGBM_naive$predAccuracy, resGBM_informed$predAccuracy
), 2, 2)
rownames(tb) <- c("Naive", "Informed")
colnames(tb) <- c("Random Forest", "Boosting")
print(xtable::xtable(tb))

# old ----
aa <- as.numeric(as.matrix(datTrain[, -4]))
dim(aa) <- dim(datTrain) - c(0, 1)
obj <- xgboost::xgboost(data = aa, label = datTrain[, 4]-1L, nrounds = 2000,
                        objective = "multi:softprob", num_class = 5, print_every_n = 50)
bb <- as.numeric(as.matrix(datTest[, -4]))
dim(bb) <- dim(datTest) - c(0, 1)
preds0 <- predict(obj, bb)
preds <- apply(matrix(preds0, nrow(datTest), 5L, TRUE), 1L, which.max)
confusionTable <- table(preds, datTest[, colnames(datTest) == "Outcome"])
predAccuracy <- sum(diag(confusionTable)) / sum(confusionTable)

getIdx <- function(nm, nms) {
  # idx <- startsWith(nms, paste0(nm, "["))
  # if (!any(idx)) # scalar
  idx <- startsWith(nms, nm)
  return(idx)
}

getEsts <- function(nms, samples) {

  cnms <- colnames(samples)
  output <- vector("list", length(nms))
  names(output) <- nms

  for (nm in nms) {
    idx <- getIdx(nm, cnms)
    output[[nm]] <- colMeans(samples[, idx, drop = FALSE])
  }
  return(output)
}

plotTrueVsEst <- function(samples, parEst, trueVals, parTrue = parEst, bayesplot = FALSE, demean = FALSE) {

  idx <- getIdx(parEst, colnames(samples))
  est <- c(colMeans(samples[, idx, drop = FALSE]))
  tru <- c(trueVals[[parTrue]])
  if (demean) {
    est <- est - mean(est)
    tru <- tru - mean(tru)
  }
  if (bayesplot || length(est) == 1) {
   return(bayesplot::mcmc_recover_intervals(x = samples[, idx, drop = FALSE], true = tru))
  } else {
  return(
    qplot(tru, est) +
      geom_abline() +
      labs(x = paste("True:", parTrue), y = paste("Estimated:", parTrue)) +
      ggtitle(paste("cor:", round(cor(tru, est), 3))) +
      theme_bw()
  )
  }
}

getDummy <- function(tb) {
  ranges <- sapply(split(tb, tb$parameter), function(x) range(pretty(c(x$true, x$estimated))))
  dummy <- data.frame(
    true = c(ranges),
    estimated = c(ranges),
    parameter = rep(colnames(ranges), each = 2)
  )
  return(dummy)
}

makeGraph <- function(stanObj, obsList, nms = c("lt", "alpha", "beta", "kappa", "eta", "reg", "zeta",
                                       "muEta", "muBeta", "sdBeta"),
                      nrow = 2) {

  nms <- match.arg(nms, several.ok = TRUE)
  samples <- as.matrix(stanObj)
  ests <- getEsts(nms, samples)

  tb <- data.frame()
  for (nm in names(ests)) {
    nm2show <- switch(
      nm,
      "lt"      = "theta", #"italic(T)",
      "reg"     = "lambda",
      "muEta"   = "mu[eta]",
      "muBeta"  = "mu[beta]",
      "sdBeta"  = "sigma[beta]",
      nm
    )
    tb2add <- data.frame(
      estimated = c(ests[[nm]]),
      true      = c(obsList[[nm]]),
      parameter = nm2show
    )
    tb <- rbind(tb, tb2add)
  }

  dummy <- getDummy(tb)

  grey  <- "grey80"
  grey2 <- "grey83"
  graph <- ggplot(tb, aes(x = true, y = estimated)) +
    geom_abline() + geom_point(shape = 21, fill = "grey70", stroke = 1.25, size = 2) +
    facet_wrap(~parameter, scales = "free", labeller = label_parsed, nrow = nrow) +
    geom_blank(data = dummy) +
    labs(x = "True value", y = "Posterior Mean") +
    theme(
      axis.ticks       = element_line(size = 1),
      panel.grid       = element_line(colour = grey),
      panel.grid.minor = element_line(colour = grey2),
      strip.text       = element_text(size = 40),
      strip.background = element_rect(colour = grey, fill = "white"),
      panel.background = element_rect(fill = "transparent", colour = grey),
      text             = element_text(size = 30)
    )

  return(graph)

}

saveGraph <- function(filename, graph, width = 20, height = 10) {

  p1 <- file.path("figures", filename)
  p2 <- file.path("paper", p1)

  pdf(p1, width = width, height = height)
  print(graph)
  dev.off()
  file.copy(p1, p2, overwrite = TRUE)

}


getSamplesOfPatient <- function(samples, idxPatient, obsList) {

  nr    <- obsList$nr
  ni    <- obsList$ni
  nc    <- obsList$nc
  np    <- obsList$np
  nl    <- obsList$nl
  ncatP <- obsList$ncatP

  idxGroup <- obsList$designP[idxPatient]

  baseNames <- c("lt", "eta", "muEta")

  idxLt    <- paste0("lt[",    seq_len(ni), ",", idxPatient, "]")
  idxEta   <- paste0("eta[",   seq_len(nl), ",", idxPatient, "]")
  idxMuEta <- paste0("muEta[", seq_len(nl), ",", idxGroup,   "]")

  return(list(
    lt    = samples[, idxLt,    drop = FALSE],
    eta   = samples[, idxEta,   drop = FALSE],
    muEta = samples[, idxMuEta, drop = FALSE]
  ))
}

plotEtaDensity <- function(e1, e2, constructs = NULL) {

  np <- 2L
  nl <- ncol(e1)
  nSamples <- nrow(e2)

  if (is.null(constructs))
    constructs <- paste("eta[", 1:ncol(nl), "]")

  df <- data.frame(
    value = c(e1, e2),
    Patient = factor(rep(1:np, each = nSamples * nl)),
    construct = factor(rep(constructs, each = nSamples), levels = constructs)
  )
  mapping <- aes(x = value, group = Patient, fill = Patient)

  grey  <- "grey80"
  grey2 <- "grey83"
  graph <- ggplot(data = df, mapping) +
    geom_density(alpha = 0.6, size = 1.05) +
    facet_wrap(~construct, scales = "free", nrow = 1, labeller = label_parsed) +
    viridis::scale_fill_viridis(option = "C", discrete = TRUE) +
    labs(x = "Value of Latent Construct", y = "Posterior Density") +
    theme(
      axis.ticks       = element_line(size = 1),
      panel.grid       = element_line(colour = grey),
      panel.grid.minor = element_line(colour = grey2),
      strip.background = element_rect(colour = grey, fill = "white"),
      panel.background = element_rect(fill = "transparent", colour = grey),
      text             = element_text(size = 30)
    )
  return(graph)
}

plotEtaDiffDensity <- function(e1, e2, constructs = NULL) {

  nl <- ncol(e1)
  nDens <- 2^11

  if (is.null(constructs))
    constructs <- rep(paste("eta[", 1:ncol(e1), "]"), nDens)

  d <- e1 - e2
  dens <- apply(d, 2, density, n = nDens)
  df <- do.call(rbind, lapply(dens, function(x) data.frame(x = x$x, y = x$y)))
  # df$construct <- factor(rep(paste("eta[", 1:nl, "]"), each = nDens))
  df$construct <- factor(rep(constructs, each = nDens), levels = constructs)

  xMax <- tapply(df$y, df$construct, which.max)
  yMax <- diag(do.call(rbind, tapply(df$x, df$construct, `[`, xMax)))

  probs <- numeric(length(yMax))
  df2 <- data.frame()
  for (i in seq_along(yMax)) {
    compare <- if (yMax[i] > 0) `>=` else `<`
    dftmp <- as.data.frame(dens[[i]][c("x", "y")])
    # dftmp$construct <- factor(rep(paste("eta[", i, "]"), each = nDens))
    dftmp$construct <- factor(rep(constructs[i], nDens), levels = constructs[i])
    idx <- compare(dftmp$x, 0)
    df2 <- rbind(df2, dftmp[idx, ])
    probs[i] <- mean(compare(d[, i], 0))
  }

  mapping <- aes(x = x, y = y)
  greyLine <- "grey20"
  grey  <- "grey80"
  grey2 <- "grey83"
  graph <- ggplot(data = df, mapping) +
    geom_ribbon(data = df2, aes(ymin = 0, ymax = y), fill = grey) +
    geom_line(alpha = 0.6, size = 1.05, color = greyLine) +
    facet_wrap(~construct, scales = "free", nrow = 1, labeller = label_parsed) +
    viridis::scale_fill_viridis(option = "C", discrete = TRUE) +
    labs(x = "Value of Latent Construct", y = "Posterior Density") +
    theme(
      axis.ticks       = element_line(size = 1),
      panel.grid       = element_line(colour = grey),
      panel.grid.minor = element_line(colour = grey2),
      strip.background = element_rect(colour = grey, fill = "white"),
      panel.background = element_rect(fill = "transparent", colour = grey),
      text             = element_text(size = 30)
    )
  return(list(graph = graph, probs = probs))
}

computeSimilarities <- function(obsList, samples) {

  np <- obsList$np
  nl <- obsList$nl
  latentMembership <- obsList$latentMembership
  uniq <- unique(latentMembership)

  # get all Eta samples of all patients
  idx <- apply(expand.grid(1:nl, 1:np), 1L, function(x) paste0("eta[", x[1], ",", x[2], "]"))
  sAll <- samples[, idx]
  # compute posterior means
  cnms <- matrix(colMeans(sAll), nl, np)

  mat <- matrix(NA, np * (np - 1L) / 2L * nl, 5L, dimnames = list(NULL, c("Patient 1", "Patient 2", "Eta", "DiffPosterior", "DiffObserved")))
  ci <- 1

  # compute for each latent variable and for all combinations of patients:
  #  - the difference in posterior means.
  #  - the difference in observed means.
  for (l in seq_len(nl)) {
    idx <- latentMembership == uniq[l]
    for (i in seq_len(np - 1L)) {
      for (j in (i + 1L):np) {

        mat[ci, ] <- c(i = i, j = j, eta = l,
          abs(cnms[l, i] - cnms[l, j]),
          abs(mean(obsList$x[, idx, i]) - mean(obsList$x[, idx, j]))
        )
        ci <- ci + 1L
      }
    }
  }
  return(mat)
}

plotDiffPostMeanVersusDiffObsMean <- function(mat, constructs = NULL) {
  greyLine <- "grey20"
  grey  <- "grey80"
  grey2 <- "grey83"

  if (is.null(constructs)) {
    colorscale <- viridis::scale_color_viridis(option = "C", discrete = TRUE)
    shapescale <- ggplot2::scale_shape()
    scaleName <- "Construct"
  } else {
    colorscale <- viridis::scale_color_viridis(option = "C", discrete = TRUE, labels = constructs)
    shapescale <- ggplot2::scale_shape(labels = constructs)
    scaleName <- ""
  }

  return(
  ggplot(data = as.data.frame(mat), aes(x = DiffObserved, y = DiffPosterior,
                                        group = factor(Eta), shape = factor(Eta), color = factor(Eta))) +
    geom_abline(color = greyLine) +
    geom_point(alpha = .6, size = 3) +
    colorscale +
    shapescale +
    labs(x = "Difference in Observed Means", y = "Difference in Posterior Means",
         shape = scaleName, color = scaleName) +
    theme(
        axis.ticks       = element_line(size = 1),
        panel.grid       = element_line(colour = grey),
        panel.grid.minor = element_line(colour = grey2),
        panel.background = element_rect(fill = "transparent", colour = grey),
        text             = element_text(size = 30)
    )
  )
}

plotPostMeanVersusObsMean <- function(obsList, samples) {

  np <- obsList$np
  nl <- obsList$nl

  # get all Eta samples of all patients
  idx <- apply(expand.grid(1:nl, 1:np), 1L, function(x) paste0("eta[", x[1], ",", x[2], "]"))
  sAll <- samples[, idx]
  # compute posterior means
  cmns <- matrix(colMeans(sAll), nl, np)

  tmp <- colMeans(obsList$x) # average over all raters
  uniq <- unique(obsList$latentMembership)
  cmnsObs <- 0*cmns
  for (u in uniq) {
    idx <- u == obsList$latentMembership
    cmnsObs[u, ] <- colMeans(tmp[idx, ])
  }

  dat <- data.frame(
    x = c(cmnsObs),
    y = c(cmns),
    eta = factor(rep(seq_len(nl), np))
  )

  greyLine <- "grey20"
  grey  <- "grey80"
  grey2 <- "grey83"

  return(
  ggplot(data = dat, aes(x = x, y = y, group = eta, shape = eta, color = eta)) +
    geom_abline(intercept = -mean(cmnsObs), color = greyLine) +
    geom_point(alpha = .6, size = 3, show.legend = FALSE) +
    viridis::scale_color_viridis(option = "C", discrete = TRUE) +
    labs(x = "Observed Mean", y = "Posterior Mean", shape = "Construct", color = "Construct") +
    theme(
        axis.ticks       = element_line(size = 1),
        panel.grid       = element_line(colour = grey),
        panel.grid.minor = element_line(colour = grey2),
        panel.background = element_rect(fill = "transparent", colour = grey),
        text             = element_text(size = 30)
    )
  )

}

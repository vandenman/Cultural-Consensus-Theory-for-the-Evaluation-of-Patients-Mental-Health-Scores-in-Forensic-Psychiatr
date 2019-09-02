plotOrderedDensity <- function(a = 1, b = 0, ncat = 5, gamma = NULL, xmin = -7, xmax = 7, ymax = .3, nPoints = 1e3) {

  require("ggplot2")
  require("tibble")

  # inside panel specs
  x1 <- xmin
  x2 <- xmin + 5 / 16 * (xmax - xmin)
  y1 <- ymax * 1/4
  y2 <- y1   + 9 / 16 * (ymax - 0)

  dat <- tibble(
    x = seq(xmin, xmax, length.out = nPoints),
    y = dlogis(x)
  )

  if (is.null(gamma)) {
    in01 <- seq(0, 1, length.out = ncat + 1)[-c(1, ncat + 1)]
    gamma <- qlogis(in01)
  } else {
    ncat <- length(gamma) + 1L
  }
  delta <- a*gamma + b

  datThresholds <- tibble(
    x = rep(delta, 2),
    # y = c(rep(0, ncat + 1), dlogis(unbiased)),
    y = rep(c(0, .95*ymax), each = length(delta)),
    g = rep(seq_along(delta), 2)
  )

  unbiased2 <- c(xmin, delta, xmax)
  datTxt1 <- tibble(
    x = (unbiased2[-(ncat + 1)] + unbiased2[-1]) / 2,
    y = .9 * ymax,
    label = as.character(1:ncat)
  )

  datTxt2 <- tibble(
    x = delta,
    y = ymax,
    label = paste0("delta[", seq_along(delta), "]")
  )

  datBars1 <- tibble(
    x = 1:ncat,
    y = dOrdered(seq_len(ncat), delta, 0, 1, log = FALSE)
  )

  lineCol <- "gray20"
  fillCol <- "grey80"

  labSize <- 12
  thm <- theme(
    axis.line        = element_blank(),
    # axis.title.x = element_text(size = 25),
    # axis.title.y     = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    text             = element_text(size = 30)
  )

  graph1 <- ggplot(data = dat, aes(x = x, y = y)) +
    geom_hline(yintercept = 0) +
    geom_ribbon(aes(x = x, ymin = 0, ymax = y), fill = fillCol) +
    geom_line(color = lineCol) +
    geom_line(data = datThresholds, aes(x = x, y = y, group = g), inherit.aes = FALSE) +
    geom_text(data = datTxt1, aes(x = x, y = y, label = label), size = labSize, inherit.aes = FALSE) +
    geom_text(data = datTxt2, aes(x = x, y = y, label = label), size = labSize, inherit.aes = FALSE, parse = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    thm + theme(axis.title = element_blank(), plot.margin = margin(r = 2))


  bar1 <- ggplot(data = datBars1, aes(x = x, y = y)) +
    geom_col() +
    scale_y_continuous(breaks = c(0, .2, .4), limits = c(0, .4)) +
    labs(x = "x", y = "P(X = x)") + thm +
    theme(axis.ticks = element_line(),
          axis.text  = element_text())

  combi1 <- graph1 + annotation_custom(
    grob = ggplotGrob(bar1),
    xmin = x1, xmax = x2, ymin = y1, ymax = y2
  )
  return(combi1)

}

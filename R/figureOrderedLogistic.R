rm(list = ls())
library(ggplot2)
library(tibble)
source("R/simulate.R")
source("R/utils.R")
source("R/plotOrderedDensity.R")

graph1 <- plotOrderedDensity(a = 1, b = 0)
graph2 <- plotOrderedDensity(a = 2, b = 0.5)

saveGraph("orderedLogisticUnbiased.pdf", graph1, 10, 10)
saveGraph("orderedLogisticBiased.pdf",   graph2, 10, 10)


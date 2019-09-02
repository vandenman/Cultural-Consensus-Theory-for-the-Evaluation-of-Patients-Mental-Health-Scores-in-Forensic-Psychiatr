rm(list = ls())
library(ggplot2)
library(tibble)
library(PearsonDS)
# library(ggridges)

set.seed(42)
nt <- 7
nprec <- 1e4
nn <- 1e3


l <- rnorm(nt, 1:nt / 2, 1)
s <- rep(.15, nt) #abs(rnorm(nt, .2, .5))
m <- 1+abs(rnorm(nt, 10, .5))
nu <- abs(rnorm(nt, 10, .5))

tb <- tibble()
for (t in seq_len(nt)) {
  xtmp <- seq(l[t] - 3*s[t], l[t] + 3*s[t], length.out = nn)
  fx <- dnorm(xtmp, l[t], s[t])
  fx <- .5 * fx / max(fx)
  # fx <- c(t + fx, t - fx )
  tmptb <- tibble(
    x = c(xtmp, xtmp),
    y = c(t + fx, t - fx),
    t = rep(c(t, t + 0.5), each = nn)
    # dd = rpearsonIV(nprec, m = m[t], nu = nu[t], location = l[t], scale = s[t])
  )
  # tmptb <- tmptb[order(tmptb$y), ]
  tb <- rbind(tb, tmptb)
}

thm <- theme(
  # axis.line        = element_blank(), 
  # axis.title.x = element_text(size = 25),
  # axis.title.y     = element_blank(),
  # axis.text        = element_blank(),
  axis.ticks       = element_line(size = 1),
  panel.grid       = element_blank(),
  panel.background = element_blank(),
  text             = element_text(size = 30)
)

br <- pretty(tb$x)

graph <- ggplot(tb, aes(x = x, y = y, group = t)) +
  geom_line() +
  # geom_violin() + 
  # geom_boxplot() + 
  labs(y = "Assessment", x = "Latent truth") +
  scale_y_continuous(breaks = 1:nt) + 
  geom_segment(y = 1, yend = nt, x = -Inf, xend = -Inf) +
  geom_segment(y = -Inf, yend = -Inf, x = br[1], xend = br[length(br)]) + 
  scale_x_continuous(breaks = br, limits = range(br)) +
  coord_flip() +
  thm
graph

# # tb %>% 
# ggplot(data = as.data.frame(tb), aes(x = dd, y = factor(tt))) +
#   geom_density_ridges(rel_min_height = 0.005) + 
#   coord_flip() +
#   theme_ridges()
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   geom_density_ridges(rel_min_height = 0.005)


fname <- "progressMonitoring.pdf"
p1 <- file.path("figures", fname)
p2 <- file.path("paper", p1)

pdf(p1, width = 10, height = 10)
print(graph)
dev.off()
file.copy(p1, p2, overwrite = TRUE)
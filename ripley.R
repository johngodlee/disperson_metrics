# Ripley's K and L function visualisations
# John Godlee (johngodlee@gmail.com)
# 2021-11-14

# Packages
library(dplyr)
library(ggplot2)
library(spatstat)
library(sf)
library(nngeo)
library(parallel)

# Import data
s <- read.csv("dat/trees.csv")
dat_even_csr <- readRDS("dat/dat_even_csr_points.rds")
dat_csr_clust <- readRDS("dat/dat_csr_clust_points.rds")

s_fil <- s[s$plot_id == 1,]

# spatstat Ripley's K and L
s_ppp <- ppp(s_fil$x_grid, s_fil$y_grid, c(0,100), c(0,100))

s_k <- Kest(s_ppp)
s_kenv <- envelope(s_ppp,Kest)

s_l <- Lest(s_ppp)
s_lenv <- envelope(s_ppp,Lest)

pdf(file = "img/spatstat_ripley_k.pdf", width = 10, height = 5)
par(mfrow = c(1,2)) 
plot(s_k, main = NULL)
plot(s_kenv, main = NULL)
dev.off()

pdf(file = "img/spatstat_ripley_l.pdf", width = 10, height = 5)
par(mfrow = c(1,2)) 
plot(s_l, main = NULL)
plot(s_lenv, main = NULL)
dev.off()

pdf(file = "img/spatstat_ripley_l_norm.pdf", width = 10, height = 5)
par(mfrow = c(1,2)) 
plot(s_l, . -r ~ r, ylab=expression(hat("L")), xlab = "d (m)", main = NULL)
plot(s_lenv, . -r ~ r, ylab=expression(hat("L")), xlab = "d (m)", legend = FALSE, main = NULL)
dev.off()

s_lenv$obs_r <- s_lenv$obs - s_lenv$r
s_lenv$theo_r <- s_lenv$theo - s_lenv$r
s_lenv$lo_r <- s_lenv$lo - s_lenv$r
s_lenv$hi_r <- s_lenv$hi - s_lenv$r

pdf(file = "img/spatstat_ripley_l_norm_ggplot.pdf", width = 8, height = 5)
ggplot() + 
  geom_ribbon(data = s_lenv, aes(x = r, ymin = lo_r, ymax = hi_r), alpha = 0.5) + 
  geom_line(data = s_lenv, aes(x = r, y = theo_r), linetype = 2, colour = "red") + 
  geom_line(data = s_lenv, aes(x = r, y = obs_r)) + 
  theme_bw() + 
  labs(y = expression(hat("L")), x = "d (m)")
dev.off()

adj_sample_even_csr <- seq(0, length(dat_even_csr[[1]]), 50)

dat_even_csr_l <- do.call(rbind, lapply(dat_even_csr[[1]][adj_sample_even_csr], function(x) {
  x_coords <- as.data.frame(st_coordinates(x))
  x_ppp <- ppp(x_coords$X, x_coords$Y, c(0,100), c(0,100))
  x_l <- envelope(x_ppp, Lest, verbose = FALSE)

  x_l$obs_r <- x_l$obs - x_l$r
  x_l$theo_r <- x_l$theo - x_l$r
  x_l$lo_r <- x_l$lo - x_l$r
  x_l$hi_r <- x_l$hi - x_l$r

  x_l_df <- as.data.frame(x_l)
  x_l_df$adj <- unique(x$adj)

  return(x_l_df)
}))

pdf(file = "img/even_csr_ripley_l.pdf", width = 6, height = 4)
ggplot() + 
  geom_line(data = dat_even_csr_l, 
    aes(x = r, y = obs, group = adj,  colour = adj)) + 
  scale_colour_scico(name = "Randomness", palette = "batlow") + 
  theme_bw() + 
  theme(
    legend.position = c(0.8,0.3), 
    legend.box.background = element_rect(colour = "black")) + 
  labs(y = expression(hat("L")), x = "d (m)")
dev.off()

adj_sample_csr_clust <- seq(0, length(dat_csr_clust[[1]]), 50)

dat_csr_clust_l <- do.call(rbind, lapply(dat_csr_clust[[1]][adj_sample_csr_clust], function(x) {
  x_coords <- as.data.frame(st_coordinates(x))
  x_ppp <- ppp(x_coords$X, x_coords$Y, c(0,100), c(0,100))
  x_l <- envelope(x_ppp, Lest, verbose = FALSE)

  x_l$obs_r <- x_l$obs - x_l$r
  x_l$theo_r <- x_l$theo - x_l$r
  x_l$lo_r <- x_l$lo - x_l$r
  x_l$hi_r <- x_l$hi - x_l$r

  x_l_df <- as.data.frame(x_l)
  x_l_df$adj <- unique(x$adj)

  return(x_l_df)
}))

pdf(file = "img/csr_clust_ripley_l.pdf", width = 6, height = 4)
ggplot() + 
  geom_line(data = dat_csr_clust_l, 
    aes(x = r, y = obs, group = adj,  colour = adj)) + 
  scale_colour_scico(name = "Clustering", palette = "batlow") + 
  theme_bw() + 
  theme(
    legend.position = c(0.8,0.3), 
    legend.box.background = element_rect(colour = "black")) + 
  labs(y = expression(hat("L")), x = "d (m)")
dev.off()


